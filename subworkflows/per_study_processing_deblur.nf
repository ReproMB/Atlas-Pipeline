/*
 * Subworkflow: per_study_processing_deblur
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Deblur variant of per_study_processing, for "false FASTQ" studies whose
 * quality scores are synthetic (flat Q30) and therefore unusable by DADA2.
 *
 * Identical front end (import → cutadapt → quality summary) and identical back
 * end (taxonomy → filter → collapse); only the denoise core changes:
 *
 *   DADA2 (quality-driven, learned error model, adaptive retry)
 *     ─────────────────────────────────────────────────────────►  REMOVED
 *   Deblur (static positional error model, fixed geometry-driven trim)
 *     ─────────────────────────────────────────────────────────►  IN
 *
 * The DADA2 QC / adaptive-retry machinery is intentionally dropped: it reads
 * DADA2's denoising-stats schema (filtered/merged/non-chimeric %) which Deblur
 * does not produce, and Deblur does not crash the way DADA2 does on bad quality,
 * so there is nothing to self-heal. The guardrail here is the min-reads
 * abundance gate, not a retry loop.
 *
 * DEMUX_SUMMARIZE is kept but its role is diagnostic only: it lets you eyeball
 * the flat-Q30 profile (confirming the study really is a false-FASTQ case) and
 * the read-length distribution (which is what you set the Deblur trim-length
 * from). Its JSON is passed to the denoise processes purely so overrides ride
 * along; the trim-length itself is resolved from geometry, never from quality.
 */

include { QIIME_IMPORT_PAIRED     } from '../modules/qiime_import/main'
include { QIIME_IMPORT_SINGLE     } from '../modules/qiime_import/main'
include { CUTADAPT_TRIM_PAIRED    } from '../modules/cutadapt/main'
include { CUTADAPT_TRIM_SINGLE    } from '../modules/cutadapt/main'
include { DEMUX_SUMMARIZE         } from '../modules/quality_summary/main'
include { DEBLUR_MERGE_AND_DENOISE } from '../modules/deblur/main'
include { DEBLUR_DENOISE_FWD       } from '../modules/deblur/main'
include { DEBLUR_VISUALIZE_STATS   } from '../modules/deblur/main'
include { CLASSIFY_TAXONOMY       } from '../modules/taxonomy/main'
include { FILTER_BACTERIA         } from '../modules/taxonomy/main'
include { COLLAPSE_TAXONOMY       } from '../modules/taxonomy/main'


// Decide join (merge overlapping pairs) vs forward-only, per study.
//   * SINGLE                              → always forward (nothing to join)
//   * per-study override (meta / sheet)   → wins
//   * region in params.deblur_forward_regions (long, non-overlapping at typical
//     read lengths, e.g. V3-V5)           → forward
//   * otherwise                           → join
// Match whatever the DADA2 studies did for the same region so the length regime
// is comparable at genus collapse.
def wantsJoin = { meta ->
    if (meta.paired_single != 'PAIRED') return false
    def ov = meta.denoise_override?.deblur_mode ?: meta.deblur_mode
    if (ov == 'join')    return true
    if (ov == 'forward') return false
    def fwdRegions = (params.deblur_forward_regions ?: []) as List
    return !(meta.region in fwdRegions)
}


workflow PER_STUDY_PROCESSING_DEBLUR {

    take:
    studies_ch         // channel: [meta, manifest]
    classifiers_ch     // channel: [primer_meta, classifier.qza]

    main:

    // ── 1. Import ───────────────────────────────────────────────────────
    studies_ch
        .branch { meta, manifest ->
            paired: meta.paired_single == 'PAIRED'
            single: meta.paired_single == 'SINGLE'
        }
        .set { import_branch }

    QIIME_IMPORT_PAIRED(import_branch.paired)
    QIIME_IMPORT_SINGLE(import_branch.single)

    demux_ch = QIIME_IMPORT_PAIRED.out.demux
        .mix(QIIME_IMPORT_SINGLE.out.demux)

    // ── 1b. Primer / adapter removal (cutadapt, sequence-based) ─────────
    // Deblur's trim-length truncates from the 3' end only; it does NOT strip a
    // 5' primer. Primers must therefore already be gone before denoise — either
    // the submission was pre-primer-trimmed (run_primer_trim=false) or cutadapt
    // removes them here (run_primer_trim=true). cutadapt keys off sequence, not
    // quality, so it is unaffected by the synthetic Q30.
    if (params.run_primer_trim) {
        demux_ch
            .branch { meta, demux ->
                paired: meta.paired_single == 'PAIRED'
                single: meta.paired_single == 'SINGLE'
            }
            .set { trim_branch }

        CUTADAPT_TRIM_PAIRED(trim_branch.paired)
        CUTADAPT_TRIM_SINGLE(trim_branch.single)

        trimmed_ch = CUTADAPT_TRIM_PAIRED.out.trimmed
            .mix(CUTADAPT_TRIM_SINGLE.out.trimmed)
        cutadapt_logs = CUTADAPT_TRIM_PAIRED.out.log
            .mix(CUTADAPT_TRIM_SINGLE.out.log)
    } else {
        trimmed_ch = demux_ch
        cutadapt_logs = Channel.empty()
    }

    // ── 2. Quality Summary (diagnostic: confirm flat-Q30 + read lengths) ─
    DEMUX_SUMMARIZE(trimmed_ch)

    // ── 3. Deblur denoise ────────────────────────────────────────────────
    denoise_input = trimmed_ch
        .join(DEMUX_SUMMARIZE.out.quality_stats)     // [meta, demux.qza, quality.json]

    denoise_input
        .branch { meta, demux, quality ->
            join:    wantsJoin(meta)
            forward: true
        }
        .set { deblur_branch }

    DEBLUR_MERGE_AND_DENOISE(deblur_branch.join)
    DEBLUR_DENOISE_FWD(deblur_branch.forward)

    table_final    = DEBLUR_MERGE_AND_DENOISE.out.table.mix(DEBLUR_DENOISE_FWD.out.table)
    rep_seqs_final = DEBLUR_MERGE_AND_DENOISE.out.rep_seqs.mix(DEBLUR_DENOISE_FWD.out.rep_seqs)
    stats_final    = DEBLUR_MERGE_AND_DENOISE.out.stats.mix(DEBLUR_DENOISE_FWD.out.stats)

    DEBLUR_VISUALIZE_STATS(stats_final)

    // ── 4. Taxonomy Classification ──────────────────────────────────────
    classifiers_keyed = classifiers_ch
        .map { primer_meta, classifier -> [primer_meta.region_id, classifier] }

    rep_seqs_keyed = rep_seqs_final
        .map { meta, rep_seqs -> [meta.primer_region_id, meta, rep_seqs] }

    classify_input = rep_seqs_keyed
        .combine(classifiers_keyed, by: 0)
        .map { region_id, meta, rep_seqs, classifier -> [meta, rep_seqs, classifier] }

    CLASSIFY_TAXONOMY(classify_input)

    // ── 5. Filter to Bacteria only (sOTU level, before collapse) ─────────
    filter_input = table_final
        .join(CLASSIFY_TAXONOMY.out.taxonomy)

    FILTER_BACTERIA(filter_input)

    // ── 6. Collapse to Genus ─────────────────────────────────────────────
    collapse_input = FILTER_BACTERIA.out.bacteria_table
        .join(CLASSIFY_TAXONOMY.out.taxonomy)

    COLLAPSE_TAXONOMY(collapse_input)

    emit:
    tables       = table_final                        // [meta, table.qza]
    rep_seqs     = rep_seqs_final                      // [meta, rep_seqs.qza]
    taxonomy     = CLASSIFY_TAXONOMY.out.taxonomy      // [meta, taxonomy.qza]
    genus_tables = COLLAPSE_TAXONOMY.out.genus_table   // [meta, genus_table.qza]
    quality_viz  = DEMUX_SUMMARIZE.out.visualization    // [meta, summary.qzv]
    stats_viz    = DEBLUR_VISUALIZE_STATS.out.stats_viz // [meta, deblur_stats.qzv]
    cutadapt_logs = cutadapt_logs                       // [meta, cutadapt.log]
}
