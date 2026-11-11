/*
 * Subworkflow: per_study_processing
 * Processes each study unit: import → quality → denoise → QC → taxonomy → collapse.
 * Handles both paired-end and single-end studies.
 */

include { QIIME_IMPORT_PAIRED  } from '../modules/qiime_import/main'
include { QIIME_IMPORT_SINGLE  } from '../modules/qiime_import/main'
include { CUTADAPT_TRIM_PAIRED } from '../modules/cutadapt/main'
include { CUTADAPT_TRIM_SINGLE } from '../modules/cutadapt/main'
include { DEMUX_SUMMARIZE      } from '../modules/quality_summary/main'
include { DADA2_DENOISE_PAIRED } from '../modules/denoise/main'
include { DADA2_DENOISE_SINGLE } from '../modules/denoise/main'
include { CHECK_DENOISE_STATS  } from '../modules/denoise_qc/main'
// Aliased copies for the adaptive second pass (a process can be called once).
include { DADA2_DENOISE_PAIRED as DADA2_DENOISE_PAIRED_RETRY } from '../modules/denoise/main'
include { DADA2_DENOISE_SINGLE as DADA2_DENOISE_SINGLE_RETRY } from '../modules/denoise/main'
include { CHECK_DENOISE_STATS  as CHECK_DENOISE_STATS_RETRY  } from '../modules/denoise_qc/main'
include { ADAPT_DENOISE_PARAMS } from '../modules/denoise_retry/main'
include { CLASSIFY_TAXONOMY    } from '../modules/taxonomy/main'
include { FILTER_BACTERIA      } from '../modules/taxonomy/main'
include { COLLAPSE_TAXONOMY    } from '../modules/taxonomy/main'
include { TABULATE_STATS       } from '../modules/taxonomy/main'


workflow PER_STUDY_PROCESSING {

    take:
    studies_ch         // channel: [meta, manifest]
    classifiers_ch     // channel: [primer_meta, classifier.qza]

    main:

    // ── 1. Import ───────────────────────────────────────────────────────
    // Split by paired/single
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

    // ── 1b. Primer / adapter removal (cutadapt) ─────────────────────────
    // Strip primers (and 3' read-through) by sequence before quality summary,
    // so downstream quality profiles and DADA2 see primer-free reads.
    // Gated on params.run_primer_trim; when disabled, pass demux through.
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

    // ── 2. Quality Summary + Auto-params ────────────────────────────────
    DEMUX_SUMMARIZE(trimmed_ch)

    // ── 3. Denoise (pass 1) ─────────────────────────────────────────────
    denoise_input = trimmed_ch
        .join(DEMUX_SUMMARIZE.out.quality_stats)
        // Results in: [meta, demux.qza, quality_stats.json]

    denoise_input
        .branch { meta, demux, quality ->
            paired: meta.paired_single == 'PAIRED'
            single: meta.paired_single == 'SINGLE'
        }
        .set { denoise_branch }

    DADA2_DENOISE_PAIRED(denoise_branch.paired)
    DADA2_DENOISE_SINGLE(denoise_branch.single)

    table1 = DADA2_DENOISE_PAIRED.out.table.mix(DADA2_DENOISE_SINGLE.out.table)
    rep1   = DADA2_DENOISE_PAIRED.out.rep_seqs.mix(DADA2_DENOISE_SINGLE.out.rep_seqs)
    stats1 = DADA2_DENOISE_PAIRED.out.stats.mix(DADA2_DENOISE_SINGLE.out.stats)

    // ── 4. Denoise QC (pass 1) ──────────────────────────────────────────
    CHECK_DENOISE_STATS(stats1.join(DEMUX_SUMMARIZE.out.quality_stats))
    report1   = CHECK_DENOISE_STATS.out.report
    warnings1 = CHECK_DENOISE_STATS.out.warnings

    // ── 4b. Adaptive one-shot retry for catastrophic studies ────────────
    if (params.retry_denoise) {
        // Flag each pass-1 report with catastrophic + failure_mode.
        report1
            .map { meta, rpt ->
                def j = new groovy.json.JsonSlurper().parseText(rpt.text)
                [meta + [catastrophic: (j.catastrophic ?: false),
                         failure_mode: j.failure_mode], rpt]
            }
            .branch { meta, rpt ->
                retry: meta.catastrophic
                ok:    true
            }
            .set { qc_split }

        // Studies whose denoise crashed (empty table → absent from stats1).
        crashed = denoise_input
            .map { meta, demux, q -> [meta.id, [meta, demux, q]] }
            .join(stats1.map { meta, s -> [meta.id, s] }, by: 0, remainder: true)
            .filter { it[1] != null && it[2] == null }
            .map { id, payload, s ->
                [payload[0] + [catastrophic: true, failure_mode: 'EMPTY'],
                 payload[1], payload[2]]
            }

        // Re-fetch demux+quality for the catastrophic (non-crashed) studies.
        retry_qc = qc_split.retry
            .map { meta, rpt -> [meta.id, meta] }
            .join(denoise_input.map { meta, demux, q -> [meta.id, demux, q] }, by: 0)
            .map { id, meta, demux, q -> [meta, demux, q] }

        retry_input = retry_qc.mix(crashed)   // [meta(+failure_mode), demux, quality]

        ADAPT_DENOISE_PARAMS(retry_input.map { meta, demux, q -> [meta, q] })

        // Turn adapted.json into meta.denoise_override + paired/single mode,
        // then re-attach demux+quality.
        retry_ready = ADAPT_DENOISE_PARAMS.out.adapted
            .map { meta, aj ->
                def a = new groovy.json.JsonSlurper().parseText(aj.text)
                def ov = [:]
                if (a.trunc_len_f != null) ov['trunc_len_f'] = a.trunc_len_f.toString()
                if (a.trunc_len_r != null) ov['trunc_len_r'] = a.trunc_len_r.toString()
                if (a.min_overlap != null) ov['min_overlap'] = a.min_overlap.toString()
                if (a.max_ee     != null) ov['max_ee']      = a.max_ee.toString()
                if (a.trunc_q    != null) ov['trunc_q']     = a.trunc_q.toString()
                def m2 = meta + [denoise_override: ov,
                                 paired_single: (a.mode == 'single' ? 'SINGLE' : 'PAIRED'),
                                 retried: true]
                [meta.id, m2]
            }
            .join(retry_input.map { meta, demux, q -> [meta.id, demux, q] }, by: 0)
            .map { id, m2, demux, q -> [m2, demux, q] }

        retry_ready
            .branch { meta, demux, q ->
                paired: meta.paired_single == 'PAIRED'
                single: meta.paired_single == 'SINGLE'
            }
            .set { retry_branch }

        DADA2_DENOISE_PAIRED_RETRY(retry_branch.paired)
        DADA2_DENOISE_SINGLE_RETRY(retry_branch.single)

        table2 = DADA2_DENOISE_PAIRED_RETRY.out.table.mix(DADA2_DENOISE_SINGLE_RETRY.out.table)
        rep2   = DADA2_DENOISE_PAIRED_RETRY.out.rep_seqs.mix(DADA2_DENOISE_SINGLE_RETRY.out.rep_seqs)
        stats2 = DADA2_DENOISE_PAIRED_RETRY.out.stats.mix(DADA2_DENOISE_SINGLE_RETRY.out.stats)

        // QC pass 2 (reuse the original quality JSON; key by id since meta changed).
        CHECK_DENOISE_STATS_RETRY(
            stats2.map { meta, s -> [meta.id, meta, s] }
                  .join(DEMUX_SUMMARIZE.out.quality_stats.map { meta, q -> [meta.id, q] }, by: 0)
                  .map { id, meta, s, q -> [meta, s, q] }
        )

        // Combine, FAIL-SAFE: pass-1 is the baseline for every (non-crashed)
        // study; pass-2 OVERRIDES it only where it actually materialised. So if
        // ADAPT or the pass-2 denoise fails, the study keeps its pass-1 output
        // instead of vanishing. A crashed study (no pass-1) appears only via a
        // successful pass-2. Per artifact: remainder-join pass-1 (left) with
        // pass-2 (right), preferring pass-2 when present (it[3] != null).
        table_final = table1.map { meta, x -> [meta.id, meta, x] }
            .join(table2.map { meta, x -> [meta.id, meta, x] }, by: 0, remainder: true)
            .map { it -> (it[3] != null) ? [it[3], it[4]] : [it[1], it[2]] }
        rep_seqs_final = rep1.map { meta, x -> [meta.id, meta, x] }
            .join(rep2.map { meta, x -> [meta.id, meta, x] }, by: 0, remainder: true)
            .map { it -> (it[3] != null) ? [it[3], it[4]] : [it[1], it[2]] }
        stats_final = stats1.map { meta, x -> [meta.id, meta, x] }
            .join(stats2.map { meta, x -> [meta.id, meta, x] }, by: 0, remainder: true)
            .map { it -> (it[3] != null) ? [it[3], it[4]] : [it[1], it[2]] }
        warnings_final = warnings1.map { meta, x -> [meta.id, meta, x] }
            .join(CHECK_DENOISE_STATS_RETRY.out.warnings.map { meta, x -> [meta.id, meta, x] }, by: 0, remainder: true)
            .map { it -> (it[3] != null) ? [it[3], it[4]] : [it[1], it[2]] }
        report_final = report1.map { meta, x -> [meta.id, meta, x] }
            .join(CHECK_DENOISE_STATS_RETRY.out.report.map { meta, x -> [meta.id, meta, x] }, by: 0, remainder: true)
            .map { it -> (it[3] != null) ? [it[3], it[4]] : [it[1], it[2]] }
    } else {
        table_final    = table1
        rep_seqs_final = rep1
        stats_final    = stats1
        report_final   = report1
        warnings_final = warnings1
    }

    TABULATE_STATS(stats_final)

    // ── 5. Taxonomy Classification ──────────────────────────────────────
    classifiers_keyed = classifiers_ch
        .map { primer_meta, classifier -> [primer_meta.region_id, classifier] }

    rep_seqs_keyed = rep_seqs_final
        .map { meta, rep_seqs -> [meta.primer_region_id, meta, rep_seqs] }

    classify_input = rep_seqs_keyed
        .combine(classifiers_keyed, by: 0)
        .map { region_id, meta, rep_seqs, classifier ->
            [meta, rep_seqs, classifier]
        }

    CLASSIFY_TAXONOMY(classify_input)

    // ── 6. Filter to Bacteria only (ASV level, before collapse) ────────────
    filter_input = table_final
        .join(CLASSIFY_TAXONOMY.out.taxonomy)

    FILTER_BACTERIA(filter_input)

    // ── 7. Collapse to Genus ────────────────────────────────────────────────
    collapse_input = FILTER_BACTERIA.out.bacteria_table
        .join(CLASSIFY_TAXONOMY.out.taxonomy)

    COLLAPSE_TAXONOMY(collapse_input)

    emit:
    tables       = table_final                            // [meta, table.qza]
    rep_seqs     = rep_seqs_final                          // [meta, rep_seqs.qza]
    taxonomy     = CLASSIFY_TAXONOMY.out.taxonomy          // [meta, taxonomy.qza]
    genus_tables = COLLAPSE_TAXONOMY.out.genus_table       // [meta, genus_table.qza]
    qc_reports   = report_final                            // [meta, qc_report.json]
    qc_warnings  = warnings_final                          // [meta, qc_warnings.txt]
    quality_viz  = DEMUX_SUMMARIZE.out.visualization        // [meta, summary.qzv]
    stats_viz    = TABULATE_STATS.out.stats_viz             // [meta, stats.qzv]
    cutadapt_logs = cutadapt_logs                           // [meta, cutadapt.log]
}
