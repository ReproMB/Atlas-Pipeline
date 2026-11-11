/*
 * Module: deblur
 * ~~~~~~~~~~~~~~~
 * Denoising path for "false FASTQ" studies — submissions whose per-base quality
 * scores are synthetic (e.g. a flat Q30 line) rather than real instrument Phred
 * values. DADA2 cannot be trusted on these because its error model is *learned
 * from the quality scores*; a constant Q30 either crashes the fit or produces a
 * meaningless posterior. Deblur is used instead: its error model is a STATIC
 * positional profile applied above a fixed trim length, so it never reads the
 * quality scores and cannot be corrupted by them.
 *
 * Two entry points, chosen per region (see per_study_processing_deblur.nf):
 *
 *   DEBLUR_MERGE_AND_DENOISE  (paired, overlapping region → join)
 *       vsearch merge-pairs → quality-filter q-score → deblur denoise-16S
 *       Merging is alignment-driven, so a constant Q30 does not break the
 *       overlap resolution (the merged posterior quality is meaningless but
 *       unused).
 *
 *   DEBLUR_DENOISE_FWD        (single-end, or long region processed forward-only)
 *       quality-filter q-score → deblur denoise-16S
 *
 * NOTE on the quality-filter step: at a uniform Q30 it passes every read — it is
 * INERT here, exactly as DADA2's maxEE was. It is kept for recipe fidelity and
 * to emit filter-stats, but it provides NO front-end junk removal for these
 * studies. Compensation lives entirely in (a) deblur's positive 16S recruitment
 * filter and (b) a raised --p-min-reads (params.deblur_min_reads, default 50, vs
 * the QIIME2 default of 10).
 *
 * TRIM LENGTH is the parameter that governs cross-study comparability and it
 * CANNOT be read off a quality profile here (the profiles are flat lines). It is
 * resolved, per study, from geometry — see resolveTrimLength() — and MUST be a
 * single value within a region group so that identical sequences collapse to
 * identical features. Reads shorter than the trim are dropped; longer ones are
 * truncated to it.
 *
 * Outputs mirror the DADA2 module exactly (table / rep_seqs / stats) so the
 * downstream taxonomy → filter → collapse chain is byte-for-byte unchanged.
 * Deblur sOTUs and DADA2 ASVs are the same kind of object: denoised exact
 * sequences classified against the same SILVA classifier.
 */


// ── Trim-length resolution (geometry-driven, NOT quality-driven) ─────────────
// Precedence: per-study manual override  →  per-region map  →  global default.
// A non-positive result is a hard error: silently disabling the trim (-1) would
// let each study define its own feature space and destroy comparability, which
// is the whole point of using a fixed trim here.
def resolveTrimLength = { meta ->
    def ov     = meta.denoise_override
    def byReg  = (params.deblur_trim_length_by_region ?: [:])
    def picked = (ov?.deblur_trim_length ?: meta.deblur_trim_length
                  ?: byReg[meta.region] ?: params.deblur_trim_length)
    def n = (picked != null && "${picked}".isInteger()) ? "${picked}".toInteger() : 0
    if (n <= 0) {
        error """Study ${meta.id}: no Deblur trim-length resolved for region '${meta.region}'.
Set one (single value per region group) via, in order of precedence:
  * manual_params column  deblur_trim_length
  * params.deblur_trim_length_by_region = [ '${meta.region}': <bp> , ... ]
  * params.deblur_trim_length = <bp>
Pick it from the read-length distribution / amplicon geometry, not from quality."""
    }
    return n
}


process DEBLUR_MERGE_AND_DENOISE {
    tag "${meta.id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(demux_qza), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_table.qza"),    emit: table
    tuple val(meta), path("${meta.id}_rep_seqs.qza"), emit: rep_seqs
    tuple val(meta), path("${meta.id}_stats.qza"),    emit: stats

    script:
    def trim = resolveTrimLength(meta)
    // minovlen falls back to the DADA2 value; maxdiffs has its own default (10,
    // see nextflow.config). Explicit null checks so a 0 override is honoured.
    def merge_minovlen = (params.deblur_merge_minovlen != null) ? params.deblur_merge_minovlen : params.min_overlap
    def merge_maxdiffs = (params.deblur_merge_maxdiffs != null) ? params.deblur_merge_maxdiffs : params.max_merge_mismatch
    """
    #!/bin/bash
    set -euo pipefail

    echo "Deblur[join] ${meta.id}: region=${meta.region} trim_length=${trim} \
merge_minovlen=${merge_minovlen} merge_maxdiffs=${merge_maxdiffs} \
min_reads=${params.deblur_min_reads} min_size=${params.deblur_min_size} \
min_quality=${params.deblur_min_quality} (quality scores intentionally ignored)"

    # ── 1. Merge overlapping pairs (alignment-driven; quality-agnostic) ──────
    # minovlen is DADA2's min_overlap (12). maxdiffs is relaxed to 10 because
    # Deblur, being single-read, must merge RAW reads, whereas DADA2 merges after
    # error correction; zero mismatches rejected genuine overlaps.
    qiime vsearch merge-pairs \\
        --i-demultiplexed-seqs ${demux_qza} \\
        --p-minovlen ${merge_minovlen} \\
        --p-maxdiffs ${merge_maxdiffs} \\
        --p-allowmergestagger \\
        --p-threads ${task.cpus} \\
        --o-merged-sequences ${meta.id}_merged.qza \\
        --o-unmerged-sequences ${meta.id}_unmerged.qza \\
        --verbose

    # ── 2. quality-filter q-score (INERT at flat Q30; kept for recipe/stats) ─
    qiime quality-filter q-score \\
        --i-demux ${meta.id}_merged.qza \\
        --p-min-quality ${params.deblur_min_quality} \\
        --o-filtered-sequences ${meta.id}_qfilt.qza \\
        --o-filter-stats ${meta.id}_qfilt_stats.qza

    # ── 3. Deblur denoise (static positional error model; no quality read) ───
    qiime deblur denoise-16S \\
        --i-demultiplexed-seqs ${meta.id}_qfilt.qza \\
        --p-trim-length ${trim} \\
        --p-min-reads ${params.deblur_min_reads} \\
        --p-min-size ${params.deblur_min_size} \\
        --p-sample-stats \\
        --p-jobs-to-start ${task.cpus} \\
        --o-representative-sequences ${meta.id}_rep_seqs.qza \\
        --o-table ${meta.id}_table.qza \\
        --o-stats ${meta.id}_stats.qza \\
        --verbose
    """
}


process DEBLUR_DENOISE_FWD {
    tag "${meta.id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(demux_qza), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_table.qza"),    emit: table
    tuple val(meta), path("${meta.id}_rep_seqs.qza"), emit: rep_seqs
    tuple val(meta), path("${meta.id}_stats.qza"),    emit: stats

    script:
    def trim = resolveTrimLength(meta)
    // deblur denoise-16S reads only the forward strand of a paired artifact, so a
    // paired study forced forward-only (long, non-overlapping region) needs no
    // re-import — the same demux artifact is passed straight through.
    """
    #!/bin/bash
    set -euo pipefail

    echo "Deblur[fwd] ${meta.id}: region=${meta.region} trim_length=${trim} \
min_reads=${params.deblur_min_reads} min_size=${params.deblur_min_size} \
min_quality=${params.deblur_min_quality} (quality scores intentionally ignored)"

    # ── 1. quality-filter q-score (INERT at flat Q30; kept for recipe/stats) ─
    qiime quality-filter q-score \\
        --i-demux ${demux_qza} \\
        --p-min-quality ${params.deblur_min_quality} \\
        --o-filtered-sequences ${meta.id}_qfilt.qza \\
        --o-filter-stats ${meta.id}_qfilt_stats.qza

    # ── 2. Deblur denoise (static positional error model; no quality read) ───
    qiime deblur denoise-16S \\
        --i-demultiplexed-seqs ${meta.id}_qfilt.qza \\
        --p-trim-length ${trim} \\
        --p-min-reads ${params.deblur_min_reads} \\
        --p-min-size ${params.deblur_min_size} \\
        --p-sample-stats \\
        --p-jobs-to-start ${task.cpus} \\
        --o-representative-sequences ${meta.id}_rep_seqs.qza \\
        --o-table ${meta.id}_table.qza \\
        --o-stats ${meta.id}_stats.qza \\
        --verbose
    """
}


// Deblur emits DeblurStats, which `qiime metadata tabulate` cannot render.
// Use the dedicated visualizer instead of the DADA2 TABULATE_STATS process.
process DEBLUR_VISUALIZE_STATS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(stats_qza)

    output:
    tuple val(meta), path("${meta.id}_deblur_stats.qzv"), emit: stats_viz

    script:
    """
    qiime deblur visualize-stats \\
        --i-deblur-stats ${stats_qza} \\
        --o-visualization ${meta.id}_deblur_stats.qzv
    """
}
