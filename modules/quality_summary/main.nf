/*
 * Module: quality_summary
 * Runs demux summarize and extracts quality stats for auto-parametrization.
 */

process DEMUX_SUMMARIZE {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/quality", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.qzv') ? fn : null }

    input:
    tuple val(meta), path(demux_qza)

    output:
    tuple val(meta), path("${meta.id}_summary.qzv"),          emit: visualization
    tuple val(meta), path("${meta.id}_quality_stats.json"),    emit: quality_stats

    script:
    // When cutadapt already removed primers (run_primer_trim), reads are
    // primer-free: pass 0 so the truncation floor and trim_left output are 0
    // and DADA2 does not trim again.
    def primer_f_len = params.run_primer_trim ? 0 : (meta.primer_len_f ?: 0)
    def primer_r_len = params.run_primer_trim ? 0 : (meta.primer_len_r ?: 0)
    def region = meta.region ?: ''
    """
    qiime demux summarize \\
        --i-data ${demux_qza} \\
        --p-n ${params.quality_sample_n} \\
        --o-visualization ${meta.id}_summary.qzv

    parse_quality_stats.py \\
        --qzv ${meta.id}_summary.qzv \\
        --output ${meta.id}_quality_stats.json \\
        --min-quality ${params.quality_drop_threshold} \\
        --min-quality-p9 ${params.quality_p9_threshold} \\
        --primer-len-f ${primer_f_len} \\
        --primer-len-r ${primer_r_len} \\
        --region '${region}' \\
        --percentile ${params.quality_percentile} \\
        --window ${params.quality_window} \\
        --tail-check-window ${params.tail_check_window} \\
        --retry-quality-floor ${params.retry_quality_floor}
    """
}
