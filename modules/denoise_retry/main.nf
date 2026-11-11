/*
 * Module: denoise_retry
 * Computes second-pass DADA2 parameters for a study whose first-pass denoise
 * was catastrophic. Emits a small JSON (adapted.json) that the subworkflow turns
 * into meta.denoise_override and uses to route the retry to paired vs single.
 */

process ADAPT_DENOISE_PARAMS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy'

    input:
    tuple val(meta), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_adapted_params.json"), emit: adapted

    script:
    def paired_flag = (meta.paired_single == 'PAIRED') ? '--paired' : ''
    def region = meta.region ?: ''
    """
    adapt_denoise_params.py \\
        --quality-json ${quality_json} \\
        --failure-mode ${meta.failure_mode} \\
        ${paired_flag} \\
        --region '${region}' \\
        --min-overlap ${params.min_overlap} \\
        --output ${meta.id}_adapted_params.json
    """
}
