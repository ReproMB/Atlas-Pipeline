/*
 * Module: cutadapt
 * Sequence-based primer removal (and 3' read-through removal) with cutadapt,
 * run BEFORE quality summary / DADA2. Replaces DADA2 positional trim-left.
 *
 * 5' primers are removed via --p-front-{f,r}. 3' read-through (a read running
 * into the reverse complement of the opposite primer) is removed via
 * --p-adapter-{f,r}, fed the reverse-complement primers carried on meta.
 *
 * Reads lacking the primer are, by default, passed through untrimmed rather
 * than discarded (params.cutadapt_discard_untrimmed) so already-primer-trimmed
 * submissions are not wiped. The cutadapt log is emitted for QC.
 */

process CUTADAPT_TRIM_PAIRED {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/primer_trim", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? fn : null }

    input:
    tuple val(meta), path(demux_qza)

    output:
    tuple val(meta), path("${meta.id}_trimmed.qza"), emit: trimmed
    tuple val(meta), path("${meta.id}_cutadapt.log"), emit: log

    script:
    def discard = params.cutadapt_discard_untrimmed ? '--p-discard-untrimmed' : '--p-no-discard-untrimmed'
    """
    qiime cutadapt trim-paired \\
        --i-demultiplexed-sequences ${demux_qza} \\
        --p-front-f ${meta.fwd_primer} \\
        --p-front-r ${meta.rev_primer} \\
        --p-adapter-f ${meta.rev_primer_rc} \\
        --p-adapter-r ${meta.fwd_primer_rc} \\
        --p-error-rate ${params.cutadapt_error_rate} \\
        --p-match-read-wildcards \\
        --p-minimum-length ${params.cutadapt_min_length} \\
        ${discard} \\
        --p-cores ${task.cpus} \\
        --o-trimmed-sequences ${meta.id}_trimmed.qza \\
        --verbose 2>&1 | tee ${meta.id}_cutadapt.log
    """
}


process CUTADAPT_TRIM_SINGLE {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/primer_trim", mode: 'copy',
        saveAs: { fn -> fn.endsWith('.log') ? fn : null }

    input:
    tuple val(meta), path(demux_qza)

    output:
    tuple val(meta), path("${meta.id}_trimmed.qza"), emit: trimmed
    tuple val(meta), path("${meta.id}_cutadapt.log"), emit: log

    script:
    def discard = params.cutadapt_discard_untrimmed ? '--p-discard-untrimmed' : '--p-no-discard-untrimmed'
    """
    qiime cutadapt trim-single \\
        --i-demultiplexed-sequences ${demux_qza} \\
        --p-front ${meta.fwd_primer} \\
        --p-adapter ${meta.rev_primer_rc} \\
        --p-error-rate ${params.cutadapt_error_rate} \\
        --p-match-read-wildcards \\
        --p-minimum-length ${params.cutadapt_min_length} \\
        ${discard} \\
        --p-cores ${task.cpus} \\
        --o-trimmed-sequences ${meta.id}_trimmed.qza \\
        --verbose 2>&1 | tee ${meta.id}_cutadapt.log
    """
}
