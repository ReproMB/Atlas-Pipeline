/*
 * Module: classifier_build
 * Builds QIIME2 Naive Bayes classifiers from SILVA reference sequences.
 * One classifier per unique (deduplicated) primer combination.
 * Supports caching: pre-built classifiers in classifier_cache_dir are reused.
 */

// Expected min/max amplicon lengths by region for extract-reads
def REGION_LENGTHS = [
    'V1_V2': [min: 200, max: 400],
    'V1_V3': [min: 400, max: 600],
    'V1_V5': [min: 800, max: 1400],
    'V2_V3': [min: 200, max: 400],
    'V3':    [min: 150, max: 250],
    'V3_V4': [min: 350, max: 500],
    'V3_V5': [min: 450, max: 700],
    'V4':    [min: 200, max: 310],
    'V4_V5': [min: 250, max: 400],
]


process EXTRACT_READS {
    tag "${primer_meta.region_id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/classifiers/ref_seqs", mode: 'copy'

    input:
    tuple val(primer_meta), path(silva_seqs)

    output:
    tuple val(primer_meta), path("ref_seqs_${primer_meta.region_id}.qza"), emit: ref_reads

    script:
    def region_key = primer_meta.region?.replaceAll('-', '_') ?: 'V3_V4'
    def lengths = REGION_LENGTHS[region_key] ?: [min: 100, max: 600]
    def min_len = primer_meta.min_length ?: lengths.min
    def max_len = primer_meta.max_length ?: lengths.max
    """
    qiime feature-classifier extract-reads \\
        --i-sequences ${silva_seqs} \\
        --p-f-primer '${primer_meta.fwd_primer}' \\
        --p-r-primer '${primer_meta.rev_primer}' \\
        --p-trim-right 0 \\
        --p-trunc-len 0 \\
        --p-trim-left 0 \\
        --p-identity 0.8 \\
        --p-min-length ${min_len} \\
        --p-max-length ${max_len} \\
        --p-n-jobs ${task.cpus} \\
        --p-read-orientation both \\
        --o-reads ref_seqs_${primer_meta.region_id}.qza
    """
}


process FIT_CLASSIFIER {
    tag "${primer_meta.region_id}"
    label 'process_high'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/classifiers", mode: 'copy'

    input:
    tuple val(primer_meta), path(ref_reads), path(silva_tax)

    output:
    tuple val(primer_meta), path("classifier_${primer_meta.region_id}.qza"), emit: classifier

    script:
    """
    qiime feature-classifier fit-classifier-naive-bayes \\
        --i-reference-reads ${ref_reads} \\
        --i-reference-taxonomy ${silva_tax} \\
        --p-classify--alpha ${params.classifier_alpha} \\
        --o-classifier classifier_${primer_meta.region_id}.qza
    """
}

