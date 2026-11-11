/*
 * Module: taxonomy
 * Taxonomic classification using sklearn classifier and collapse to genus level.
 */

process CLASSIFY_TAXONOMY {
    tag "${meta.id}"
    label 'process_medium'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/taxonomy", mode: 'copy'

    input:
    tuple val(meta), path(rep_seqs), path(classifier)

    output:
    tuple val(meta), path("${meta.id}_taxonomy.qza"),    emit: taxonomy
    tuple val(meta), path("${meta.id}_taxonomy.qzv"),    emit: taxonomy_viz

    script:
    """
    qiime feature-classifier classify-sklearn \\
        --i-reads ${rep_seqs} \\
        --i-classifier ${classifier} \\
        --p-confidence ${params.classify_confidence} \\
        --o-classification ${meta.id}_taxonomy.qza

    qiime metadata tabulate \\
        --m-input-file ${meta.id}_taxonomy.qza \\
        --o-visualization ${meta.id}_taxonomy.qzv
    """
}


process FILTER_BACTERIA {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/taxonomy", mode: 'copy'

    input:
    tuple val(meta), path(table_qza), path(taxonomy_qza)

    output:
    tuple val(meta), path("${meta.id}_bacteria_table.qza"), emit: bacteria_table

    script:
    """
    # Filter at ASV level (before collapse) where feature IDs match taxonomy QZA
    qiime taxa filter-table \\
        --i-table ${table_qza} \\
        --i-taxonomy ${taxonomy_qza} \\
        --p-include 'd__Bacteria' \\
        --o-filtered-table ${meta.id}_bacteria_table.qza
    """
}


process COLLAPSE_TAXONOMY {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/taxonomy", mode: 'copy'

    input:
    tuple val(meta), path(table_qza), path(taxonomy_qza)

    output:
    tuple val(meta), path("${meta.id}_genus_table.qza"), emit: genus_table

    script:
    """
    qiime taxa collapse \\
        --i-table ${table_qza} \\
        --i-taxonomy ${taxonomy_qza} \\
        --p-level ${params.collapse_level} \\
        --o-collapsed-table ${meta.id}_genus_table.qza
    """
}


process TABULATE_STATS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/denoise", mode: 'copy'

    input:
    tuple val(meta), path(stats_qza)

    output:
    tuple val(meta), path("${meta.id}_stats.qzv"), emit: stats_viz

    script:
    """
    qiime metadata tabulate \\
        --m-input-file ${stats_qza} \\
        --o-visualization ${meta.id}_stats.qzv
    """
}
