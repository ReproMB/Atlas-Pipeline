/*
 * Module: downstream
 * Post-processing: export of the genus tables used by the Atlas, and per-study
 * summaries (alpha/beta diversity on a feature-filtered copy, Lactobacillus
 * abundance). The summaries are not used to build the Atlas tables.
 */

process EXPORT_BIOM {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/exports", mode: 'copy'

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("${meta.id}_feature_table.biom"),  emit: biom
    tuple val(meta), path("${meta.id}_feature_table.tsv"),   emit: tsv

    script:
    """
    # Table is already bacteria-only (filtered upstream at ASV level)
    qiime tools export \\
        --input-path ${table_qza} \\
        --output-path export_dir/

    mv export_dir/feature-table.biom ${meta.id}_feature_table.biom

    biom convert \\
        -i ${meta.id}_feature_table.biom \\
        -o ${meta.id}_feature_table.tsv \\
        --to-tsv
    """
}

process EXPORT_TAXONOMY {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/exports", mode: 'copy'

    input:
    tuple val(meta), path(taxonomy_qza)

    output:
    tuple val(meta), path("${meta.id}_taxonomy.tsv"), emit: taxonomy_tsv

    script:
    """
    qiime tools export \\
        --input-path ${taxonomy_qza} \\
        --output-path tax_export/

    mv tax_export/taxonomy.tsv ${meta.id}_taxonomy.tsv
    """
}

process FILTER_BY_PREVALENCE {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/filtered", mode: 'copy'

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("${meta.id}_filtered_table.qza"), emit: filtered_table

    script:
    """
    qiime feature-table filter-features \\
        --i-table ${table_qza} \\
        --p-min-samples 2 \\
        --p-min-frequency 10 \\
        --o-filtered-table ${meta.id}_filtered_table.qza
    """
}

process ALPHA_DIVERSITY {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    when:
    params.run_diversity

    publishDir "${params.outdir}/${meta.id}/diversity/alpha", mode: 'copy'

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("${meta.id}_shannon.qza"),  emit: shannon
    tuple val(meta), path("${meta.id}_observed.qza"), emit: observed
    tuple val(meta), path("${meta.id}_simpson.qza"),  emit: simpson

    script:
    """
    qiime diversity alpha \\
        --i-table ${table_qza} \\
        --p-metric shannon \\
        --o-alpha-diversity ${meta.id}_shannon.qza

    qiime diversity alpha \\
        --i-table ${table_qza} \\
        --p-metric observed_features \\
        --o-alpha-diversity ${meta.id}_observed.qza

    qiime diversity alpha \\
        --i-table ${table_qza} \\
        --p-metric simpson \\
        --o-alpha-diversity ${meta.id}_simpson.qza
    """
}

process ALPHA_GROUP_SIGNIFICANCE {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    when:
    params.run_diversity

    publishDir "${params.outdir}/${meta.id}/diversity/alpha_significance", mode: 'copy'

    input:
    tuple val(meta), path(shannon), path(observed), path(simpson), path(metadata_tsv)

    output:
    tuple val(meta), path("${meta.id}_shannon_sig.qzv"),  emit: shannon_sig
    tuple val(meta), path("${meta.id}_observed_sig.qzv"), emit: observed_sig
    tuple val(meta), path("${meta.id}_simpson_sig.qzv"),  emit: simpson_sig

    script:
    """
    qiime diversity alpha-group-significance \\
        --i-alpha-diversity ${shannon} \\
        --m-metadata-file ${metadata_tsv} \\
        --o-visualization ${meta.id}_shannon_sig.qzv || true

    qiime diversity alpha-group-significance \\
        --i-alpha-diversity ${observed} \\
        --m-metadata-file ${metadata_tsv} \\
        --o-visualization ${meta.id}_observed_sig.qzv || true

    qiime diversity alpha-group-significance \\
        --i-alpha-diversity ${simpson} \\
        --m-metadata-file ${metadata_tsv} \\
        --o-visualization ${meta.id}_simpson_sig.qzv || true
    """
}

process BETA_DIVERSITY {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    when:
    params.run_diversity

    publishDir "${params.outdir}/${meta.id}/diversity/beta", mode: 'copy'

    input:
    tuple val(meta), path(table_qza)

    output:
    tuple val(meta), path("${meta.id}_braycurtis.qza"),  emit: braycurtis
    tuple val(meta), path("${meta.id}_jaccard.qza"),     emit: jaccard

    script:
    """
    qiime diversity beta \\
        --i-table ${table_qza} \\
        --p-metric braycurtis \\
        --o-distance-matrix ${meta.id}_braycurtis.qza

    qiime diversity beta \\
        --i-table ${table_qza} \\
        --p-metric jaccard \\
        --o-distance-matrix ${meta.id}_jaccard.qza
    """
}

process BETA_GROUP_SIGNIFICANCE {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    when:
    params.run_diversity

    publishDir "${params.outdir}/${meta.id}/diversity/beta_significance", mode: 'copy'

    input:
    tuple val(meta), path(braycurtis), path(jaccard), path(metadata_tsv)

    output:
    tuple val(meta), path("${meta.id}_braycurtis_cond_sig.qzv"), emit: bray_cond_sig
    tuple val(meta), path("${meta.id}_braycurtis_site_sig.qzv"), emit: bray_site_sig
    tuple val(meta), path("${meta.id}_jaccard_cond_sig.qzv"),    emit: jac_cond_sig
    tuple val(meta), path("${meta.id}_jaccard_site_sig.qzv"),    emit: jac_site_sig

    script:
    """
    # Bray-Curtis
    qiime diversity beta-group-significance \\
        --i-distance-matrix ${braycurtis} \\
        --m-metadata-file ${metadata_tsv} \\
        --m-metadata-column ConditionAtlas \\
        --p-method permanova \\
        --o-visualization ${meta.id}_braycurtis_cond_sig.qzv || true

    qiime diversity beta-group-significance \\
        --i-distance-matrix ${braycurtis} \\
        --m-metadata-file ${metadata_tsv} \\
        --m-metadata-column SiteAtlas \\
        --p-method permanova \\
        --o-visualization ${meta.id}_braycurtis_site_sig.qzv || true

    # Jaccard
    qiime diversity beta-group-significance \\
        --i-distance-matrix ${jaccard} \\
        --m-metadata-file ${metadata_tsv} \\
        --m-metadata-column ConditionAtlas \\
        --p-method permanova \\
        --o-visualization ${meta.id}_jaccard_cond_sig.qzv || true

    qiime diversity beta-group-significance \\
        --i-distance-matrix ${jaccard} \\
        --m-metadata-file ${metadata_tsv} \\
        --m-metadata-column SiteAtlas \\
        --p-method permanova \\
        --o-visualization ${meta.id}_jaccard_site_sig.qzv || true
    """
}

process LACTOBACILLUS_ABUNDANCE {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}/lactobacillus", mode: 'copy'

    input:
    tuple val(meta), path(feature_table), path(taxonomy), path(metadata_tsv)

    output:
    tuple val(meta), path("${meta.id}_lactobacillus_summary.csv"), emit: lacto_csv
    tuple val(meta), path("${meta.id}_lactobacillus_boxplot.png"), emit: lacto_plot

    script:
    """
    summarize_lactobacillus.py \\
        --table ${feature_table} \\
        --tax ${taxonomy} \\
        --meta ${metadata_tsv} \\
        --prefix ${meta.id}
    """
}
