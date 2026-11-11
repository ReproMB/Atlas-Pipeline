/*
 * Module: input_parsing
 * Parses the main study CSV, cleans primers, groups reproductive tract sites,
 * and generates QIIME2 manifests for each processing unit.
 */

process PARSE_AND_CLEAN_INPUT {
    tag "parse_input"
    label 'process_low'

    input:
    path(input_csv)

    output:
    path("cleaned_studies.csv"),        emit: cleaned_csv
    path("primer_dedup_report.txt"),    emit: dedup_report

    script:
    """
    clean_primers.py \\
        --input-csv ${input_csv} \\
        --output-csv cleaned_studies.csv \\
        --dedup-report primer_dedup_report.txt
    """
}

process PARSE_METADATA {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), val(sites_list)
    path(metadata_dir)

    output:
    tuple val(meta), path("${meta.study_project}_metadata.tsv"), emit: metadata_tsv
    tuple val(meta), path("${meta.study_project}_filter.txt"),    emit: sample_filter

    script:
    def sites_arg = sites_list.join(',')
    """
    parse_study_metadata.py \\
        --study ${meta.study_project} \\
        --metadata-dir ${metadata_dir} \\
        --sites '${sites_arg}' \\
        --out-tsv ${meta.study_project}_metadata.tsv \\
        --out-filter ${meta.study_project}_filter.txt
    """
}

process GENERATE_MANIFEST {
    tag "${meta.id}"
    label 'process_low'

    input:
    tuple val(meta), val(sites_list), path(sample_filter)
    path(data_dir)

    output:
    tuple val(meta), path("${meta.id}_manifest.csv"), path("${meta.id}.mode"), emit: manifest

    script:
    def sites_arg = sites_list.join(',')
    def paired_flag = meta.paired_single == 'PAIRED' ? '--paired' : ''
    """
    generate_manifest.py \\
        --study ${meta.study_project} \\
        --sites '${sites_arg}' \\
        --data-dir ${data_dir} \\
        ${paired_flag} \\
        --sample-filter ${sample_filter} \\
        --output ${meta.id}_manifest.csv \\
        --mode-file ${meta.id}.mode
    """
}
