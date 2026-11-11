/*
 * Module: qiime_import
 * Imports demultiplexed FASTQ files into QIIME2 artifacts (.qza).
 * Supports both paired-end and single-end formats.
 */

process QIIME_IMPORT_PAIRED {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    input:
    tuple val(meta), path(manifest)

    output:
    tuple val(meta), path("${meta.id}_demux.qza"), emit: demux

    script:
    """
    qiime tools import \\
        --type 'SampleData[PairedEndSequencesWithQuality]' \\
        --input-format PairedEndFastqManifestPhred${params.phred_offset}V2 \\
        --input-path ${manifest} \\
        --output-path ${meta.id}_demux.qza
    """
}


process QIIME_IMPORT_SINGLE {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    input:
    tuple val(meta), path(manifest)

    output:
    tuple val(meta), path("${meta.id}_demux.qza"), emit: demux

    script:
    """
    qiime tools import \\
        --type 'SampleData[SequencesWithQuality]' \\
        --input-format SingleEndFastqManifestPhred${params.phred_offset}V2 \\
        --input-path ${manifest} \\
        --output-path ${meta.id}_demux.qza
    """
}
