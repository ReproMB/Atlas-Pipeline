/*
 * Module: denoise_qc
 * Evaluates DADA2 denoising statistics and generates QC alerts.
 * Warns users when processing quality is poor.
 */

process CHECK_DENOISE_STATS {
    tag "${meta.id}"
    label 'process_low'
    container "${params.qiime_container}"

    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy'

    input:
    tuple val(meta), path(stats_qza), path(quality_json)

    output:
    tuple val(meta), path("${meta.id}_qc_report.json"),   emit: report
    tuple val(meta), path("${meta.id}_qc_warnings.txt"),  emit: warnings

    script:
    """
    #!/bin/bash
    set -euo pipefail

    # Export DADA2 stats to TSV
    qiime tools export \\
        --input-path ${stats_qza} \\
        --output-path stats_export/

    # Run QC check
    check_denoise_stats.py \\
        --stats stats_export/stats.tsv \\
        --study-id ${meta.id} \\
        --quality-json ${quality_json} \\
        --min-passing-pct ${params.min_reads_passing_pct} \\
        --min-merged-pct ${params.min_reads_merged_pct} \\
        --min-nonchimeric-pct ${params.min_reads_nonchimeric_pct} \\
        --min-reads-output ${params.min_reads_output} \\
        --max-cv ${params.heterogeneity_cv_max} \\
        --retry-filter-pct ${params.retry_filter_pct} \\
        --retry-merged-pct ${params.retry_merged_pct} \\
        --retry-zero-frac ${params.retry_zero_frac} \\
        --output-report ${meta.id}_qc_report.json \\
        --output-warnings ${meta.id}_qc_warnings.txt
    """
}


process AGGREGATE_QC_REPORTS {
    tag "aggregate_qc"
    label 'process_low'

    publishDir "${params.outdir}/qc_summary", mode: 'copy'

    input:
    path(reports)
    path(warnings)

    output:
    path("qc_summary.json"),    emit: summary
    path("qc_alerts.txt"),      emit: alerts

    script:
    """
    #!/usr/bin/env python3
    import json, glob, os

    reports = glob.glob('*_qc_report.json')
    warnings_files = glob.glob('*_qc_warnings.txt')

    summary = {'total_studies': len(reports), 'studies': []}
    alerts_text = []
    has_issues = False

    for rpath in sorted(reports):
        with open(rpath) as f:
            r = json.load(f)
        summary['studies'].append({
            'study_id': r['study_id'],
            'level': r['overall_level'],
            'n_samples': r['n_samples'],
            'median_pct_filtered': r['median_pct_filtered'],
            'median_pct_nonchimeric': r['median_pct_nonchimeric'],
        })
        if r['overall_level'] in ('WARNING', 'ERROR', 'CRITICAL'):
            has_issues = True

    with open('qc_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    with open('qc_alerts.txt', 'w') as f:
        if has_issues:
            f.write('DENOISE QC SUMMARY: issues detected\\n')
            f.write('=' * 60 + '\\n\\n')
            for wpath in sorted(warnings_files):
                with open(wpath) as wf:
                    content = wf.read().strip()
                    if 'ALERT' in content or 'CRITICAL' in content or 'ERROR' in content:
                        f.write(content + '\\n\\n')
        else:
            f.write('All studies passed QC checks.\\n')
    """
}
