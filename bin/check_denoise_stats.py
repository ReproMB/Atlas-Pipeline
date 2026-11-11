#!/usr/bin/env python3
"""
check_denoise_stats.py — Evaluate DADA2 denoising statistics and generate QC alerts.

Reads the exported stats.tsv from DADA2 and evaluates:
  1. % reads passing quality filter
  2. % reads merged (paired-end only)
  3. % reads non-chimeric
  4. Samples with 0 output reads
  5. Read count heterogeneity (CV)
  6. Suitability for downstream rarefaction

Alert levels:
  INFO     — All metrics within healthy range
  WARNING  — Suboptimal but usable
  ERROR    — Processing likely failed; manual review needed
  CRITICAL — Samples with zero output

Usage:
    check_denoise_stats.py \
        --stats stats.tsv \
        --study-id PRJNA424567_Vagina \
        --min-passing-pct 50 \
        --min-merged-pct 60 \
        --min-nonchimeric-pct 60 \
        --min-reads-output 100 \
        --output-report qc_report.json \
        --output-warnings qc_warnings.txt
"""

import argparse
import csv
import json
import math
import sys


def read_stats(stats_path: str) -> list:
    """Read DADA2 denoising stats TSV."""
    samples = []
    with open(stats_path, 'r') as f:
        # Skip comment lines
        lines = [l for l in f if not l.startswith('#')]

    reader = csv.DictReader(lines, delimiter='\t')
    for row in reader:
        try:
            sample = {
                'sample_id': row.get('sample-id', row.get('', '')),
                'input': int(float(row.get('input', 0))),
                'filtered': int(float(row.get('filtered', 0))),
                'denoised': int(float(row.get('denoised', 0))),
                'merged': int(float(row.get('merged', 0))) if 'merged' in row else None,
                'non_chimeric': int(float(row.get('non-chimeric', 0))),
            }
            # Compute percentages
            inp = sample['input']
            if inp > 0:
                sample['pct_filtered'] = 100.0 * sample['filtered'] / inp
                sample['pct_nonchimeric'] = 100.0 * sample['non_chimeric'] / inp
                if sample['merged'] is not None:
                    sample['pct_merged'] = 100.0 * sample['merged'] / inp
                else:
                    sample['pct_merged'] = None
            else:
                sample['pct_filtered'] = 0.0
                sample['pct_nonchimeric'] = 0.0
                sample['pct_merged'] = 0.0 if sample['merged'] is not None else None

            samples.append(sample)
        except (ValueError, KeyError) as e:
            print(f"WARNING: Could not parse row: {e}", file=sys.stderr)
            continue

    return samples


def evaluate_stats(samples: list, min_passing_pct: float, min_merged_pct: float,
                   min_nonchimeric_pct: float, min_reads_output: int,
                   study_id: str, max_cv: float = 2.0) -> dict:
    """Evaluate stats and generate QC report with alerts."""
    alerts = []
    overall_level = 'INFO'
    is_paired = any(s['merged'] is not None for s in samples)

    # ── Per-sample alerts ────────────────────────────────────────────────
    zero_output = []
    low_output = []
    low_filter = []
    low_merge = []
    low_chimera = []

    for s in samples:
        sid = s['sample_id']

        # Zero output check
        if s['non_chimeric'] == 0:
            zero_output.append(sid)
            continue

        # Low absolute output
        if s['non_chimeric'] < min_reads_output:
            low_output.append(sid)

        # Filter pass rate
        if s['pct_filtered'] < min_passing_pct:
            low_filter.append(sid)

        # Merge rate (paired only)
        if is_paired and s['pct_merged'] is not None:
            if s['pct_merged'] < min_merged_pct:
                low_merge.append(sid)

        # Non-chimeric rate
        if s['pct_nonchimeric'] < min_nonchimeric_pct:
            low_chimera.append(sid)

    # ── Aggregate alerts ─────────────────────────────────────────────────
    n = len(samples)

    if zero_output:
        alerts.append({
            'level': 'CRITICAL',
            'message': f"{len(zero_output)}/{n} samples have 0 reads output",
            'samples': zero_output[:20],  # Limit to first 20
        })
        overall_level = 'CRITICAL'

    if low_output:
        alerts.append({
            'level': 'ERROR',
            'message': f"{len(low_output)}/{n} samples have < {min_reads_output} reads output",
            'samples': low_output[:20],
        })
        if overall_level not in ('CRITICAL',):
            overall_level = 'ERROR'

    if low_filter:
        pct_affected = 100.0 * len(low_filter) / n
        level = 'ERROR' if pct_affected > 50 else 'WARNING'
        alerts.append({
            'level': level,
            'message': f"{len(low_filter)}/{n} samples have < {min_passing_pct}% reads passing filter",
            'samples': low_filter[:20],
        })
        if level == 'ERROR' and overall_level not in ('CRITICAL',):
            overall_level = 'ERROR'
        elif level == 'WARNING' and overall_level == 'INFO':
            overall_level = 'WARNING'

    if low_merge and is_paired:
        pct_affected = 100.0 * len(low_merge) / n
        level = 'ERROR' if pct_affected > 50 else 'WARNING'
        alerts.append({
            'level': level,
            'message': f"{len(low_merge)}/{n} samples have < {min_merged_pct}% reads merged",
            'samples': low_merge[:20],
        })
        if level == 'ERROR' and overall_level not in ('CRITICAL',):
            overall_level = 'ERROR'
        elif level == 'WARNING' and overall_level == 'INFO':
            overall_level = 'WARNING'

    if low_chimera:
        pct_affected = 100.0 * len(low_chimera) / n
        level = 'ERROR' if pct_affected > 50 else 'WARNING'
        alerts.append({
            'level': level,
            'message': f"{len(low_chimera)}/{n} samples have < {min_nonchimeric_pct}% non-chimeric reads",
            'samples': low_chimera[:20],
        })
        if level == 'ERROR' and overall_level not in ('CRITICAL',):
            overall_level = 'ERROR'
        elif level == 'WARNING' and overall_level == 'INFO':
            overall_level = 'WARNING'

    # ── Heterogeneity check ──────────────────────────────────────────────
    output_counts = [s['non_chimeric'] for s in samples if s['non_chimeric'] > 0]
    if output_counts:
        mean_out = sum(output_counts) / len(output_counts)
        if mean_out > 0:
            variance = sum((c - mean_out) ** 2 for c in output_counts) / len(output_counts)
            cv = math.sqrt(variance) / mean_out
            if cv > max_cv:
                alerts.append({
                    'level': 'WARNING',
                    'message': f"High heterogeneity in output read counts (CV={cv:.2f})",
                    'samples': [],
                })
                if overall_level == 'INFO':
                    overall_level = 'WARNING'

    # ── Global summary stats ─────────────────────────────────────────────
    median_pct_filter = sorted([s['pct_filtered'] for s in samples])[n // 2] if n else 0
    median_pct_nonchim = sorted([s['pct_nonchimeric'] for s in samples])[n // 2] if n else 0
    median_pct_merge = None
    if is_paired:
        merge_vals = [s['pct_merged'] for s in samples if s['pct_merged'] is not None]
        median_pct_merge = sorted(merge_vals)[len(merge_vals) // 2] if merge_vals else 0

    report = {
        'study_id': study_id,
        'n_samples': n,
        'is_paired': is_paired,
        'overall_level': overall_level,
        'median_pct_filtered': round(median_pct_filter, 2),
        'median_pct_merged': round(median_pct_merge, 2) if median_pct_merge is not None else None,
        'median_pct_nonchimeric': round(median_pct_nonchim, 2),
        'zero_output_samples': len(zero_output),
        'low_output_samples': len(low_output),
        'alerts': alerts,
    }

    return report


def main():
    parser = argparse.ArgumentParser(description='Check DADA2 denoising stats')
    parser.add_argument('--stats', required=True, help='DADA2 stats TSV')
    parser.add_argument('--study-id', required=True, help='Study identifier')
    parser.add_argument('--quality-json', default=None, help='Quality stats JSON (for overlap warnings)')
    parser.add_argument('--min-passing-pct', type=float, default=50.0)
    parser.add_argument('--min-merged-pct', type=float, default=60.0)
    parser.add_argument('--min-nonchimeric-pct', type=float, default=60.0)
    parser.add_argument('--min-reads-output', type=int, default=100)
    parser.add_argument('--max-cv', type=float, default=2.0,
                        help='CV of output read counts above which a WARNING is raised')
    # Catastrophic-failure thresholds (trigger the adaptive denoise retry).
    parser.add_argument('--retry-filter-pct', type=float, default=40.0)
    parser.add_argument('--retry-merged-pct', type=float, default=30.0)
    parser.add_argument('--retry-zero-frac', type=float, default=0.30)
    parser.add_argument('--output-report', required=True, help='Output JSON report')
    parser.add_argument('--output-warnings', required=True, help='Output warnings TXT')
    args = parser.parse_args()

    samples = read_stats(args.stats)
    if not samples:
        print(f"ERROR: No samples found in {args.stats}", file=sys.stderr)
        sys.exit(1)

    report = evaluate_stats(
        samples,
        args.min_passing_pct,
        args.min_merged_pct,
        args.min_nonchimeric_pct,
        args.min_reads_output,
        args.study_id,
        max_cv=args.max_cv,
    )

    # ── Inject overlap warnings from quality stats JSON ────────────────
    overlap_warning = None
    min_overlap_used = None
    if args.quality_json:
        try:
            with open(args.quality_json) as f:
                qstats = json.load(f)
            min_overlap_used = qstats.get('recommended_min_overlap')
            # Find overlap warnings in the quality stats warnings list
            for w in qstats.get('warnings', []):
                if 'OVERLAP' in w:
                    overlap_warning = w
                    report['alerts'].append({
                        'level': 'WARNING' if 'WARNING' in w else 'CRITICAL',
                        'message': w,
                        'samples': [],
                    })
                    if report['overall_level'] == 'INFO':
                        report['overall_level'] = 'WARNING'
                    break
            if min_overlap_used and min_overlap_used != 12:
                report['min_overlap_adjusted'] = min_overlap_used

            # ── Correctly attribute a low merge rate ────────────────────
            # When the preflight predicted insufficient overlap AND merging
            # actually came out low, make the cause explicit so it is not
            # mistaken for a read-quality problem.
            merge_feasibility = qstats.get('merge_feasibility')
            report['merge_feasibility'] = merge_feasibility
            median_merged = report.get('median_pct_merged')
            low_merge_observed = (
                report['is_paired'] and median_merged is not None
                and median_merged < args.min_merged_pct
            )
            if merge_feasibility in ('infeasible', 'marginal') and low_merge_observed:
                level = 'ERROR' if merge_feasibility == 'infeasible' else 'WARNING'
                report['alerts'].append({
                    'level': level,
                    'message': (
                        f"MERGE_INFEASIBLE: low merge rate (median {median_merged}%) is "
                        f"caused by insufficient read overlap "
                        f"(merge_feasibility={merge_feasibility}), not read quality. "
                        f"Increase truncation lengths or process forward reads only."
                    ),
                    'samples': [],
                })
                if level == 'ERROR' and report['overall_level'] not in ('CRITICAL',):
                    report['overall_level'] = 'ERROR'
                elif level == 'WARNING' and report['overall_level'] == 'INFO':
                    report['overall_level'] = 'WARNING'
        except Exception as e:
            print(f"WARNING: Could not read quality JSON: {e}", file=sys.stderr)

    # ── Classify catastrophic failure (triggers the adaptive retry) ────────
    # Cause-ordered so the retry knows what to fix. Low non-chimeric alone is
    # NOT catastrophic (chimeras are usually real); only filter/merge collapse
    # or mass zero-output qualify.
    n = report['n_samples']
    zero_frac = (report['zero_output_samples'] / n) if n else 1.0
    med_filt = report['median_pct_filtered']
    med_merged = report['median_pct_merged']
    catastrophic = False
    failure_mode = None
    if med_filt is not None and med_filt < args.retry_filter_pct:
        catastrophic, failure_mode = True, 'FILTER_COLLAPSE'
    elif report['is_paired'] and med_merged is not None and med_merged < args.retry_merged_pct:
        catastrophic, failure_mode = True, 'MERGE_COLLAPSE'
    elif zero_frac > args.retry_zero_frac:
        catastrophic, failure_mode = True, 'ZERO_OUTPUT'
    report['catastrophic'] = catastrophic
    report['failure_mode'] = failure_mode
    if catastrophic:
        print(f"  CATASTROPHIC: {args.study_id} failure_mode={failure_mode} "
              f"(filt={med_filt}, merged={med_merged}, zero_frac={zero_frac:.2f})",
              file=sys.stderr)

    # Write JSON report
    with open(args.output_report, 'w') as f:
        json.dump(report, f, indent=2)

    # Write human-readable warnings
    with open(args.output_warnings, 'w') as f:
        if report['overall_level'] in ('WARNING', 'ERROR', 'CRITICAL'):
            f.write(f"{'='*60}\n")
            f.write(f"DENOISE QC ALERT — {args.study_id}\n")
            f.write(f"Overall level: {report['overall_level']}\n")
            f.write(f"{'='*60}\n\n")
            f.write(f"Samples: {report['n_samples']}\n")
            f.write(f"Median % filtered: {report['median_pct_filtered']}\n")
            if report['median_pct_merged'] is not None:
                f.write(f"Median % merged: {report['median_pct_merged']}\n")
            f.write(f"Median % non-chimeric: {report['median_pct_nonchimeric']}\n")
            if min_overlap_used and min_overlap_used != 12:
                f.write(f"\nmin_overlap adjusted: 12 → {min_overlap_used}\n")
            f.write('\n')

            for alert in report['alerts']:
                f.write(f"[{alert['level']}] {alert['message']}\n")
                if alert['samples']:
                    f.write(f"  Affected: {', '.join(alert['samples'][:10])}")
                    if len(alert['samples']) > 10:
                        f.write(f" ... and {len(alert['samples'])-10} more")
                    f.write('\n')
                f.write('\n')

            f.write(f"{'='*60}\n")
            f.write("ACTION REQUIRED: Review quality profiles and denoise parameters.\n")
        else:
            f.write(f"Study {args.study_id}: All QC checks passed (INFO level).\n")
            if min_overlap_used and min_overlap_used != 12:
                f.write(f"Note: min_overlap was adjusted to {min_overlap_used} based on expected overlap.\n")

    # Print summary to stderr
    level = report['overall_level']
    print(f"{args.study_id}: {level} "
          f"(filter={report['median_pct_filtered']:.1f}%, "
          f"nonchim={report['median_pct_nonchimeric']:.1f}%)",
          file=sys.stderr)


if __name__ == '__main__':
    main()
