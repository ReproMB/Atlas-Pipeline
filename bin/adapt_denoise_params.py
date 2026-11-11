#!/usr/bin/env python3
"""
adapt_denoise_params.py — Compute second-pass DADA2 parameters for a study whose
first-pass denoise was catastrophic (filter/merge collapse or empty-table crash).

Decides, from the quality_stats.json written by parse_quality_stats.py, the failure mode
+ region/amplicon + the observed first-pass merged-ASV length, EITHER a better
paired-end parameter set OR — only when merging is physically impossible — a
single-end (forward-only) fallback.

Writes a small JSON the subworkflow turns into `meta.denoise_override` and uses
to route the retry to paired vs single denoise:
    {"mode": "paired"|"single", "trunc_len_f": N, "trunc_len_r": N,
     "min_overlap": N, "max_ee": float|null, "trunc_q": int|null, "reason": "..."}

Usage:
    adapt_denoise_params.py --quality-json q.json --failure-mode MERGE_COLLAPSE \
        --paired --region V3-V4 --output adapted.json
"""

import argparse
import json
import sys


def main():
    p = argparse.ArgumentParser(description='Adapt denoise params for a 2nd pass')
    p.add_argument('--quality-json', required=True)
    p.add_argument('--failure-mode', required=True,
                   choices=['FILTER_COLLAPSE', 'MERGE_COLLAPSE', 'ZERO_OUTPUT', 'EMPTY'])
    p.add_argument('--paired', action='store_true', help='Study was paired-end')
    p.add_argument('--region', default='')
    p.add_argument('--min-overlap', type=int, default=12)
    p.add_argument('--output', required=True)
    args = p.parse_args()

    with open(args.quality_json) as f:
        q = json.load(f)

    trim_f = q.get('trim_left_f', 0) or 0
    trim_r = q.get('trim_left_r', 0) or 0
    # Ceilings computed by parse_quality_stats.py (longest good-quality length, capped by
    # the read-length distribution). Fall back to read length, then trunc.
    maxq_f = q.get('max_qual_len_f') or q.get('cap_f') or q.get('readlen_f') or q.get('trunc_len_f', 0)
    maxq_r = q.get('max_qual_len_r') or q.get('cap_r') or q.get('readlen_r') or q.get('trunc_len_r', 0)
    cap_f = q.get('cap_f') or q.get('readlen_f') or maxq_f

    out = {'mode': 'paired' if args.paired else 'single',
           'trunc_len_f': None, 'trunc_len_r': None,
           'min_overlap': None, 'max_ee': None, 'trunc_q': None, 'reason': ''}

    # ── FILTER_COLLAPSE / ZERO_OUTPUT (and single-end EMPTY) ───────────────
    # Almost always trunc longer than the reads / variable length → pull trunc
    # back to the read-length cap and relax max_ee so more reads pass.
    if args.failure_mode in ('FILTER_COLLAPSE', 'ZERO_OUTPUT') and not (
            args.paired and q.get('merge_feasibility') in ('marginal', 'infeasible')):
        out['trunc_len_f'] = int(cap_f)
        if args.paired:
            out['trunc_len_r'] = int(q.get('cap_r') or q.get('readlen_r') or maxq_r)
        out['max_ee'] = 3.0
        out['reason'] = (f"{args.failure_mode}: pulled trunc to read-length cap "
                         f"and relaxed max_ee→3 to recover reads.")
        write(out, args)
        return

    # ── MERGE_COLLAPSE / EMPTY / merge-driven ZERO_OUTPUT ──────────────────
    # Decide paired (maximize coverage) vs single (physically can't span).
    # Feasibility uses the amplicon MINIMUM (as parse_quality_stats.py does): give the
    # paired retry one shot at max coverage if it can span the shortest amplicon;
    # otherwise forward-only. The observed merged length is NOT used here — for a
    # merge collapse it is biased short (only the short fragments merged). When
    # the region is unknown there is no reliable amplicon → single is the safe call.
    amp_min = q.get('amplicon_min')
    max_cov = (int(maxq_f) - trim_f) + (int(maxq_r) - trim_r) if args.paired else 0
    target = (amp_min + args.min_overlap) if amp_min else None

    infeasible = (
        not args.paired
        or q.get('merge_infeasible_hard')
        or amp_min is None
        or max_cov < target
    )

    if infeasible:
        out['mode'] = 'single'
        out['trunc_len_f'] = int(maxq_f) if maxq_f else int(cap_f)
        why = ('region unknown' if amp_min is None
               else f"max coverage {max_cov}bp < amplicon_min {amp_min} + overlap")
        out['reason'] = (f"{args.failure_mode}: merge not feasible ({why}) "
                         f"→ forward-only.")
    else:
        out['mode'] = 'paired'
        out['trunc_len_f'] = int(maxq_f)
        out['trunc_len_r'] = int(maxq_r)
        out['min_overlap'] = q.get('recommended_min_overlap') or args.min_overlap
        out['reason'] = (
            f"{args.failure_mode}: maximized coverage to {max_cov}bp "
            f"(trunc {maxq_f}/{maxq_r}) to span amplicon_min {amp_min}.")

    # A filter-driven failure routed here (paired study whose merge was already
    # marginal) still needs the filter relaxed; otherwise the retry can repeat
    # the first pass exactly when the rescue had already set trunc to maxq.
    if args.failure_mode in ('FILTER_COLLAPSE', 'ZERO_OUTPUT'):
        out['max_ee'] = 3.0
        out['reason'] += ' max_ee relaxed to 3 (filter-driven failure).'
    write(out, args)


def write(out, args):
    with open(args.output, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"  ADAPT [{args.failure_mode}] → mode={out['mode']} "
          f"trunc_f={out['trunc_len_f']} trunc_r={out['trunc_len_r']} "
          f"max_ee={out['max_ee']} :: {out['reason']}", file=sys.stderr)


if __name__ == '__main__':
    main()
