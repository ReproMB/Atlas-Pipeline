#!/usr/bin/env python3
"""
consensus.py — Consensus of repeated Lactobacillaceae species-group assignments.

DADA2 assignTaxonomy draws its bootstrap subsamples at random and
assign_taxonomy_gtdb.R does not fix a seed, so an ASV close to the minBoot
threshold can change group between runs. For release 2026.11 the species step
was run 10 times per processing unit and the runs were combined per ASV:

  * same group in every run                      → that group
  * group changes, but always within one GTDB genus
    (or between a species group and Unclassified) → unresolved (empty group)
  * GTDB genus changes between runs              → genus_change (empty group)
    (not possible with the release training set, which holds only
    Lactobacillus; kept for training sets with several genera)

An empty group means the ASV keeps its SILVA genus label in the resolved table
(export_resolved_table.py), exactly as an ASV below minBoot does.

Input:  the <prefix>_lacto_asv_classification.tsv written by each run of
        subclassify_lactobacillus.py.
Output: a classification TSV in the same format (ASV_ID, GTDB_species_group, …)
        plus the per-ASV category, the median genus bootstrap and the group
        assigned in each run. Pass it to export_resolved_table.py with
        --lacto-classification.

Usage:
    consensus.py runs/run_*/PREFIX_lacto_asv_classification.tsv \\
        --output taxonomy/species/PREFIX_lacto_asv_classification_consensus.tsv
"""

import argparse
import sys

import pandas as pd

UNCLASSIFIED = 'Unclassified'


def gtdb_genus(group):
    """'Lactobacillus crispatus' → 'Lactobacillus'; None for Unclassified/empty."""
    if not isinstance(group, str) or not group or group == UNCLASSIFIED:
        return None
    return group.split()[0]


def consensus(runs):
    """runs: list of DataFrames indexed by ASV_ID. Returns the consensus table."""
    ids = runs[0].index
    groups = pd.DataFrame({
        f'run_{k + 1:02d}': r.reindex(ids)['GTDB_species_group'].fillna(UNCLASSIFIED)
        for k, r in enumerate(runs)
    })
    boot = pd.DataFrame({
        f'run_{k + 1:02d}': r.reindex(ids)['boot_Genus'] for k, r in enumerate(runs)
    })

    keep = [c for c in ('SILVA_genus', 'mean_count', 'max_count') if c in runs[0].columns]
    out = runs[0][keep].copy()

    group, category = [], []
    for asv in ids:
        values = list(groups.loc[asv])
        distinct = set(values)
        if len(distinct) == 1:
            category.append('stable')
            group.append('' if values[0] == UNCLASSIFIED else values[0])
        elif len({gtdb_genus(v) for v in distinct} - {None}) > 1:
            category.append('genus_change')
            group.append('')
        else:
            category.append('unresolved')
            group.append('')

    out.insert(0, 'GTDB_species_group', group)
    out['category'] = category
    out['boot_Genus_median'] = boot.median(axis=1)
    out['groups_per_run'] = groups.apply(lambda r: ';'.join(r.values), axis=1)
    out.index.name = 'ASV_ID'
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('runs', nargs='+', help='lacto_asv_classification.tsv of each run')
    p.add_argument('--output', required=True)
    args = p.parse_args()

    if len(args.runs) < 2:
        sys.exit('ERROR: a consensus needs at least two runs')

    runs = [pd.read_csv(f, sep='\t', index_col=0) for f in sorted(args.runs)]
    ids = set(runs[0].index)
    if any(set(r.index) != ids for r in runs[1:]):
        sys.exit('ERROR: the runs do not classify the same ASVs (different inputs?)')

    out = consensus(runs)
    out.to_csv(args.output, sep='\t')
    n = out['category'].value_counts()
    print(f"{len(out)} Lactobacillaceae ASVs from {len(runs)} runs: "
          f"{n.get('stable', 0)} stable, {n.get('unresolved', 0)} unresolved, "
          f"{n.get('genus_change', 0)} genus change", file=sys.stderr)


if __name__ == '__main__':
    main()
