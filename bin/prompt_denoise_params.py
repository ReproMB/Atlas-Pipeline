#!/usr/bin/env python3
"""
prompt_denoise_params.py — Interactively collect per-study DADA2 denoise
parameters and write a manual_params.csv consumed by the pipeline's
`--denoise_mode manual`.

Because a Nextflow process cannot reliably prompt mid-run, this standalone
helper is run by the user BEFORE the pipeline. It asks, for each study unit in
the input CSV, for the denoise parameters to use instead of the automatic ones.

Study-unit identity (the `id` column) is produced exactly as main.nf builds
`meta.id`:
  1. clean_primers.py is run on the input CSV (same canonicalization the
     pipeline uses) to obtain ForwardPrimerCanonical / ReversePrimerCanonical /
     PrimerRegionID.
  2. Rows are grouped by (study, site_group, region, paired, fwd, rev,
     primer_region_id); reproductive sites of the same study collapse into one
     unit (id = study_<sorted sites>), gut sites use id = study_<first site>.

Any cell left blank falls back downstream to the automatic value (then global
config), so you only need to type the parameters you actually want to override.

Usage:
    prompt_denoise_params.py TEST.csv -o manual_params.csv
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

# Must stay in sync with main.nf (repro_sites). Single source of the grouping rule.
REPRO_SITES = ['Endometrium', 'Vagina', 'Cervix', 'Cervicovaginal']

# Fields collected per study, in prompt order, with (label, kind, hint).
# kind: 'int' | 'float'. Reverse/overlap fields are skipped for single-end.
FIELDS = [
    ('trunc_len_f', 'int',   'forward 3\' truncation length',     False),
    ('trunc_len_r', 'int',   'reverse 3\' truncation length',     True),
    ('trim_left_f', 'int',   'forward 5\' trim (e.g. primer len, 0 if none)', False),
    ('trim_left_r', 'int',   'reverse 5\' trim',                  True),
    ('min_overlap', 'int',   'min overlap for paired merge (default 12)', True),
    ('max_ee',      'float', 'max expected errors (default 2.0)', False),
    ('trunc_q',     'int',   'truncation quality threshold (default 2)', False),
]
OUT_COLUMNS = ['id'] + [f[0] for f in FIELDS]


def run_clean_primers(input_csv: str, script_dir: str) -> str:
    """Run clean_primers.py on input_csv, return path to a temp cleaned CSV."""
    clean_script = os.path.join(script_dir, 'clean_primers.py')
    if not os.path.exists(clean_script):
        sys.exit(f"ERROR: clean_primers.py not found next to this script ({clean_script})")
    fd, out_path = tempfile.mkstemp(suffix='_cleaned.csv', prefix='denoise_units_')
    os.close(fd)
    try:
        subprocess.run(
            [sys.executable, clean_script,
             '--input-csv', input_csv, '--output-csv', out_path],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        sys.exit(f"ERROR: clean_primers.py failed ({e})")
    return out_path


def build_study_units(cleaned_csv: str) -> list:
    """
    Replicate main.nf grouping → ordered list of study-unit dicts.
    Each: {id, study, site_group, region, paired, sites}
    """
    units = {}          # group_key -> unit dict
    order = []          # preserve first-seen order of group keys
    with open(cleaned_csv, newline='') as f:
        for row in csv.DictReader(f):
            study = (row.get('StudyProject') or '').strip()
            if not study:
                continue
            site = (row.get('Site') or '').strip()
            region = (row.get('Region') or '').strip()
            paired = (row.get('Paired/Single') or '').strip()
            fwd = (row.get('ForwardPrimerCanonical') or row.get('ForwardPrimerClean') or '').strip()
            rev = (row.get('ReversePrimerCanonical') or row.get('ReversePrimerClean') or '').strip()
            prid = (row.get('PrimerRegionID') or '').strip()
            site_group = 'reproductive' if site in REPRO_SITES else 'gut'

            key = (study, site_group, region, paired, fwd, rev, prid)
            if key not in units:
                units[key] = {
                    'study': study, 'site_group': site_group,
                    'region': region, 'paired': paired, 'sites': [],
                }
                order.append(key)
            # Collect sites without deduplication, mirroring main.nf groupTuple
            # (which keeps every row's site, in CSV order).
            units[key]['sites'].append(site)

    result = []
    for key in order:
        u = units[key]
        sites = u['sites'] or ['']
        if u['site_group'] == 'reproductive':
            site_part = '_'.join(sorted(sites))
        else:
            site_part = sites[0]
        u['id'] = f"{u['study']}_{site_part}"
        result.append(u)
    return result


def ask(field_name: str, kind: str, hint: str) -> str:
    """Prompt until a blank or a valid int/float is entered. Returns the raw string."""
    while True:
        raw = input(f"    {field_name:<12} [{hint}] (Enter=auto): ").strip()
        if raw == '':
            return ''
        try:
            if kind == 'int':
                int(raw)
            else:
                float(raw)
            return raw
        except ValueError:
            print(f"      ! not a valid {kind}, try again", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description='Collect per-study DADA2 denoise params')
    parser.add_argument('input_csv', help='Pipeline input study table CSV (e.g. TEST.csv)')
    parser.add_argument('-o', '--output', default='manual_params.csv',
                        help='Output manual params CSV (default: manual_params.csv)')
    args = parser.parse_args()

    if not os.path.exists(args.input_csv):
        sys.exit(f"ERROR: input CSV not found: {args.input_csv}")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    cleaned = run_clean_primers(args.input_csv, script_dir)
    try:
        units = build_study_units(cleaned)
    finally:
        try:
            os.remove(cleaned)
        except OSError:
            pass

    if not units:
        sys.exit("ERROR: no study units found in the input CSV")

    print(f"\nFound {len(units)} study unit(s). For each, enter denoise params "
          f"(blank = use automatic value).\n", file=sys.stderr)

    rows = []
    for i, u in enumerate(units, 1):
        is_paired = u['paired'].upper() == 'PAIRED'
        print(f"[{i}/{len(units)}] {u['id']}  "
              f"(region={u['region'] or '?'}, {u['paired'] or '?'})", file=sys.stderr)
        row = {'id': u['id']}
        for name, kind, hint, reverse_only in FIELDS:
            if reverse_only and not is_paired:
                row[name] = ''
                continue
            row[name] = ask(name, kind, hint)
        rows.append(row)
        print('', file=sys.stderr)

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} study unit(s) → {args.output}", file=sys.stderr)
    print(f"Run the pipeline with: "
          f"--denoise_mode manual --manual_params {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
