#!/usr/bin/env python3
"""
clean_primers.py — Strip adapter/barcode sequences from primers
and deduplicate primer combinations using IUPAC generalization.

Usage:
    clean_primers.py --input-csv studies.csv --output-csv cleaned_primers.csv
"""

import argparse
import csv
import re
import sys
from collections import defaultdict

# ── IUPAC ambiguity codes ────────────────────────────────────────────────────
IUPAC = {
    'A': frozenset('A'),   'C': frozenset('C'),
    'G': frozenset('G'),   'T': frozenset('T'),
    'R': frozenset('AG'),  'Y': frozenset('CT'),
    'S': frozenset('GC'),  'W': frozenset('AT'),
    'K': frozenset('GT'),  'M': frozenset('AC'),
    'B': frozenset('CGT'), 'D': frozenset('AGT'),
    'H': frozenset('ACT'), 'V': frozenset('ACG'),
    'N': frozenset('ACGT'),
}

# ── Known adapter sequences (sorted longest first for greedy match) ──────────
ADAPTERS = sorted([
    'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG',    # Nextera Read1
    'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG',   # Nextera Read2
    'AATGATACGGCGACCACCGAGATCTACAC',         # P5 adapter
    'CAAGCAGAAGACGGCATACGAGAT',              # P7 adapter
    'CCTATCCCCTGTGTGCCTTGGCAGTCTCAG',        # 454 Adapter A
    'CCATCTCATCCCTGCGTGTCTCCGACTCAG',        # 454 Adapter B
    'GCCTTGCCAGCCCGCTCAGTC',                 # Ion Torrent A
    'GCCTCCCTCGCGCCATCAG',                   # Ion Torrent trP1
], key=len, reverse=True)


def strip_adapters(primer: str) -> str:
    """Remove known adapter sequences and inline barcodes from a primer."""
    p = primer.upper().strip()
    # Strip known adapters from the 5' end
    for adapter in ADAPTERS:
        if p.startswith(adapter):
            p = p[len(adapter):]
            break
    # Strip inline barcode placeholders (X or N runs at the start)
    p = re.sub(r'^X{2,}', '', p)
    p = re.sub(r'^N{2,}', '', p)
    # Strip trailing artifacts like "-3" notation
    p = re.sub(r'-\d+$', '', p)
    return p


def generalizes(a: str, b: str) -> bool:
    """
    Check if primer 'a' is an IUPAC generalization of primer 'b'.
    Returns True if 'a' is equal to or more degenerate than 'b' at every position.
    Primers must be of equal length.
    """
    if len(a) != len(b):
        return False
    for ca, cb in zip(a.upper(), b.upper()):
        bases_a = IUPAC.get(ca)
        bases_b = IUPAC.get(cb)
        if bases_a is None or bases_b is None:
            return False
        if not bases_a.issuperset(bases_b):
            return False
    return True


def find_canonical_primer(primers: list) -> tuple:
    """
    Given a list of (fwd, rev) pairs of the same length, find the one that
    generalizes all others (or return the most general one).
    """
    if len(primers) == 1:
        return primers[0]

    best = primers[0]
    for candidate in primers[1:]:
        fwd_c, rev_c = candidate
        fwd_b, rev_b = best
        if generalizes(fwd_c, fwd_b) and generalizes(rev_c, rev_b):
            best = candidate
    return best


def deduplicate_primers(primer_pairs: list) -> dict:
    """
    Group primer pairs by IUPAC generalization and return a mapping
    from each original pair to its canonical (most general) representative.

    Args:
        primer_pairs: list of (fwd, rev) tuples

    Returns:
        dict mapping (fwd, rev) -> (canonical_fwd, canonical_rev)
    """
    unique = list(set(primer_pairs))
    used = set()
    groups = []  # list of [member_indices]

    for i in range(len(unique)):
        if i in used:
            continue
        group = [i]
        f1, r1 = unique[i]
        for j in range(i + 1, len(unique)):
            if j in used:
                continue
            f2, r2 = unique[j]
            if len(f1) == len(f2) and len(r1) == len(r2):
                if (generalizes(f1, f2) and generalizes(r1, r2)) or \
                   (generalizes(f2, f1) and generalizes(r2, r1)):
                    group.append(j)
                    used.add(j)
        used.add(i)
        groups.append(group)

    # Build mapping
    mapping = {}
    for group in groups:
        members = [unique[i] for i in group]
        canonical = find_canonical_primer(members)
        for member in members:
            mapping[member] = canonical

    return mapping


def generate_region_id(fwd: str, rev: str, region: str) -> str:
    """Generate a stable identifier for a primer+region combination."""
    # Use a short hash of the primers to create a unique but readable ID
    import hashlib
    primer_hash = hashlib.md5(f"{fwd}_{rev}".encode()).hexdigest()[:8]
    safe_region = region.replace('-', '_')
    return f"{safe_region}_{primer_hash}"


def main():
    parser = argparse.ArgumentParser(description='Clean and deduplicate 16S primers')
    parser.add_argument('--input-csv', required=True, help='Input study table CSV')
    parser.add_argument('--output-csv', required=True, help='Output cleaned CSV')
    parser.add_argument('--dedup-report', default=None, help='Optional dedup report file')
    args = parser.parse_args()

    # ── Read and clean ───────────────────────────────────────────────────────
    rows = []
    with open(args.input_csv, 'r') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        for row in reader:
            if not row.get('StudyProject') or row.get('StudyProject').strip() == '':
                continue
            fwd = row.get('ForwardPrimer', '') or ''
            rev = row.get('ReversePrimer', '') or ''
            row['ForwardPrimerClean'] = strip_adapters(fwd)
            row['ReversePrimerClean'] = strip_adapters(rev)
            rows.append(row)

    # ── Deduplicate per region ───────────────────────────────────────────────
    primers_by_region = defaultdict(list)
    for row in rows:
        key = (row['ForwardPrimerClean'], row['ReversePrimerClean'])
        region = (row.get('Region', '') or '').strip()
        primers_by_region[region].append(key)

    all_mappings = {}
    report_lines = []
    for region, pairs in sorted(primers_by_region.items()):
        mapping = deduplicate_primers(pairs)
        for orig, canon in mapping.items():
            all_mappings[(region, orig)] = canon
        # Report
        unique_canons = set(mapping.values())
        unique_origs = set(mapping.keys())
        report_lines.append(
            f"{region}: {len(unique_origs)} unique → {len(unique_canons)} classifier(s)"
        )
        if len(unique_origs) > len(unique_canons):
            for orig, canon in mapping.items():
                if orig != canon:
                    report_lines.append(f"  MERGED: {orig[0][:30]}... → {canon[0][:30]}...")

    # ── Assign canonical primers and region IDs ──────────────────────────────
    out_fieldnames = fieldnames + [
        'ForwardPrimerClean', 'ReversePrimerClean',
        'ForwardPrimerCanonical', 'ReversePrimerCanonical',
        'PrimerRegionID'
    ]
    # Remove duplicates in fieldnames (ForwardPrimerClean already added)
    seen = set()
    unique_fieldnames = []
    for fn in out_fieldnames:
        if fn not in seen:
            unique_fieldnames.append(fn)
            seen.add(fn)

    for row in rows:
        region = (row.get('Region', '') or '').strip()
        key = (row['ForwardPrimerClean'], row['ReversePrimerClean'])
        canon = all_mappings.get((region, key), key)
        row['ForwardPrimerCanonical'] = canon[0]
        row['ReversePrimerCanonical'] = canon[1]
        row['PrimerRegionID'] = generate_region_id(canon[0], canon[1], region)

    # ── Write output ─────────────────────────────────────────────────────────
    with open(args.output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=unique_fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

    # ── Write dedup report ───────────────────────────────────────────────────
    if args.dedup_report:
        with open(args.dedup_report, 'w') as f:
            f.write("Primer Deduplication Report\n")
            f.write("=" * 40 + "\n")
            for line in report_lines:
                f.write(line + "\n")
            total_unique = len(set(
                (r['ForwardPrimerCanonical'], r['ReversePrimerCanonical'], (r.get('Region', '') or '').strip())
                for r in rows if r['ForwardPrimerCanonical'] and r['ReversePrimerCanonical']
            ))
            f.write(f"\nTotal classifiers needed: {total_unique}\n")

    print(f"Cleaned {len(rows)} entries → {args.output_csv}", file=sys.stderr)


if __name__ == '__main__':
    main()
