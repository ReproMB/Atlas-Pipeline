#!/usr/bin/env python3
"""
generate_manifest.py — Generate QIIME2 manifest files from fetchngs output.

Scans a data directory for FASTQ files and generates CSV manifests
compatible with QIIME2 import (PairedEnd or SingleEnd Phred33V2).

Usage:
    generate_manifest.py \
        --study PRJNA424567 \
        --sites Vagina,Gut \
        --data-dir /path/to/data \
        --paired \
        --output manifest.csv
"""

import argparse
import csv
import re
import sys
from pathlib import Path


# ── Regex patterns for read-number detection ──────────────────────────────
# These require an explicit delimiter (_  or .) before R1/R2 or 1/2,
# preventing greedy stripping of sample-name characters like "C100b2".
_RE_READ2 = re.compile(r'(?:[_.]R2|_2)(?:_\d+)?\.(?:fastq|fq)\.gz$')
_RE_STRIP_READNUM = re.compile(r'(?:[_.]R[12]|_[12])(?:_\d+)?\.(?:fastq|fq)\.gz$')
_RE_STRIP_EXT = re.compile(r'\.(?:fastq|fq)\.gz$')


def _extract_sample_id(filename: str) -> tuple:
    """
    Extract (base_name, read_number) from a FASTQ filename.

    Returns the sample base name and '1' or '2' for read number.
    """
    # Try to strip a read-number suffix first
    base_name = _RE_STRIP_READNUM.sub('', filename)
    if base_name == filename:
        # No read-number suffix found — strip extension only
        base_name = _RE_STRIP_EXT.sub('', filename)

    # Determine read number
    read_num = '2' if _RE_READ2.search(filename) else '1'

    return base_name, read_num


def find_fastq_files(data_dir: str, study_id: str, sites: list) -> dict:
    """
    Find FASTQ files for a study across given body sites.

    Searches in order:
        {data_dir}/{site}/{study_id}/fastq/
        {data_dir}/{site}/{study_id}/
        {data_dir}/{study_id}/fastq/
        {data_dir}/{study_id}/

    Returns dict: {sample_id: {'forward': path, 'reverse': path, 'site': site}}
    """
    samples = {}
    data_path = Path(data_dir)

    for site in sites:
        search_paths = [
            data_path / site.lower() / study_id / 'fastq',
            data_path / site.lower() / study_id,
            data_path / study_id / 'fastq',
            data_path / study_id,
        ]

        for search_path in search_paths:
            if not search_path.exists():
                continue

            fastq_files = sorted(
                list(search_path.glob('*.fastq.gz')) +
                list(search_path.glob('*.fq.gz'))
            )
            if not fastq_files:
                continue

            for fq in fastq_files:
                base_name, read_num = _extract_sample_id(fq.name)

                sample_id = base_name
                # Try to extract the Run ID specifically (e.g. from SRX..._SRR...)
                run_match = re.search(r'([SED]RR\d+)', base_name)
                if run_match:
                    sample_id = run_match.group(1)

                if sample_id not in samples:
                    samples[sample_id] = {
                        'forward': None, 'reverse': None, 'site': site
                    }

                abs_path = str(fq.resolve())
                if read_num == '2':
                    samples[sample_id]['reverse'] = abs_path
                else:
                    samples[sample_id]['forward'] = abs_path

            if samples:
                break  # found files in this path, skip fallbacks

    return samples


def write_manifest_paired(samples: dict, output_path: str, sample_filter=None):
    """Write PairedEndFastqManifestPhred33V2 CSV."""
    written = 0
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample-id', 'forward-absolute-filepath', 'reverse-absolute-filepath'])
        for sid in sorted(samples):
            if sample_filter and sid not in sample_filter:
                continue
            info = samples[sid]
            if info['forward'] and info['reverse']:
                writer.writerow([sid, info['forward'], info['reverse']])
                written += 1
            else:
                print(f"WARNING: {sid} missing read pair, skipping", file=sys.stderr)
    return written


def write_manifest_single(samples: dict, output_path: str, sample_filter=None):
    """Write SingleEndFastqManifestPhred33V2 CSV."""
    written = 0
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['sample-id', 'absolute-filepath'])
        for sid in sorted(samples):
            if sample_filter and sid not in sample_filter:
                continue
            info = samples[sid]
            if info['forward']:
                writer.writerow([sid, info['forward']])
                written += 1
            else:
                print(f"WARNING: {sid} has no forward read, skipping", file=sys.stderr)
    return written


def main():
    parser = argparse.ArgumentParser(description='Generate QIIME2 manifest')
    parser.add_argument('--study', required=True, help='Project accession')
    parser.add_argument('--sites', required=True, help='Comma-separated body sites')
    parser.add_argument('--data-dir', required=True, help='Base data directory')
    parser.add_argument('--paired', action='store_true', help='Paired-end manifest')
    parser.add_argument('--output', required=True, help='Output manifest CSV')
    parser.add_argument('--sample-filter', default=None,
                        help='File with sample IDs to include (one per line)')
    parser.add_argument('--mode-file', default=None,
                        help='Output file to record the actual mode used (PAIRED or SINGLE)')
    args = parser.parse_args()

    sites = [s.strip() for s in args.sites.split(',')]
    sample_filter = None
    if args.sample_filter:
        with open(args.sample_filter) as f:
            sample_filter = set(l.strip() for l in f if l.strip())

    samples = find_fastq_files(args.data_dir, args.study.strip(), sites)
    if not samples:
        print(f"ERROR: No FASTQ files found for {args.study} in {sites}", file=sys.stderr)
        sys.exit(1)

    # Determine actual mode: auto-fallback from paired to single if no pairs exist
    actual_mode = 'SINGLE'
    if args.paired:
        has_pairs = any(s['forward'] and s['reverse'] for s in samples.values())
        if has_pairs:
            actual_mode = 'PAIRED'
            n = write_manifest_paired(samples, args.output, sample_filter)
        else:
            print(f"WARNING: --paired specified but no paired reads found for {args.study}. "
                  f"Falling back to single-end mode.", file=sys.stderr)
            n = write_manifest_single(samples, args.output, sample_filter)
    else:
        n = write_manifest_single(samples, args.output, sample_filter)

    print(f"Manifest: {n} samples → {args.output} (mode={actual_mode})", file=sys.stderr)

    # Write mode file for downstream pipeline coordination
    if args.mode_file:
        with open(args.mode_file, 'w') as f:
            f.write(actual_mode)

    if n == 0:
        print(f"ERROR: Manifest is empty — no samples survived filtering for {args.study}. "
              f"Check that sample IDs in the filter file match the FASTQ filenames.",
              file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
