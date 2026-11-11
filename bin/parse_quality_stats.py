#!/usr/bin/env python3
"""
parse_quality_stats.py — Extract quality statistics from QIIME2 demux summarize QZV
and auto-determine DADA2 denoise parameters (trunc/trim).

The QZV is a ZIP file containing quality summary TSVs:
  - forward-seven-number-summaries.tsv
  - reverse-seven-number-summaries.tsv (paired-end only)

Algorithm:
  1. Extract per-position quality at multiple percentiles (25%, 9%).
  2. PRIMARY check  — sliding-window average of the 25th percentile drops
     below --min-quality.  When triggered the cut point is walked LEFT on
     raw values to find the last base still ≥ threshold.
  3. SECONDARY (tail) check — even when the primary check finds no truncation,
     scan the last --tail-check-window positions (default 20) for collapse of
     the 9th percentile below --min-quality-p9 (default 15).  If the collapse
     is found, the truncation point is set to the START of that collapse —
     capturing the 'looks-ok-on-median-but-9th-pct-is-junk' pattern.
  4. trim_left = primer length (if primers not already stripped).
  5. Emit JSON with recommended params + warnings.

Usage:
    parse_quality_stats.py \
        --qzv summary.qzv \
        --output quality_stats.json \
        --min-quality 20 \
        --min-quality-p9 15 \
        --primer-len-f 17 \
        --primer-len-r 21 \
        --tail-check-window 20
"""

import argparse
import csv
import json
import os
import sys
import tempfile
import zipfile


# Expected amplicon sizes (without primers) by 16S region.
# Values are (min_bp, max_bp).
AMPLICON_SIZES = {
    'V1-V2': (290, 330),
    'V1-V3': (450, 510),
    'V3-V4': (420, 460),
    'V3-V5': (550, 590),
    'V4':    (250, 260),
    'V4-V5': (350, 390),
    'V5-V6': (260, 300),
    'V6-V8': (400, 450),
}

# Auto-truncation is never allowed to exceed the read length that still retains
# this fraction of reads. DADA2 discards every read shorter than --p-trunc-len,
# so a quality-based truncation longer than the reads themselves (common with
# variable-length / pyrosequencing data) silently wipes whole samples. This cap
# clamps trunc to the actual read-length distribution (the per-position read
# counts already present in the QIIME2 quality table). Not a Nextflow param.
READLEN_RETENTION = 0.8


def extract_quality_table(qzv_path: str, direction: str = 'forward') -> list:
    """
    Extract per-position quality summaries from a QZV file.

    QIIME2 seven-number-summaries format (transposed):
        Rows = percentile labels (2%, 9%, 25%, 50%, 75%, 91%, 98%)
        Columns = base positions (0, 1, 2, ...)

    Returns:
        List of dicts, one per position, each with keys like '50%', '25%', etc.
    """
    # Try both .tsv and .csv variants
    candidates = [
        f'{direction}-seven-number-summaries.tsv',
        f'{direction}-seven-number-summaries.csv',
    ]
    data = []

    with zipfile.ZipFile(qzv_path, 'r') as zf:
        matching = []
        for candidate in candidates:
            matching = [n for n in zf.namelist() if n.endswith(candidate)]
            if matching:
                break

        if not matching:
            print(f"  WARNING: No {direction} quality summary found in QZV",
                  file=sys.stderr)
            # List what IS in the QZV for debugging
            quality_files = [n for n in zf.namelist()
                           if 'seven-number' in n or 'quality' in n.lower()]
            if quality_files:
                print(f"  Available quality files: {quality_files}", file=sys.stderr)
            return data

        with zf.open(matching[0]) as f:
            content = f.read().decode('utf-8').strip()

        # Auto-detect delimiter
        first_line = content.split('\n')[0]
        delimiter = '\t' if '\t' in first_line else ','

        lines = content.split('\n')
        if len(lines) < 2:
            print(f"  WARNING: Quality file has < 2 lines", file=sys.stderr)
            return data

        # Parse header — these are the position indices
        header = lines[0].split(delimiter)
        # First element is empty or a label, rest are position numbers
        positions = []
        for h in header[1:]:
            try:
                positions.append(int(float(h.strip())))
            except ValueError:
                positions.append(h.strip())

        # Initialize per-position dicts
        position_data = [{} for _ in positions]

        # Parse each row — row label is the percentile, values are per-position
        for line in lines[1:]:
            if not line.strip():
                continue
            parts = line.split(delimiter)
            if len(parts) < 2:
                continue

            percentile_label = parts[0].strip()  # e.g. "50%", "2%", etc.

            for j, val_str in enumerate(parts[1:]):
                if j >= len(positions):
                    break
                try:
                    position_data[j][percentile_label] = float(val_str.strip())
                except (ValueError, TypeError):
                    continue

        data = position_data

    # Diagnostic output
    if data:
        print(f"  Parsed {len(data)} positions for {direction}", file=sys.stderr)
        # Show first and last few medians
        medians = [d.get('50%', d.get('median', None)) for d in data]
        medians_valid = [m for m in medians if m is not None]
        if medians_valid:
            print(f"  Median Q at pos 0: {medians_valid[0]:.1f}, "
                  f"pos {len(medians_valid)-1}: {medians_valid[-1]:.1f}",
                  file=sys.stderr)
    else:
        print(f"  WARNING: No quality data parsed for {direction}", file=sys.stderr)

    return data


def extract_sample_counts(qzv_path: str) -> dict:
    """Extract per-sample sequence counts from the QZV."""
    counts = {}
    with zipfile.ZipFile(qzv_path, 'r') as zf:
        matching = [n for n in zf.namelist()
                    if n.endswith('per-sample-fastq-counts.tsv')
                    or n.endswith('per-sample-sequence-counts.tsv')]
        if not matching:
            return counts

        with zf.open(matching[0]) as f:
            content = f.read().decode('utf-8')
            reader = csv.DictReader(content.strip().split('\n'), delimiter='\t')
            for row in reader:
                try:
                    # Column names vary: 'Sequence count', 'count', etc.
                    sid = row.get('Sample name', row.get('sample-id', ''))
                    count_val = row.get('Sequence count',
                                row.get('count',
                                row.get('Sequences', '0')))
                    if sid:
                        counts[sid] = int(float(count_val))
                except (ValueError, KeyError):
                    continue

    return counts


def _extract_percentile_values(quality_data: list, percentile: str) -> list:
    """
    Helper: pull raw float values for *percentile* from the per-position dicts.
    Missing values are filled with 40.0 (optimistic) so downstream smoothing
    is not skewed by sparse positions.
    """
    keys = [percentile, percentile.replace('%', ''), f"{percentile}%"]
    values = []
    for entry in quality_data:
        val = None
        for k in keys:
            if k in entry:
                val = entry[k]
                break
        if val is None:
            # Fallback: try common synonyms
            for k in ['50%', 'median', '50']:
                if k in entry:
                    val = entry[k]
                    break
        values.append(val if val is not None else 40.0)
    return values


def _sliding_average(values: list, window: int) -> list:
    """Return causal sliding-window (trailing) average."""
    out = []
    for i in range(len(values)):
        start = max(0, i - window + 1)
        chunk = values[start:i + 1]
        out.append(sum(chunk) / len(chunk))
    return out


def find_trunc_position(
    quality_data: list,
    min_quality: float,
    primer_len: int = 0,
    percentile: str = '25%',
    window: int = 10,
    min_quality_p9: float = 15.0,
    tail_check_window: int = 20,
) -> int:
    """
    Find the optimal truncation position using a dual-percentile strategy.

    PRIMARY CHECK
    -------------
    Sliding-window average of *percentile* (default 25th pct) vs *min_quality*
    detects where quality degrades.  Once that crossing position is found, the
    algorithm walks LEFT through the raw (unsmoothed) per-position values to
    find the last base whose raw value is still ≥ min_quality.  That position
    becomes the actual cut point, so we never keep a base that is already below
    threshold, but we also never cut further than strictly necessary.

    SECONDARY / TAIL CHECK
    ----------------------
    Even when the primary check reports "no truncation needed", the 9th
    percentile at the read tail often collapses well below acceptable quality
    (visible as the lowest whisker/dot nearly reaching 0 in QIIME2 plots).
    This check scans the last *tail_check_window* positions (after the primer
    region) for the 9th percentile dropping below *min_quality_p9*.  The
    truncation is set to exactly that first collapse position — no buffer
    is applied here because the position itself is already the first bad base.

    Returns 0 only when both checks report no issues — which is now much harder
    to trigger for typical Illumina 16S read quality profiles.
    """
    if not quality_data:
        print("  WARNING: Empty quality data — cannot auto-determine trunc position. "
              "Returning 0 (no truncation).", file=sys.stderr)
        return 0

    n_pos = len(quality_data)

    # ── PRIMARY CHECK: smoothed primary percentile ───────────────────────────
    raw_primary = _extract_percentile_values(quality_data, percentile)
    smoothed_primary = _sliding_average(raw_primary, window)

    primary_cut = None
    for i, avg_q in enumerate(smoothed_primary):
        position = i + 1  # 1-indexed
        if position <= primer_len:
            continue
        if avg_q < min_quality:
            # Walk LEFT on the raw values to find the last base still >= threshold.
            # Start one position to the left of the crossing (i-1, 0-indexed: i-1).
            walk_idx = i - 1
            while walk_idx > primer_len and raw_primary[walk_idx] < min_quality:
                walk_idx -= 1
            # walk_idx is now the 0-indexed position of the last acceptable base.
            primary_cut = walk_idx + 1  # convert to 1-indexed
            print(
                f"  PRIMARY: smoothed {percentile} ({window}bp window) drops below "
                f"{min_quality} at position {position}; walking left on raw values "
                f"→ last acceptable base at position {primary_cut} "
                f"(raw val={raw_primary[walk_idx]:.1f})",
                file=sys.stderr,
            )
            break

    if primary_cut is None:
        print(
            f"  PRIMARY: {percentile} quality stays ≥ {min_quality} for all "
            f"{n_pos} positions.",
            file=sys.stderr,
        )

    # ── SECONDARY CHECK: 9th-percentile tail collapse ────────────────────────
    raw_p9 = _extract_percentile_values(quality_data, '9%')
    # scan only the terminal tail_check_window positions (after primer region)
    tail_start_idx = max(primer_len, n_pos - tail_check_window)
    p9_collapse_idx = None
    for i in range(tail_start_idx, n_pos):
        if raw_p9[i] < min_quality_p9:
            p9_collapse_idx = i  # 0-indexed; first collapse position in tail
            break

    secondary_cut = None
    if p9_collapse_idx is not None:
        secondary_cut = p9_collapse_idx + 1  # 1-indexed; exact collapse position, no buffer
        print(
            f"  TAIL CHECK: 9th-pct drops below {min_quality_p9} first at position "
            f"{secondary_cut} (raw val={raw_p9[p9_collapse_idx]:.1f})",
            file=sys.stderr,
        )
    else:
        print(
            f"  TAIL CHECK: 9th-pct ≥ {min_quality_p9} in all tail positions checked.",
            file=sys.stderr,
        )

    # ── Decide final cut position ─────────────────────────────────────────────
    # Take the MOST CONSERVATIVE (earliest) cut among checks that fired.
    candidates = [c for c in [primary_cut, secondary_cut] if c is not None]

    if not candidates:
        print(
            f"  DECISION: Both checks pass — no truncation required. "
            f"(total positions={n_pos})",
            file=sys.stderr,
        )
        return 0

    final_cut = min(candidates)  # most conservative (earliest) among fired checks

    print(
        f"  DECISION: primary_cut={primary_cut}, secondary_cut={secondary_cut} → "
        f"final trunc={final_cut}  (primer_len={primer_len})",
        file=sys.stderr,
    )
    return final_cut


def readlen_cap_position(quality_data: list, retention: float = READLEN_RETENTION) -> int:
    """
    Longest 1-indexed read length that still retains *retention* of the reads,
    derived from the per-position 'count' row of the QIIME2 quality table
    (count[p] = number of sampled reads reaching position p; non-increasing).
    Returns 0 when counts are unavailable.
    """
    counts = [e.get('count') for e in quality_data]
    valid = [c for c in counts if c is not None]
    if not valid:
        return 0
    cmax = max(valid)
    if cmax <= 0:
        return 0
    cap = 0
    for i, c in enumerate(counts):
        if c is not None and c >= retention * cmax:
            cap = i + 1  # 1-indexed length
    return cap


def apply_readlen_cap(trunc: int, quality_data: list, direction: str,
                      retention: float = READLEN_RETENTION) -> tuple:
    """
    Clamp a quality-based truncation length so it never exceeds the read-length
    distribution. Returns (trunc, warning_or_None). trunc == 0 (no truncation)
    is left untouched — DADA2 then applies no length filter, which is already
    safe for variable-length reads.
    """
    if not trunc:
        return trunc, None
    cap = readlen_cap_position(quality_data, retention)
    if cap and trunc > cap:
        warning = (
            f"READ-LENGTH CAP ({direction}): quality-based trunc {trunc} exceeds the "
            f"length that retains {int(retention * 100)}% of reads ({cap}bp); clamped "
            f"to {cap}. (A longer trunc would discard most reads — variable read "
            f"lengths.)"
        )
        print(f"  {warning}", file=sys.stderr)
        return cap, warning
    return trunc, None


def max_quality_length(quality_data: list, percentile: str = '25%',
                       floor: float = 20.0, window: int = 10,
                       primer_len: int = 0) -> int:
    """
    Longest 1-indexed position the read can be truncated to while the chosen
    *percentile* (smoothed, default 25th) stays >= *floor* — i.e. the primary
    quality limit WITHOUT the conservative 9th-percentile tail check.

    This is the ceiling for the overlap-rescue extension: when paired coverage
    is too short to span the amplicon, truncation may be extended up to here
    (overriding the tail cut) because the typical read quality is still good.
    Bounded separately by the read-length cap. Returns len(quality_data) when
    the percentile never drops below the floor.
    """
    if not quality_data:
        return 0
    raw = _extract_percentile_values(quality_data, percentile)
    smoothed = _sliding_average(raw, window)
    for i, avg_q in enumerate(smoothed):
        position = i + 1
        if position <= primer_len:
            continue
        if avg_q < floor:
            walk_idx = i - 1
            while walk_idx > primer_len and raw[walk_idx] < floor:
                walk_idx -= 1
            return walk_idx + 1
    return len(quality_data)


def extend_trunc_for_overlap(trunc_f: int, trunc_r: int,
                             fwd_quality: list, rev_quality: list,
                             trim_f: int, trim_r: int, region: str,
                             percentile: str, floor: float, window: int,
                             default_min_overlap: int = 12) -> dict:
    """
    When paired-end coverage is too short to span the region's amplicon, extend
    the truncations toward the longest still-good-quality length (25th-pct >=
    floor), bounded by the read-length cap, to recover the merge overlap.

    Returns a dict: new trunc_len_f/r, the extension ceilings (max_qual_len_f/r,
    cap_f/r, readlen_f/r, amplicon_min/max), merge_infeasible_hard (bool: True
    when even the longest good-quality coverage cannot span the amplicon →
    single-end is the only option), and an optional 'warning'.
    """
    amp_range = AMPLICON_SIZES.get(region)
    readlen_f = len(fwd_quality)
    readlen_r = len(rev_quality)
    cap_f = readlen_cap_position(fwd_quality)
    cap_r = readlen_cap_position(rev_quality)

    def ceiling(maxq, cap):
        return min(maxq, cap) if cap else maxq

    max_qual_f = ceiling(
        max_quality_length(fwd_quality, percentile, floor, window, trim_f), cap_f)
    max_qual_r = ceiling(
        max_quality_length(rev_quality, percentile, floor, window, trim_r), cap_r)

    info = {
        'readlen_f': readlen_f, 'readlen_r': readlen_r,
        'cap_f': cap_f, 'cap_r': cap_r,
        'max_qual_len_f': max_qual_f, 'max_qual_len_r': max_qual_r,
        'amplicon_min': amp_range[0] if amp_range else None,
        'amplicon_max': amp_range[1] if amp_range else None,
        'trunc_len_f': trunc_f, 'trunc_len_r': trunc_r,
        'merge_infeasible_hard': False,
        'warning': None,
    }
    if not amp_range:
        return info

    amp_min, amp_max = amp_range
    # Only extend (never shorten) toward the good-quality ceiling.
    ext_f = max(trunc_f, max_qual_f)
    ext_r = max(trunc_r, max_qual_r)
    max_coverage = (ext_f - trim_f) + (ext_r - trim_r)

    if max_coverage >= amp_min + default_min_overlap:
        # Extension recovers a feasible merge → adopt the longer truncations.
        if ext_f != trunc_f or ext_r != trunc_r:
            info['warning'] = (
                f"OVERLAP RESCUE: extended trunc to {ext_f}/{ext_r} (from "
                f"{trunc_f}/{trunc_r}) — coverage {max_coverage}bp now spans "
                f"amplicon≈{amp_min}-{amp_max}bp (region={region}). Overrode the "
                f"9th-pct tail cut; 25th-pct stays >= {floor:.0f}."
            )
        info['trunc_len_f'] = ext_f
        info['trunc_len_r'] = ext_r
    else:
        # Even at the longest good-quality length the reads cannot span the
        # amplicon → forward-only is the only way (extreme case).
        info['merge_infeasible_hard'] = True
        info['trunc_len_f'] = ext_f  # best-effort (still used if forced paired)
        info['trunc_len_r'] = ext_r
        info['warning'] = (
            f"MERGE INFEASIBLE (hard): max good-quality coverage {max_coverage}bp "
            f"< amplicon≈{amp_min}-{amp_max}bp + overlap (region={region}). "
            f"Reads cannot span the amplicon — recommend single-end (forward)."
        )
    return info


def compute_stats(
    quality_data: list,
    sample_counts: dict,
    min_quality: float,
    primer_len_f: int,
    primer_len_r: int,
    has_reverse: bool,
    percentile: str = '25%',
    window: int = 10,
    min_quality_p9: float = 15.0,
    tail_check_window: int = 20,
) -> dict:
    """Compute all auto-parameters and warnings."""
    warnings = []

    # ── Forward parameters ───────────────────────────────────────────────
    trunc_f = find_trunc_position(
        quality_data, min_quality, primer_len_f,
        percentile, window, min_quality_p9, tail_check_window
    )
    trunc_f, cap_warning = apply_readlen_cap(trunc_f, quality_data, 'forward')
    if cap_warning:
        warnings.append(cap_warning)
    trim_f = primer_len_f

    # ── Reverse parameters (if paired) ──────────────────────────────────
    trunc_r = 0
    trim_r = 0

    # ── Sequence count statistics ────────────────────────────────────────
    if sample_counts:
        counts_list = list(sample_counts.values())
        total_seqs = sum(counts_list)
        min_seqs = min(counts_list)
        max_seqs = max(counts_list)
        n_samples = len(counts_list)
        median_seqs = sorted(counts_list)[n_samples // 2]
        mean_seqs = total_seqs / n_samples if n_samples > 0 else 0

        # Check for problematic samples
        if min_seqs < 100:
            warnings.append(
                f"CRITICAL: {sum(1 for c in counts_list if c < 100)} samples "
                f"have < 100 sequences (min={min_seqs})"
            )
        if min_seqs == 0:
            warnings.append(
                f"CRITICAL: {sum(1 for c in counts_list if c == 0)} samples "
                f"have 0 sequences"
            )

        # Heterogeneity check
        if mean_seqs > 0:
            import math
            variance = sum((c - mean_seqs) ** 2 for c in counts_list) / n_samples
            cv = math.sqrt(variance) / mean_seqs
            if cv > 2.0:
                warnings.append(
                    f"WARNING: High heterogeneity in sequence counts (CV={cv:.2f})"
                )
    else:
        total_seqs = 0
        min_seqs = 0
        max_seqs = 0
        median_seqs = 0
        n_samples = 0

    # ── Truncation quality checks ────────────────────────────────────────
    if trunc_f > 0:
        usable_len_f = trunc_f - trim_f
        if usable_len_f < 100:
            warnings.append(
                f"WARNING: Usable forward read length after trim+trunc is only "
                f"{usable_len_f}bp (trunc_f={trunc_f}, trim_f={trim_f})"
            )
        if usable_len_f < 50:
            warnings.append(
                f"ERROR: Usable forward read length is critically short "
                f"({usable_len_f}bp). Review quality profiles manually."
            )

    result = {
        'trunc_len_f': trunc_f,
        'trunc_len_r': trunc_r,
        'trim_left_f': trim_f,
        'trim_left_r': trim_r,
        'total_seqs': total_seqs,
        'n_samples': n_samples,
        'min_seqs_per_sample': min_seqs,
        'max_seqs_per_sample': max_seqs,
        'median_seqs_per_sample': median_seqs,
        'auto_params_used': True,
        'warnings': warnings,
    }

    return result


def compute_recommended_min_overlap(
    trunc_f: int, trunc_r: int, trim_f: int, trim_r: int,
    region: str, default_min_overlap: int = 12,
) -> dict:
    """
    Compute the expected overlap and recommend min_overlap for DADA2.

    If the expected overlap is < default_min_overlap (12) but >= 4, the
    recommendation steps down to 6, then 4.  Returns a dict with:
      - expected_overlap_min/max
      - recommended_min_overlap
      - overlap_warning (str or None)
      - merge_feasibility: 'ok' | 'marginal' | 'infeasible' | 'unknown'

    'merge_feasibility' is a machine-readable flag so the denoise QC step can
    correctly attribute a low merge rate to insufficient overlap rather than to
    read quality.  When the region is unknown (not in AMPLICON_SIZES) the overlap
    cannot be predicted: a visible 'unknown' note is emitted instead of silently
    skipping the check.
    """
    amp_range = AMPLICON_SIZES.get(region)
    if not amp_range or not (trunc_f and trunc_r):
        # Unknown region or single-end → cannot predict overlap. Make it visible.
        reason = ('single-end (no reverse read)' if not (trunc_f and trunc_r)
                  else f"region '{region}' not in known amplicon sizes "
                       f"({', '.join(sorted(AMPLICON_SIZES))})")
        note = (
            f"OVERLAP UNKNOWN: expected overlap not checked — {reason}. "
            f"A low merge rate downstream cannot be attributed to overlap vs quality."
        )
        print(f"  {note}", file=sys.stderr)
        return {
            'expected_overlap_min': None,
            'expected_overlap_max': None,
            'recommended_min_overlap': default_min_overlap,
            'overlap_warning': None if not (trunc_f and trunc_r) else note,
            'merge_feasibility': 'unknown',
        }

    amp_min, amp_max = amp_range
    eff_f = trunc_f - trim_f
    eff_r = trunc_r - trim_r
    total_coverage = eff_f + eff_r

    # Worst-case overlap (longest amplicon)
    overlap_min = total_coverage - amp_max
    # Best-case overlap (shortest amplicon)
    overlap_max = total_coverage - amp_min

    warning = None
    feasibility = 'ok'

    if overlap_min >= default_min_overlap:
        # Comfortable overlap → no change
        rec = default_min_overlap
    elif overlap_min >= 6:
        rec = 6
        feasibility = 'marginal'
        warning = (
            f"OVERLAP WARNING: Expected overlap is {overlap_min}-{overlap_max}bp "
            f"(region={region}, coverage={eff_f}+{eff_r}={total_coverage}bp, "
            f"amplicon≈{amp_min}-{amp_max}bp). "
            f"Reduced --p-min-overlap from {default_min_overlap} to {rec}."
        )
    elif overlap_min >= 4:
        rec = 4
        feasibility = 'marginal'
        warning = (
            f"OVERLAP WARNING: Expected overlap is {overlap_min}-{overlap_max}bp "
            f"(region={region}, coverage={eff_f}+{eff_r}={total_coverage}bp, "
            f"amplicon≈{amp_min}-{amp_max}bp). "
            f"Reduced --p-min-overlap from {default_min_overlap} to {rec}. "
            f"Merge quality may be compromised."
        )
    else:
        # Worst-case (longest amplicon) overlap is too small to merge at all, but
        # if the best-case (shortest amplicon) still has >=4bp, step min_overlap
        # down to 4 anyway so that fraction of the amplicon-length distribution
        # isn't rejected purely by an overlap floor it already clears — those
        # amplicons only fail if the actual overlap ends up too erroneous, not
        # because the floor is set higher than their overlap.
        feasibility = 'infeasible'
        if overlap_max >= 4:
            rec = 4
            warning = (
                f"OVERLAP CRITICAL: Expected overlap is {overlap_min}-{overlap_max}bp — "
                f"insufficient for paired-end merging for the longer amplicons in range "
                f"(region={region}, coverage={eff_f}+{eff_r}={total_coverage}bp, "
                f"amplicon≈{amp_min}-{amp_max}bp). "
                f"Reduced --p-min-overlap from {default_min_overlap} to {rec} so the "
                f"shorter amplicons (overlap>={rec}bp) can still merge; longer ones will fail."
            )
        else:
            rec = default_min_overlap  # keep default; reads will simply fail to merge
            warning = (
                f"OVERLAP CRITICAL: Expected overlap is {overlap_min}-{overlap_max}bp — "
                f"insufficient for paired-end merging "
                f"(region={region}, coverage={eff_f}+{eff_r}={total_coverage}bp, "
                f"amplicon≈{amp_min}-{amp_max}bp). "
                f"Most reads will fail to merge."
            )

    if warning:
        print(f"  {warning}", file=sys.stderr)

    return {
        'expected_overlap_min': overlap_min,
        'expected_overlap_max': overlap_max,
        'recommended_min_overlap': rec,
        'overlap_warning': warning,
        'merge_feasibility': feasibility,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Extract quality stats and auto-determine denoise params'
    )
    parser.add_argument('--qzv', required=True, help='Path to demux summary QZV')
    parser.add_argument('--output', required=True, help='Output JSON path')
    parser.add_argument('--min-quality', type=float, default=20,
                        help='Quality threshold for primary (25th pct) truncation (default: 20)')
    parser.add_argument('--min-quality-p9', type=float, default=15,
                        help='Quality threshold for 9th-pct tail check (default: 15)')
    parser.add_argument('--primer-len-f', type=int, default=0,
                        help='Forward primer length for trim-left')
    parser.add_argument('--primer-len-r', type=int, default=0,
                        help='Reverse primer length for trim-left')
    parser.add_argument('--region', type=str, default='',
                        help='16S region (e.g. V3-V4) for overlap estimation')
    parser.add_argument('--percentile', type=str, default='25%',
                        help='Primary percentile for quality check (default: "25%%")')
    parser.add_argument('--window', type=int, default=10,
                        help='Sliding window size to smooth quality spikes (default: 10)')
    parser.add_argument('--tail-check-window', type=int, default=20,
                        help='Last N positions to scan for 9th-pct collapse (default: 20)')
    parser.add_argument('--retry-quality-floor', type=float, default=20,
                        help='25th-pct floor for overlap-rescue trunc extension (default: 20)')
    args = parser.parse_args()

    if not os.path.exists(args.qzv):
        print(f"ERROR: QZV file not found: {args.qzv}", file=sys.stderr)
        sys.exit(1)

    # Extract quality data
    fwd_quality = extract_quality_table(args.qzv, 'forward')
    rev_quality = extract_quality_table(args.qzv, 'reverse')
    sample_counts = extract_sample_counts(args.qzv)

    has_reverse = len(rev_quality) > 0

    # Compute auto-parameters
    result = compute_stats(
        fwd_quality, sample_counts,
        args.min_quality,
        args.primer_len_f, args.primer_len_r,
        has_reverse,
        args.percentile,
        args.window,
        args.min_quality_p9,
        args.tail_check_window,
    )

    # If reverse data exists, compute reverse truncation
    if has_reverse:
        trunc_r = find_trunc_position(
            rev_quality, args.min_quality, args.primer_len_r,
            args.percentile, args.window,
            args.min_quality_p9, args.tail_check_window,
        )
        trunc_r, cap_warning_r = apply_readlen_cap(trunc_r, rev_quality, 'reverse')
        if cap_warning_r:
            result['warnings'].append(cap_warning_r)
        result['trunc_len_r'] = trunc_r
        result['trim_left_r'] = args.primer_len_r

        if trunc_r > 0:
            usable_len_r = trunc_r - args.primer_len_r
            if usable_len_r < 100:
                result['warnings'].append(
                    f"WARNING: Usable reverse read length is {usable_len_r}bp"
                )

    # ── Overlap-aware min_overlap recommendation (paired-end only) ────────
    def store_overlap(info):
        result['expected_overlap_min'] = info['expected_overlap_min']
        result['expected_overlap_max'] = info['expected_overlap_max']
        result['recommended_min_overlap'] = info['recommended_min_overlap']
        result['merge_feasibility'] = info['merge_feasibility']

    if has_reverse and args.region:
        overlap_info = compute_recommended_min_overlap(
            result['trunc_len_f'], result['trunc_len_r'],
            result['trim_left_f'], result['trim_left_r'],
            args.region,
        )
        store_overlap(overlap_info)
        if overlap_info['overlap_warning']:
            result['warnings'].append(overlap_info['overlap_warning'])

        # ── Overlap-rescue truncation extension ──────────────────────────────
        # When coverage is too short to span the amplicon, extend truncation
        # toward the longest still-good-quality length (overriding the 9th-pct
        # tail cut) so the merge becomes feasible — or flag it as physically
        # impossible (→ single-end is the only option).
        result['merge_infeasible_hard'] = False
        if overlap_info['merge_feasibility'] in ('marginal', 'infeasible'):
            ext = extend_trunc_for_overlap(
                result['trunc_len_f'], result['trunc_len_r'],
                fwd_quality, rev_quality,
                result['trim_left_f'], result['trim_left_r'], args.region,
                args.percentile, args.retry_quality_floor, args.window,
            )
            for k in ('readlen_f', 'readlen_r', 'cap_f', 'cap_r',
                      'max_qual_len_f', 'max_qual_len_r',
                      'amplicon_min', 'amplicon_max'):
                result[k] = ext[k]
            result['merge_infeasible_hard'] = ext['merge_infeasible_hard']
            if ext['warning']:
                result['warnings'].append(ext['warning'])
            if (ext['trunc_len_f'], ext['trunc_len_r']) != \
               (result['trunc_len_f'], result['trunc_len_r']):
                result['trunc_len_f'] = ext['trunc_len_f']
                result['trunc_len_r'] = ext['trunc_len_r']
                overlap_info = compute_recommended_min_overlap(
                    result['trunc_len_f'], result['trunc_len_r'],
                    result['trim_left_f'], result['trim_left_r'], args.region,
                )
                store_overlap(overlap_info)
        else:
            result['readlen_f'] = len(fwd_quality)
            result['readlen_r'] = len(rev_quality)
            result['cap_f'] = readlen_cap_position(fwd_quality)
            result['cap_r'] = readlen_cap_position(rev_quality)
            amp = AMPLICON_SIZES.get(args.region, (None, None))
            result['amplicon_min'], result['amplicon_max'] = amp
    else:
        result['recommended_min_overlap'] = 12  # default
        # Single-end, or paired-end without a known region → not predictable.
        result['merge_feasibility'] = 'na' if not has_reverse else 'unknown'
        result['merge_infeasible_hard'] = False
        result['readlen_f'] = len(fwd_quality)
        result['readlen_r'] = len(rev_quality) if has_reverse else 0
        result['cap_f'] = readlen_cap_position(fwd_quality)
        result['cap_r'] = readlen_cap_position(rev_quality) if has_reverse else 0

    # Write output
    with open(args.output, 'w') as f:
        json.dump(result, f, indent=2)

    # Print summary
    print(f"Quality stats: trunc_f={result['trunc_len_f']}, "
          f"trunc_r={result['trunc_len_r']}, "
          f"trim_f={result['trim_left_f']}, "
          f"trim_r={result['trim_left_r']}, "
          f"min_overlap={result['recommended_min_overlap']}, "
          f"n_samples={result['n_samples']}, "
          f"warnings={len(result['warnings'])}",
          file=sys.stderr)

    if result['warnings']:
        for w in result['warnings']:
            print(f"  - {w}", file=sys.stderr)


if __name__ == '__main__':
    main()
