#!/usr/bin/env python3
"""
subclassify_lactobacillus.py

Lactobacillus species subclassification using the GTDB-derived custom
training set from Lebeer et al. 2023 (doi:10.1038/s41564-023-01500-0).

Takes DADA2 outputs from the ATLAS pipeline and produces a per-sample
species-group abundance table distinguishing the 9 Lactobacillaceae subgenera
(crispatus, iners, gasseri, jensenii + 5 Limosilactobacillus groups).

Inputs (auto-detected from --study-dir/--prefix, or override individually):
  rep_seqs.qza         DADA2 representative sequences (QZA → FASTA)
  bacteria_table.qza   Bacteria-only ASV count table  (QZA → BIOM)
    fallback: denoise/XXX_table.qza if bacteria_table not present
  taxonomy.tsv         SILVA ASV-level taxonomy (for Lactobacillaceae filter)

Usage:
  python subclassify_lactobacillus.py \\
      --study-dir PRJNA1431678_Endometrium_Vagina \\
      --prefix    PRJNA1431678_Endometrium_Vagina \\
      --database  database/trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz \\
      [--metadata PRJNA1431678_Endometrium_Vagina/lactobacillus/XXX_lactobacillus_summary.csv] \\
      [--threads 4] [--min-boot 50]
"""

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    import biom
except ImportError:
    sys.exit(
        "ERROR: 'biom' package not found.\n"
        "Activate the QIIME2 environment first:\n"
        "  conda activate qiime2-amplicon-2025.7"
    )

# ── Rscript with DADA2 (release 2026.11: the qiime2-amplicon-2025.7 conda env) ──
# Taken from $RSCRIPT, else the first Rscript on PATH; --rscript-path overrides.
RSCRIPT_PATH = os.environ.get("RSCRIPT") or shutil.which("Rscript") or "Rscript"

R_SCRIPT = Path(__file__).parent / "assign_taxonomy_gtdb.R"

SUBGENERA_ORDER = [
    "Lactobacillus crispatus",
    "Lactobacillus iners",
    "Lactobacillus gasseri",
    "Lactobacillus jensenii",
    "Limosilactobacillus reuteri",
    "Limosilactobacillus coleohominis",
    "Limosilactobacillus fermentum",
    "Limosilactobacillus vaginalis",
    "Limosilactobacillus",
]

SUBGENERA_COLORS = [
    "#2166ac", "#d6604d", "#4dac26", "#762a83",
    "#f4a582", "#92c5de", "#a6dba0", "#e7d4e8",
    "#cccccc",
]


# ── Logging setup ─────────────────────────────────────────────────────────────

def setup_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("lacto_species")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


# ── QZA extraction helpers ────────────────────────────────────────────────────

def extract_from_qza(qza_path: Path, suffix: str, dest_dir: Path, log) -> Path:
    """Extract first entry ending with `suffix` from a QZA (ZIP) to dest_dir."""
    with zipfile.ZipFile(qza_path) as zf:
        matches = [n for n in zf.namelist() if n.endswith(suffix)]
        if not matches:
            raise FileNotFoundError(
                f"No file ending with '{suffix}' found inside {qza_path.name}\n"
                f"Contents: {zf.namelist()[:10]}"
            )
        entry = matches[0]
        dest = dest_dir / Path(entry).name
        log.debug(f"Extracting {entry} → {dest}")
        with zf.open(entry) as src, open(dest, "wb") as dst:
            dst.write(src.read())
    return dest


# ── Input file resolution ─────────────────────────────────────────────────────

def resolve_inputs(args, log):
    study = Path(args.study_dir)
    pfx = args.prefix

    def pick(override, *candidates):
        if override:
            p = Path(override)
            if not p.exists():
                raise FileNotFoundError(f"Override path not found: {p}")
            return p
        for c in candidates:
            p = study / c
            if p.exists():
                return p
        tried = [str(study / c) for c in candidates]
        raise FileNotFoundError(
            f"Could not find input file. Tried:\n  " + "\n  ".join(tried)
        )

    rep_seqs = pick(
        args.rep_seqs,
        f"denoise/{pfx}_rep_seqs.qza",
    )
    table = pick(
        args.table,
        f"taxonomy/{pfx}_bacteria_table.qza",
        f"denoise/{pfx}_table.qza",
    )
    tax_tsv = pick(
        args.taxonomy_tsv,
        f"exports/{pfx}_taxonomy.tsv",
    )

    log.info(f"rep_seqs  : {rep_seqs}")
    log.info(f"table     : {table}")
    log.info(f"taxonomy  : {tax_tsv}")
    return rep_seqs, table, tax_tsv


# ── FASTA filtering ───────────────────────────────────────────────────────────

def filter_fasta_to_lactobacillaceae(all_fasta: Path, lacto_ids: set, out_fasta: Path, log) -> int:
    """Write a new FASTA keeping only sequences whose header is in lacto_ids."""
    written = 0
    with open(all_fasta) as fh, open(out_fasta, "w") as out:
        keep = False
        for line in fh:
            if line.startswith(">"):
                seq_id = line[1:].split()[0].strip()
                keep = seq_id in lacto_ids
                if keep:
                    out.write(line)
                    written += 1
            elif keep:
                out.write(line)
    log.info(f"Filtered FASTA: {written} Lactobacillaceae ASVs → {out_fasta.name}")
    return written


# ── BIOM loading ──────────────────────────────────────────────────────────────

def load_biom_as_dataframe(biom_path: Path, log) -> pd.DataFrame:
    """Load a BIOM file and return a DataFrame (features × samples)."""
    log.debug(f"Loading BIOM: {biom_path}")
    bt = biom.load_table(str(biom_path))
    df = bt.to_dataframe(dense=True)
    log.info(f"BIOM loaded: {df.shape[0]} features × {df.shape[1]} samples")
    return df


# ── R subprocess ──────────────────────────────────────────────────────────────

def run_r_classification(filtered_fasta: Path, db_path: Path,
                         output_tsv: Path, min_boot: int, threads: int,
                         log_path: Path, log, rscript: str) -> Path:
    cmd = [rscript, str(R_SCRIPT),
           str(filtered_fasta), str(db_path),
           str(output_tsv), str(min_boot), str(threads)]
    log.info(f"Running R: {' '.join(cmd)}")

    with open(log_path, "a") as lf:
        result = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"R classification failed (exit {result.returncode}). "
            f"See log: {log_path}"
        )
    log.info(f"R classification complete → {output_tsv.name}")
    return output_tsv


# ── Count aggregation ─────────────────────────────────────────────────────────

def aggregate_by_species(gtdb_df: pd.DataFrame, counts_df: pd.DataFrame, log) -> pd.DataFrame:
    """
    Join GTDB genus assignments onto ASV count rows,
    then sum counts per species-group per sample.
    """
    # gtdb_df has ASV_ID + Genus (+ other levels + bootstraps)
    # counts_df index = Feature ID, columns = sample IDs
    lacto_ids = set(gtdb_df["ASV_ID"])
    counts_sub = counts_df.loc[counts_df.index.intersection(lacto_ids)]

    if counts_sub.empty:
        log.warning("No Lactobacillaceae ASV IDs found in count table — empty species table.")
        return pd.DataFrame()

    merged = counts_sub.join(
        gtdb_df.set_index("ASV_ID")[["Genus"]],
        how="left"
    )
    merged["Genus"] = merged["Genus"].fillna("Unclassified")

    species_counts = merged.groupby("Genus").sum()
    log.info(f"Species groups detected: {list(species_counts.index)}")
    return species_counts


def compute_relabund(species_counts: pd.DataFrame, all_counts_df: pd.DataFrame, log) -> pd.DataFrame:
    """Relative abundance: denominator = total bacteria per sample."""
    total = all_counts_df.sum(axis=0)
    relabund = species_counts.div(total, axis=1).fillna(0)
    log.info("Relative abundance computed (denominator: all bacteria per sample)")
    return relabund


# ── Plots ─────────────────────────────────────────────────────────────────────

def _order_rows(df: pd.DataFrame) -> pd.DataFrame:
    ordered = [g for g in SUBGENERA_ORDER if g in df.index]
    rest = [g for g in df.index if g not in SUBGENERA_ORDER]
    return df.loc[ordered + rest]


def plot_stacked_bar(relabund: pd.DataFrame, out_path: Path, prefix: str, log):
    df = _order_rows(relabund)
    # Sort samples by L. crispatus descending
    if "Lactobacillus crispatus" in df.index:
        sample_order = df.loc["Lactobacillus crispatus"].sort_values(ascending=False).index
    else:
        sample_order = df.columns
    df = df[sample_order]

    colors = {g: SUBGENERA_COLORS[i % len(SUBGENERA_COLORS)]
              for i, g in enumerate(_order_rows(relabund).index)}

    fig, ax = plt.subplots(figsize=(max(12, len(df.columns) * 0.25), 5))
    bottom = np.zeros(len(df.columns))
    for genus in df.index:
        ax.bar(range(len(df.columns)), df.loc[genus], bottom=bottom,
               label=genus, color=colors.get(genus, "#aaaaaa"), width=1.0, linewidth=0)
        bottom += df.loc[genus].values

    ax.set_xlim(-0.5, len(df.columns) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Samples", fontsize=11)
    ax.set_ylabel("Relative Abundance", fontsize=11)
    ax.set_title(f"{prefix}\nLactobacillaceae Species Groups", fontsize=12)
    ax.set_xticks([])
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"Saved stacked bar chart: {out_path.name}")


def plot_boxplot(relabund: pd.DataFrame, metadata: pd.DataFrame, out_path: Path, prefix: str, log):
    df = _order_rows(relabund).T.reset_index()
    df.rename(columns={"index": "sample-id"}, inplace=True)

    id_col = "sample-id" if "sample-id" in metadata.columns else metadata.columns[0]
    meta = metadata.rename(columns={id_col: "sample-id"})

    df_merged = df.merge(meta[["sample-id"] + [c for c in ["ConditionAtlas", "SiteAtlas"] if c in meta.columns]],
                         on="sample-id", how="inner")
    if df_merged.empty:
        log.warning("No metadata–sample overlap for boxplot; skipping.")
        return

    species_cols = [c for c in df.columns if c != "sample-id"]
    df_long = df_merged.melt(id_vars=["sample-id"] + [c for c in ["ConditionAtlas", "SiteAtlas"] if c in df_merged.columns],
                              value_vars=species_cols,
                              var_name="Species_group", value_name="RelAbund")

    n_species = len(species_cols)
    fig, ax = plt.subplots(figsize=(max(10, n_species * 2), 5))
    hue = "SiteAtlas" if "SiteAtlas" in df_long.columns else None
    sns.boxplot(data=df_long, x="Species_group", y="RelAbund", hue=hue,
                ax=ax, fliersize=2)
    ax.set_xlabel("Species group", fontsize=11)
    ax.set_ylabel("Relative Abundance", fontsize=11)
    ax.set_title(f"{prefix} — Lactobacillaceae Species", fontsize=12)
    plt.xticks(rotation=35, ha="right", fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    log.info(f"Saved boxplot: {out_path.name}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Lactobacillus species subclassification (Lebeer et al. 2023 method)"
    )
    parser.add_argument("--study-dir",    required=True,  help="Study result folder")
    parser.add_argument("--prefix",       required=True,  help="Study ID prefix (e.g. PRJNA1431678_Endometrium_Vagina)")
    parser.add_argument("--database",     required=True,  help="Path to GTDB Lactobacillaceae training FASTA (.gz)")
    parser.add_argument("--metadata",     default=None,   help="Optional metadata TSV/CSV (for boxplot)")
    parser.add_argument("--threads",      type=int, default=4)
    parser.add_argument("--min-boot",     type=int, default=50)
    parser.add_argument("--rep-seqs",     default=None,   help="Override rep_seqs.qza path")
    parser.add_argument("--table",        default=None,   help="Override table.qza / bacteria_table.qza path")
    parser.add_argument("--taxonomy-tsv", default=None,   help="Override taxonomy.tsv path")
    parser.add_argument("--output-dir",   default=None,   help="Override output directory")
    parser.add_argument("--rscript-path", default=RSCRIPT_PATH, help="Path to Rscript executable")
    args = parser.parse_args()

    # ── Output directory ──────────────────────────────────────────────────────
    if args.output_dir:
        out_dir = Path(args.output_dir)
    else:
        out_dir = Path(args.study_dir) / "taxonomy" / "species"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / ".tmp"
    tmp_dir.mkdir(exist_ok=True)

    log_path = out_dir / f"{args.prefix}_lacto_classify_log.txt"
    log = setup_logging(log_path)
    log.info(f"Output directory: {out_dir}")
    log.info(f"Prefix: {args.prefix}")

    # ── Validate inputs ───────────────────────────────────────────────────────
    db_path = Path(args.database)
    if not db_path.exists():
        log.error(f"Database not found: {db_path}")
        log.error("Run:  python build_database.py")
        sys.exit(1)

    rscript = shutil.which(args.rscript_path) or args.rscript_path
    if not Path(rscript).exists():
        log.error(f"Rscript not found: {rscript}")
        log.error("Activate the qiime2-amplicon-2025.7 conda environment.")
        sys.exit(1)

    if not R_SCRIPT.exists():
        log.error(f"R script not found: {R_SCRIPT}")
        sys.exit(1)

    rep_seqs_qza, table_qza, tax_tsv_path = resolve_inputs(args, log)

    # ── Step 1: Extract FASTA from rep_seqs.qza ───────────────────────────────
    log.info("Extracting rep_seqs FASTA from QZA...")
    all_fasta = extract_from_qza(rep_seqs_qza, "dna-sequences.fasta", tmp_dir, log)

    # ── Step 2: Filter to Lactobacillaceae ASVs ───────────────────────────────
    log.info("Loading SILVA taxonomy for Lactobacillaceae filter...")
    tax_df = pd.read_csv(tax_tsv_path, sep="\t")
    # Column names can vary: "Feature ID" / "#OTU ID" / first column
    id_col = tax_df.columns[0]
    taxon_col = [c for c in tax_df.columns if "taxon" in c.lower() or "tax" in c.lower()][-1]

    lacto_mask = tax_df[taxon_col].str.contains("f__Lactobacillaceae", na=False)
    lacto_ids = set(tax_df.loc[lacto_mask, id_col].str.strip())
    log.info(f"Lactobacillaceae ASVs in taxonomy: {len(lacto_ids)}")

    if len(lacto_ids) == 0:
        log.warning("No Lactobacillaceae ASVs found in taxonomy.tsv — writing empty outputs.")
        _write_empty_outputs(out_dir, args.prefix, log)
        return

    filtered_fasta = tmp_dir / f"{args.prefix}_lacto_filtered.fasta"
    n_written = filter_fasta_to_lactobacillaceae(all_fasta, lacto_ids, filtered_fasta, log)

    if n_written < 5:
        log.warning(f"Only {n_written} Lactobacillaceae sequences found (<5). "
                    "Results may be sparse.")
    if n_written == 0:
        log.warning("0 Lactobacillaceae sequences in FASTA — writing empty outputs.")
        _write_empty_outputs(out_dir, args.prefix, log)
        return

    # ── Step 3: R classification ──────────────────────────────────────────────
    gtdb_tsv = tmp_dir / f"{args.prefix}_lacto_gtdb_raw.tsv"
    run_r_classification(
        filtered_fasta, db_path, gtdb_tsv,
        args.min_boot, args.threads,
        log_path, log, rscript
    )

    # ── Step 4: Extract BIOM count matrix ────────────────────────────────────
    log.info("Extracting BIOM count table from QZA...")
    biom_path = extract_from_qza(table_qza, "feature-table.biom", tmp_dir, log)
    counts_df = load_biom_as_dataframe(biom_path, log)

    # ── Step 5: Merge classifications with counts ─────────────────────────────
    gtdb_df = pd.read_csv(gtdb_tsv, sep="\t")
    log.info(f"GTDB classifications loaded: {len(gtdb_df)} rows")

    species_counts = aggregate_by_species(gtdb_df, counts_df, log)

    if species_counts.empty:
        log.warning("Empty species count table. Writing placeholder outputs.")
        _write_empty_outputs(out_dir, args.prefix, log)
        return

    # ── Step 6: Relative abundance ────────────────────────────────────────────
    species_relabund = compute_relabund(species_counts, counts_df, log)

    # ── Step 7: Write ASV-level classification table ──────────────────────────
    lacto_counts_sub = counts_df.loc[counts_df.index.intersection(set(gtdb_df["ASV_ID"]))]
    silva_genus = tax_df.set_index(id_col)[taxon_col].str.extract(r"g__([^;]+)$")[0].str.strip()

    asv_out = gtdb_df[["ASV_ID", "Genus"]].copy()
    asv_out.rename(columns={"Genus": "GTDB_species_group"}, inplace=True)
    asv_out["SILVA_genus"] = asv_out["ASV_ID"].map(silva_genus)
    boot_cols = [c for c in gtdb_df.columns if c.startswith("boot_")]
    for bc in boot_cols:
        asv_out[bc] = gtdb_df[bc].values
    asv_out = asv_out.join(
        lacto_counts_sub.T.describe().T[["mean", "max"]].rename(
            columns={"mean": "mean_count", "max": "max_count"}
        ),
        on="ASV_ID"
    )

    asv_path = out_dir / f"{args.prefix}_lacto_asv_classification.tsv"
    asv_out.to_csv(asv_path, sep="\t", index=False)
    log.info(f"Saved ASV classification: {asv_path.name}")

    # ── Step 8: Write count and relabund tables ───────────────────────────────
    counts_path = out_dir / f"{args.prefix}_lacto_species_counts.tsv"
    relabund_path = out_dir / f"{args.prefix}_lacto_species_relabund.tsv"
    species_counts.to_csv(counts_path, sep="\t")
    species_relabund.to_csv(relabund_path, sep="\t")
    log.info(f"Saved species counts   : {counts_path.name}")
    log.info(f"Saved species relabund : {relabund_path.name}")

    # ── Step 9: Plots ─────────────────────────────────────────────────────────
    stackedbar_path = out_dir / f"{args.prefix}_lacto_species_stackedbar.png"
    plot_stacked_bar(species_relabund, stackedbar_path, args.prefix, log)

    if args.metadata:
        meta_path = Path(args.metadata)
        if meta_path.exists():
            sep = "\t" if meta_path.suffix == ".tsv" else ","
            meta_df = pd.read_csv(meta_path, sep=sep)
            boxplot_path = out_dir / f"{args.prefix}_lacto_species_boxplot.png"
            plot_boxplot(species_relabund, meta_df, boxplot_path, args.prefix, log)
        else:
            log.warning(f"Metadata file not found: {meta_path}")

    # ── Summary ───────────────────────────────────────────────────────────────
    log.info("=" * 60)
    log.info("DONE")
    log.info(f"  Species groups: {len(species_counts)}")
    log.info(f"  Samples       : {species_counts.shape[1]}")
    log.info(f"  Outputs in    : {out_dir}")
    log.info("=" * 60)


def _write_empty_outputs(out_dir: Path, prefix: str, log):
    for suffix in ["_lacto_species_counts.tsv", "_lacto_species_relabund.tsv",
                   "_lacto_asv_classification.tsv"]:
        p = out_dir / f"{prefix}{suffix}"
        pd.DataFrame().to_csv(p, sep="\t")
    log.info(f"Empty placeholder outputs written to {out_dir}")


if __name__ == "__main__":
    main()
