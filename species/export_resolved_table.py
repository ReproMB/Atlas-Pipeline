#!/usr/bin/env python3
"""
export_resolved_table.py

Exports a genus-level feature table identical in format to the pipeline's
feature_table.tsv, but with Lactobacillaceae ASVs resolved to GTDB subgroups
(L. crispatus / L. iners / L. gasseri / L. jensenii) instead of the flat
SILVA g__Lactobacillus label.

All other taxa are unchanged (SILVA genus label as-is).

Requires that subclassify_lactobacillus.py has already been run for the study.

Usage:
    python export_resolved_table.py \\
        --study-dir PRJNA397378_Cervix \\
        --prefix    PRJNA397378_Cervix

Output: {study_dir}/exports/{prefix}_resolved_genus_table.tsv
"""

import argparse
import logging
import zipfile
from pathlib import Path

import biom
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("export_resolved")

SUBGROUPS = {
    "Lactobacillus crispatus",
    "Lactobacillus iners",
    "Lactobacillus gasseri",
    "Lactobacillus jensenii",
}


def extract_biom_from_qza(qza_path: Path) -> biom.Table:
    with zipfile.ZipFile(qza_path) as zf:
        biom_entry = next(e for e in zf.namelist() if e.endswith("feature-table.biom"))
        with zf.open(biom_entry) as fh:
            data = fh.read()
    tmp = qza_path.parent / ".tmp_biom"
    tmp.write_bytes(data)
    table = biom.load_table(str(tmp))
    tmp.unlink()
    return table


def load_silva_taxonomy(tsv_path: Path) -> dict:
    df = pd.read_csv(tsv_path, sep="\t", comment="#")
    df.columns = ["Feature ID", "Taxon", "Confidence"][: len(df.columns)]
    return dict(zip(df["Feature ID"], df["Taxon"]))


def normalize_silva_label(taxon: str) -> str:
    return ";".join(p.strip() for p in taxon.split(";"))


def build_resolved_label(asv_id: str, silva_taxon: dict, gtdb_group: dict) -> str:
    raw = silva_taxon.get(asv_id, "")
    label = normalize_silva_label(raw) if raw else "Unclassified"

    group = gtdb_group.get(asv_id, "")
    if group in SUBGROUPS:
        # Replace trailing g__Lactobacillus with g__{subgroup}
        parts = label.split(";")
        for i, part in enumerate(parts):
            if part.strip() == "g__Lactobacillus":
                parts[i] = f"g__{group}"
                break
        label = ";".join(parts)

    return label


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study-dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--table", help="Override bacteria_table.qza path")
    parser.add_argument("--taxonomy-tsv", help="Override taxonomy.tsv path")
    parser.add_argument("--lacto-classification", help="Override lacto_asv_classification.tsv path")
    parser.add_argument("--output-dir", help="Output directory (default: study_dir/exports)")
    args = parser.parse_args()

    study = Path(args.study_dir)
    prefix = args.prefix

    table_path = Path(args.table) if args.table else study / "taxonomy" / f"{prefix}_bacteria_table.qza"
    if not table_path.exists():
        table_path = study / "denoise" / f"{prefix}_table.qza"
    taxonomy_path = Path(args.taxonomy_tsv) if args.taxonomy_tsv else study / "exports" / f"{prefix}_taxonomy.tsv"
    lacto_path = Path(args.lacto_classification) if args.lacto_classification else \
        study / "taxonomy" / "species" / f"{prefix}_lacto_asv_classification.tsv"
    out_dir = Path(args.output_dir) if args.output_dir else study / "exports"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{prefix}_resolved_genus_table.tsv"

    log.info("table      : %s", table_path)
    log.info("taxonomy   : %s", taxonomy_path)
    log.info("lacto cls  : %s", lacto_path)

    if not lacto_path.exists():
        log.error("Lactobacillus classification not found. Run subclassify_lactobacillus.py first.")
        raise SystemExit(1)

    log.info("Loading BIOM count table...")
    table = extract_biom_from_qza(table_path)
    counts = pd.DataFrame(
        np.array(table.matrix_data.todense()),
        index=table.ids(axis="observation"),
        columns=table.ids(axis="sample"),
    )
    log.info("  %d ASVs × %d samples", *counts.shape)

    log.info("Loading SILVA taxonomy...")
    silva_taxon = load_silva_taxonomy(taxonomy_path)

    log.info("Loading GTDB Lactobacillus classification...")
    lacto_df = pd.read_csv(lacto_path, sep="\t")
    gtdb_group = dict(zip(lacto_df["ASV_ID"], lacto_df["GTDB_species_group"].fillna("")))

    log.info("Building resolved genus labels...")
    labels = {asv: build_resolved_label(asv, silva_taxon, gtdb_group) for asv in counts.index}

    counts.index = pd.Index([labels[a] for a in counts.index])
    resolved = counts.groupby(counts.index).sum()
    log.info("  %d genus/subgroup rows after collapse", len(resolved))

    subgroup_rows = [i for i in resolved.index if any(sg in i for sg in SUBGROUPS)]
    log.info("  Lactobacillus subgroup rows: %s", subgroup_rows)

    log.info("Writing: %s", out_path)
    with open(out_path, "w") as fh:
        fh.write("# Constructed from biom file\n")
        fh.write("#OTU ID\t" + "\t".join(resolved.columns) + "\n")
        for label, row in resolved.iterrows():
            vals = "\t".join(str(v) for v in row.values)
            fh.write(f"{label}\t{vals}\n")

    log.info("Done. %d rows × %d samples → %s", *resolved.shape, out_path.name)


if __name__ == "__main__":
    main()
