#!/usr/bin/env python3
"""
build_database.py

Builds the custom Lactobacillus-subgenera DADA2 training database following
the method of Lebeer et al. 2023 (doi:10.1038/s41564-023-01500-0).

The paper's original pipeline (scripts-side/01-create-amplicon-dbs/) uses three
Python helpers (reformat_refdb.py, subset_refdb.py, adapt_refdb.py) that all
rely on Bio.SeqIO. Those helpers are re-implemented here with identical logic
using pure-Python FASTA parsing (no BioPython required).

Pipeline (mirrors 03_create_refdbs.sh exactly):
  1. Download GTDB r95 SSU representative FASTA tar.gz
  2. reformat_refdb()  → strip GTDB prefixes (d__, p__…), make DADA2 headers
  3. subset_refdb()    → keep only family Lactobacillaceae  (level 4)
  4. subset_refdb()    → keep only genus  Lactobacillus     (level 5)
  5. adapt_refdb()     → relabel species with the four species groups (changes.csv)

Output: SSUrRNA_GTDB05-lactobacillus-subgenera-rep_DADA2.fna[.gz]

Usage:
  python build_database.py \\
      --output database/trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz \\
      [--ssu-tar  /path/to/bac120_ssu_reps_r95.tar.gz]
      [--changes-csv  /path/to/changes.csv]
"""

import argparse
import gzip
import logging
import re
import sys
import tarfile
import urllib.request
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("build_db")

GTDB_SSU_REPS_URL = (
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/"
    "genomic_files_reps/bac120_ssu_reps_r95.tar.gz"
)

# ── changes.csv content ────────────────────────────────────────────────────────
# Species-to-group mapping, following Lebeer et al. (2023).
# Format: (target_lvl, target_taxon, update_lvl, update_taxon)
#   target_lvl 6 = species, update_lvl 5 = genus  (0-indexed, after prefix strip)
# Species not listed here remain under genus "Lactobacillus".
CHANGES = [
    # ── Lactobacillus crispatus group (acidophilus complex) ──────────────────
    (6, "Lactobacillus crispatus",       5, "Lactobacillus crispatus"),
    (6, "Lactobacillus acidophilus",     5, "Lactobacillus crispatus"),
    (6, "Lactobacillus amylolyticus",    5, "Lactobacillus crispatus"),
    (6, "Lactobacillus amylovorus",      5, "Lactobacillus crispatus"),
    (6, "Lactobacillus gallinarum",      5, "Lactobacillus crispatus"),
    (6, "Lactobacillus helveticus",      5, "Lactobacillus crispatus"),
    (6, "Lactobacillus kefiranofaciens", 5, "Lactobacillus crispatus"),
    (6, "Lactobacillus kitasatonis",     5, "Lactobacillus crispatus"),
    # ── Lactobacillus iners (monotypic) ──────────────────────────────────────
    (6, "Lactobacillus iners",           5, "Lactobacillus iners"),
    # ── Lactobacillus gasseri group ───────────────────────────────────────────
    (6, "Lactobacillus gasseri",         5, "Lactobacillus gasseri"),
    (6, "Lactobacillus johnsonii",       5, "Lactobacillus gasseri"),
    (6, "Lactobacillus paragasseri",     5, "Lactobacillus gasseri"),
    # ── Lactobacillus jensenii group ──────────────────────────────────────────
    (6, "Lactobacillus jensenii",        5, "Lactobacillus jensenii"),
    (6, "Lactobacillus jensenii_A",      5, "Lactobacillus jensenii"),  # GTDB r95 alt name
    (6, "Lactobacillus mulieris",        5, "Lactobacillus jensenii"),
    (6, "Lactobacillus psittaci",        5, "Lactobacillus jensenii"),
]


# ── Pure-Python FASTA parser (replaces Bio.SeqIO) ────────────────────────────

class FastaRecord:
    __slots__ = ("description", "seq")

    def __init__(self, description: str, seq: str):
        self.description = description  # full header text (without '>')
        self.seq = seq


def parse_fasta(handle):
    """Yield FastaRecord objects from an open text file handle."""
    description = None
    seq_parts = []
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if description is not None:
                yield FastaRecord(description, "".join(seq_parts))
            description = line[1:]  # strip the '>'
            seq_parts = []
        else:
            seq_parts.append(line)
    if description is not None:
        yield FastaRecord(description, "".join(seq_parts))


def write_fasta_record(record: FastaRecord, handle):
    handle.write(f">{record.description}\n{record.seq}\n")


# ── The three helper functions — same logic as the paper's scripts ────────────

def reformat_refdb(fin: str, fout: str, remove_prefices: bool = True):
    """
    Reformat GTDB SSU FASTA → DADA2 format.
    Same logic as scripts-side/01-create-amplicon-dbs/reformat_refdb.py,
    with Bio.SeqIO replaced by parse_fasta/write_fasta_record.

    GTDB r95 SSU header format:
      >accession d__Bacteria;...;g__Genus;s__Genus species [genome name]
    Tokens [1] and [2] (space-split) are joined to reconstruct the full
    lineage including the two-word species name.
    """
    with open(fout, "w"):
        pass

    with open(fin, "r") as hin, open(fout, "a") as hout:
        for record in parse_fasta(hin):
            tokens = record.description.split(" ")
            classification = tokens[1] + " " + tokens[2]
            classification = classification.split(";")
            if remove_prefices:
                classification = [re.sub("[a-z]__", "", taxon) for taxon in classification]
            new_id = ";".join(classification)
            record.description = new_id
            write_fasta_record(record, hout)


def subset_refdb(fin: str, fout: str, tax_lvl: int, taxon: str):
    """
    Subset DADA2 FASTA to records matching taxon at tax_lvl.
    Same logic as scripts-side/01-create-amplicon-dbs/subset_refdb.py.
    tax_lvl: 0=domain 1=phylum 2=class 3=order 4=family 5=genus 6=species
    """
    with open(fout, "w"):
        pass

    with open(fin, "r") as hin, open(fout, "a") as hout:
        for record in parse_fasta(hin):
            classification = record.description.split(";")
            if classification[tax_lvl] == taxon:
                write_fasta_record(record, hout)


def adapt_refdb(fin: str, fout: str, changes: list):
    """
    Apply taxonomic relabeling from a changes list.
    Same logic as scripts-side/01-create-amplicon-dbs/adapt_refdb.py,
    with the changes.csv pre-loaded as a list of tuples.
    """
    with open(fout, "w"):
        pass

    with open(fin, "r") as hin, open(fout, "a") as hout:
        for record in parse_fasta(hin):
            classification = record.description.split(";")
            for target_lvl, target_taxon, update_lvl, update_taxon in changes:
                if classification[target_lvl] == target_taxon:
                    classification[update_lvl] = update_taxon
            record.description = ";".join(classification)
            write_fasta_record(record, hout)


# ─────────────────────────────────────────────────────────────────────────────

def count_fasta(path: str) -> int:
    n = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def download_ssu(dest: Path):
    log.info("Downloading GTDB r95 SSU representative sequences...")
    log.info(f"  {GTDB_SSU_REPS_URL}")

    def _progress(count, block_size, total_size):
        if total_size > 0:
            pct = min(count * block_size * 100 / total_size, 100)
            print(f"\r  {pct:5.1f}%", end="", flush=True)

    urllib.request.urlretrieve(GTDB_SSU_REPS_URL, dest, reporthook=_progress)
    print()
    log.info(f"  {dest.stat().st_size / 1e6:.0f} MB downloaded")


def extract_fna_from_tar(tar_path: Path, work_dir: Path) -> Path:
    log.info(f"Extracting FNA from {tar_path.name}...")
    with tarfile.open(tar_path, "r:gz") as tf:
        fna_members = [m for m in tf.getmembers()
                       if m.name.endswith(".fna") or m.name.endswith(".fa")]
        if not fna_members:
            raise FileNotFoundError(f"No .fna file inside {tar_path}")
        member = fna_members[0]
        tf.extract(member, path=work_dir)
        # Handle potential nested path in tar
        candidate = work_dir / Path(member.name).name
        if not candidate.exists():
            candidate = work_dir / member.name
    log.info(f"  {candidate.name}  ({candidate.stat().st_size / 1e6:.0f} MB)")
    return candidate


def load_changes_csv(path: Path) -> list:
    """Load changes.csv produced by 02_create_lactobacillus_subgenera.R.
    Format (no header): target_lvl,target_taxon,update_lvl,update_taxon"""
    changes = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",", 3)
            if len(parts) != 4:
                log.warning(f"Skipping malformed line: {line!r}")
                continue
            changes.append((int(parts[0]), parts[1], int(parts[2]), parts[3]))
    log.info(f"Loaded {len(changes)} entries from {path}")
    return changes


def write_changes_csv(changes: list, path: Path):
    with open(path, "w") as fh:
        for t_lvl, t_taxon, u_lvl, u_taxon in changes:
            fh.write(f"{t_lvl},{t_taxon},{u_lvl},{u_taxon}\n")


def compress_to_gz(src: Path, dest: Path):
    log.info(f"Compressing → {dest.name}")
    with open(src, "rb") as f_in, gzip.open(dest, "wb") as f_out:
        f_out.write(f_in.read())


def main():
    parser = argparse.ArgumentParser(
        description="Build Lactobacillus-subgenera DADA2 training set (Lebeer et al. 2023)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output path (e.g. database/trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz)"
    )
    parser.add_argument(
        "--ssu-tar", default=None,
        help="Local bac120_ssu_reps_r95.tar.gz — skips download"
    )
    parser.add_argument(
        "--changes-csv", default=None,
        help="Path to changes.csv from 02_create_lactobacillus_subgenera.R "
             "(uses embedded mapping if omitted)"
    )
    parser.add_argument(
        "--work-dir", default=None,
        help="Working directory for intermediate files"
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep intermediate FASTA files"
    )
    args = parser.parse_args()

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    work_dir = Path(args.work_dir) if args.work_dir else output_path.parent / ".build_tmp"
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. Obtain SSU tar.gz
    if args.ssu_tar:
        ssu_tar = Path(args.ssu_tar)
        if not ssu_tar.exists():
            log.error(f"SSU tar not found: {ssu_tar}")
            sys.exit(1)
    else:
        ssu_tar = work_dir / "bac120_ssu_reps_r95.tar.gz"
        if not ssu_tar.exists():
            download_ssu(ssu_tar)
        else:
            log.info(f"Using cached: {ssu_tar}")

    # 2. Extract .fna
    fna_raw = extract_fna_from_tar(ssu_tar, work_dir)

    # 3. reformat_refdb: GTDB → DADA2 headers (strip d__/p__/… prefixes)
    fna_bac = str(work_dir / "SSUrRNA_GTDB05-bac-rep_DADA2.fna")
    log.info("reformat_refdb  (strip prefixes)...")
    reformat_refdb(str(fna_raw), fna_bac)
    log.info(f"  {count_fasta(fna_bac):,} sequences")

    if not args.keep_intermediates:
        fna_raw.unlink(missing_ok=True)

    # 4. subset_refdb: family Lactobacillaceae (level 4)
    fna_lac = str(work_dir / "SSUrRNA_GTDB05-lactobacillaceae-rep_DADA2.fna")
    log.info("subset_refdb    (family=Lactobacillaceae, level 4)...")
    subset_refdb(fna_bac, fna_lac, tax_lvl=4, taxon="Lactobacillaceae")
    log.info(f"  {count_fasta(fna_lac):,} sequences")

    if not args.keep_intermediates:
        Path(fna_bac).unlink(missing_ok=True)

    # 5. subset_refdb: genus Lactobacillus (level 5)
    fna_tmp = str(work_dir / "SSUrRNA_GTDB05-lactobacillus-rep_DADA2.fna.tmp")
    log.info("subset_refdb    (genus=Lactobacillus, level 5)...")
    subset_refdb(fna_lac, fna_tmp, tax_lvl=5, taxon="Lactobacillus")
    log.info(f"  {count_fasta(fna_tmp):,} sequences")

    # 6. Load changes
    if args.changes_csv:
        changes = load_changes_csv(Path(args.changes_csv))
    else:
        changes = CHANGES
        changes_out = output_path.parent / "changes.csv"
        write_changes_csv(changes, changes_out)
        log.info(f"Embedded subgenus mapping ({len(changes)} entries) → {changes_out}")
        log.info(
            "  For the exact paper mapping, run 02_create_lactobacillus_subgenera.R\n"
            "  then re-run with: --changes-csv results/lactobacillus_subgenera/changes.csv"
        )

    # 7. adapt_refdb: apply subgenus labels
    fna_subgenera = str(work_dir / "SSUrRNA_GTDB05-lactobacillus-subgenera-rep_DADA2.fna")
    log.info("adapt_refdb     (apply subgenus labels)...")
    adapt_refdb(fna_tmp, fna_subgenera, changes)
    log.info(f"  {count_fasta(fna_subgenera):,} sequences")

    Path(fna_tmp).unlink(missing_ok=True)

    # 8. Write final output
    if str(output_path).endswith(".gz"):
        compress_to_gz(Path(fna_subgenera), output_path)
        if not args.keep_intermediates:
            Path(fna_subgenera).unlink(missing_ok=True)
    else:
        Path(fna_subgenera).rename(output_path)

    # Keep the Lactobacillaceae (genus-level) DB alongside
    lac_out = output_path.parent / "SSUrRNA_GTDB05-lactobacillaceae-rep_DADA2.fna"
    Path(fna_lac).rename(lac_out)

    if not args.keep_intermediates and work_dir.exists():
        remaining = list(work_dir.iterdir())
        if not remaining:
            work_dir.rmdir()

    log.info("=" * 60)
    log.info(f"Subgenera DB        : {output_path}  ({output_path.stat().st_size / 1e6:.1f} MB)")
    log.info(f"Lactobacillaceae DB : {lac_out}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
