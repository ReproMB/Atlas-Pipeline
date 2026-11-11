#!/usr/bin/env bash
# Lactobacillus species-group step for one processing unit, as run for release
# 2026.11: N independent runs of subclassify_lactobacillus.py, their consensus,
# and the genus table with Lactobacillus resolved to species groups.
#
# Usage:
#   run_species_step.sh <unit_dir> <prefix> <training_set.fasta.gz> [runs=10] [threads=4]
#
#   <unit_dir>  pipeline output folder of the unit (results/<prefix>), which holds
#               denoise/<prefix>_rep_seqs.qza, taxonomy/<prefix>_bacteria_table.qza
#               and exports/<prefix>_taxonomy.tsv
#
# Needs python with biom-format, pandas, matplotlib and seaborn, and Rscript with
# dada2 (release 2026.11: the qiime2-amplicon-2025.7 conda environment; set
# RSCRIPT if Rscript is not on PATH).
#
# Output:
#   <unit_dir>/taxonomy/species/runs/run_XX/            one folder per run
#   <unit_dir>/taxonomy/species/<prefix>_lacto_asv_classification_consensus.tsv
#   <unit_dir>/exports/<prefix>_resolved_genus_table.tsv

set -euo pipefail

if [[ $# -lt 3 ]]; then
    sed -n '6,12p' "$0"
    exit 1
fi

UNIT_DIR=$1
PREFIX=$2
DATABASE=$3
RUNS=${4:-10}
THREADS=${5:-4}
HERE=$(cd "$(dirname "$0")" && pwd)
SPECIES_DIR="${UNIT_DIR}/taxonomy/species"

for k in $(seq -f '%02g' 1 "${RUNS}"); do
    RUN_DIR="${SPECIES_DIR}/runs/run_${k}"
    python "${HERE}/subclassify_lactobacillus.py" \
        --study-dir "${UNIT_DIR}" \
        --prefix "${PREFIX}" \
        --database "${DATABASE}" \
        --threads "${THREADS}" \
        --output-dir "${RUN_DIR}"
    rm -rf "${RUN_DIR}/.tmp"
done

CONSENSUS="${SPECIES_DIR}/${PREFIX}_lacto_asv_classification_consensus.tsv"
python "${HERE}/consensus.py" "${SPECIES_DIR}"/runs/run_*/"${PREFIX}"_lacto_asv_classification.tsv \
    --output "${CONSENSUS}"

python "${HERE}/export_resolved_table.py" \
    --study-dir "${UNIT_DIR}" \
    --prefix "${PREFIX}" \
    --lacto-classification "${CONSENSUS}"
