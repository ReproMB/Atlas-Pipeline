# ReproMB Atlas 16S rRNA processing pipeline

Nextflow pipeline used to reprocess the raw 16S rRNA amplicon reads of the
studies included in the [ReproMB Atlas](https://repromb.ugr.es/repromb), an
atlas of the female reproductive tract and gut microbiome. Every dataset goes
from FASTQ files to a genus-level count table with the same software versions,
reference database and decision rules. The parameters that vary between datasets
(truncation lengths, overlap, paired-end versus forward-only denoising) are
derived by fixed rules from properties measured in each dataset.

This repository contains the pipeline, the *Lactobacillus* species-group step,
and the record of how release 2026.11 of the Atlas was produced.

## Processing units

A *processing unit* is one study, one primer pair and one body-site group.
Samples of the same study taken from different reproductive-tract sites with the
same primers are denoised together (unit `PRJEB16013_Cervix_Endometrium_Vagina`,
for example); gut samples always form their own unit. The site of each sample
comes from the curated metadata, so a unit only contains the samples of its own
sites.

## Software

| Component | Version used for release 2026.11 |
|---|---|
| Nextflow | 25.10.4 (build 11173) |
| QIIME 2 | amplicon distribution 2025.7, container `quay.io/qiime2/amplicon:2025.7` (digest `sha256:889c6f646dd6fbe7a6bc350fe0fbfa07922efb9ab3bbc5349d8273cbd8e802e5`) |
| q2-dada2 / DADA2 | 2025.7.0 / DADA2 1.30.0 on R 4.3.3 |
| cutadapt | 5.1 |
| Deblur | 1.1.1 |
| vsearch | 2.22.1 |
| scikit-learn, biom-format | 1.4.2, 2.1.16 |
| Reference | SILVA 138.2 SSURef NR99, prepared with RESCRIPt ([`docs/reference_database.md`](docs/reference_database.md)) |
| *Lactobacillus* training set | The 38 GTDB r95 *Lactobacillus* sensu stricto representatives, labelled with four species groups following Lebeer et al. (2023) (`species/build_database.py`) |

All QIIME 2 steps run in the container. Four small steps (input parsing,
metadata selection, manifest writing and the QC summary) run on the host and
need Python 3 with pandas and openpyxl. The species step needs the QIIME 2
conda environment of the same distribution, which provides R and DADA2.

## Inputs

| Parameter | Content |
|---|---|
| `--input_csv` | Study table, one row per study and body site. Columns used: `StudyProject` (accession), `Site`, `Region`, `Paired/Single`, `ForwardPrimer`, `ReversePrimer`, `ReadLength`. Other columns are carried along and ignored. Examples in [`examples/`](examples/). |
| `--data_dir` | FASTQ files under `<data_dir>/<accession>/`, as downloaded by nf-core/fetchngs. |
| `--study_metadata_dir` | One sample sheet per study (`<accession>.xlsx`, `.csv` or `.tsv`) with a run identifier and the curated columns `ConditionAtlas` and `SiteAtlas` (optionally `AgeAtlas`, `BMIAtlas`, `EthnicityAtlas`, `MenstrualCycle`). Samples without condition or site are not processed. |
| `--silva_seqs`, `--silva_tax` | SILVA 138.2 sequences and taxonomy as QIIME 2 artifacts. |
| `--manual_params` | Only with `--denoise_mode manual`: one row per processing unit; any blank cell falls back to the automatic value. |

The curated per-sample metadata of release 2026.11 can be downloaded from the
Atlas (Download page). Raw reads are public under the accessions listed in
`release/2026.11/units.csv`.

## Running

DADA2 units, truncation chosen automatically:

```bash
nextflow run main.nf -profile docker \
    --input_csv studies.csv \
    --data_dir fastq/ \
    --study_metadata_dir metadata/ \
    --silva_seqs silva-138.2-ssu-nr99-seqs.qza \
    --silva_tax silva-138.2-ssu-nr99-tax.qza \
    --outdir results
```

Units with per-unit values take `--denoise_mode manual --manual_params
<sheet.csv>`; a unit without a row in the sheet keeps the automatic values, and
any blank cell falls back to the automatic value.

Units with synthetic quality scores, denoised with Deblur:

```bash
nextflow run main.nf -profile docker ... \
    --input_csv examples/samplesheet_deblur.csv \
    --denoiser deblur --run_primer_trim true \
    --denoise_mode manual --manual_params examples/manual_params_deblur.csv
```

Species groups for one unit, after the pipeline:

```bash
species/run_species_step.sh results/PRJNA1029035_Vagina PRJNA1029035_Vagina \
    trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz
```

If the FASTQ files are symbolic links into another file system, pass the bind
mount with `--docker_extra_mounts '-v /path:/path'`. Classifiers are the
slowest step; `--classifier_cache_dir results/classifiers` reuses the ones built
in a previous run.

## Outputs

One folder per processing unit, `results/<unit>/`. The files loaded into the
Atlas are in `exports/`: `<unit>_feature_table.tsv` (genus counts, bacteria
only), `<unit>_taxonomy.tsv` (ASV taxonomy) and, after the species step,
`<unit>_resolved_genus_table.tsv`. The other folders hold intermediate QIIME 2
artifacts and per-unit reports.

## Parameters

Defaults are those of release 2026.11 (`nextflow.config`). Run
`nextflow run main.nf --help` for the full list.

## Repository layout

```
main.nf, nextflow.config, conf/   workflow and configuration
modules/, subworkflows/           Nextflow processes
bin/                              scripts called by the processes
species/                          Lactobacillus species-group step
examples/                         input sheets for three release datasets
docs/                             processing methods and reference database
release/2026.11/                  per-dataset parameters of release 2026.11
```

## Use of AI tools

Claude Opus 5.5 (Anthropic) was used to help document the code and improve its
readability. All design decisions and every line of code were reviewed by the
authors, who take full responsibility for the content of this repository.

## Citation and licence

If you use this pipeline, please cite the ReproMB Atlas article (reference to be
added on publication) and this repository ([`CITATION.cff`](CITATION.cff)). The
code is released under the MIT licence ([`LICENSE`](LICENSE)).
