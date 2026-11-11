# Examples

Input rows of three datasets of release 2026.11.

| Unit | Processing | Files |
|---|---|---|
| PRJNA1029035_Vagina | DADA2, paired-end, V4, automatic mode; 25 samples, 14 MB of FASTQ | `samplesheet_dada2.csv` |
| PRJNA700829_Gut | DADA2, single-end, V3–V4, manual mode with an empty row (every value falls back to the automatic one); 60 samples, 301 MB | `samplesheet_dada2.csv`, `manual_params_dada2.csv` |
| PRJNA734262_Gut | Deblur, forward-only, V3–V4, trim length 250 bp; 87 samples, 152 MB | `samplesheet_deblur.csv`, `manual_params_deblur.csv` |

The two DADA2 units can be run together in manual mode: PRJNA1029035_Vagina
has no row in `manual_params_dada2.csv`, so it keeps the automatic values, as
in its original run. The free-text `Notes` column, which the pipeline does not
read, has been translated to English.

Commands (paths shortened):

```bash
# PRJNA1029035_Vagina
nextflow run main.nf -profile docker \
    --input_csv samplesheet_dada2.csv --data_dir fastq/ --study_metadata_dir metadata/ \
    --silva_seqs silva-138.2-ssu-nr99-seqs.qza --silva_tax silva-138.2-ssu-nr99-tax.qza \
    --outdir results

# PRJNA700829_Gut
nextflow run main.nf -profile docker \
    --input_csv samplesheet_dada2.csv --data_dir fastq/ --study_metadata_dir metadata/ \
    --silva_seqs silva-138.2-ssu-nr99-seqs.qza --silva_tax silva-138.2-ssu-nr99-tax.qza \
    --denoise_mode manual --manual_params manual_params_dada2.csv \
    --classifier_cache_dir previous_run/classifiers/ \
    --outdir results

# PRJNA734262_Gut
nextflow run main.nf -profile docker \
    --input_csv samplesheet_deblur.csv --data_dir fastq/ --study_metadata_dir metadata/ \
    --silva_seqs silva-138.2-ssu-nr99-seqs.qza --silva_tax silva-138.2-ssu-nr99-tax.qza \
    --denoiser deblur --run_primer_trim true \
    --denoise_mode manual --manual_params manual_params_deblur.csv \
    --outdir results
```

Running them needs the FASTQ files (public under the three accessions), the
curated sample sheets of the studies (`metadata/<accession>.xlsx`), SILVA 138.2
and the QIIME 2 container.
