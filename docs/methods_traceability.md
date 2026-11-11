# Processing methods

Summary of the processing steps described in the Methods of the ReproMB Atlas
article and of where each one is implemented.

## 1. Processing units

Samples of one study taken from different reproductive-tract sites with the
same primers are processed together; gut samples form their own unit. Only
samples with a curated condition and site enter the pipeline (`main.nf`,
`bin/parse_study_metadata.py`).

## 2. Primer removal

Primers are removed by position (DADA2 trim-left equal to the primer length)
or, when they occur at variable positions, by sequence with cutadapt
(`modules/cutadapt`, `--run_primer_trim`).

## 3. Truncation

Truncation lengths are chosen per dataset from the per-position quality
profile: reads are cut where the 25th-percentile quality drops below Q20, or
the 9th-percentile quality below Q15, without exceeding the length retained by
80% of the reads. For paired-end data the truncation is extended when needed to
keep enough overlap between forward and reverse reads for the expected amplicon
size (`bin/parse_quality_stats.py`).

## 4. Denoising

Amplicon sequence variants are inferred with DADA2 (maxEE 2, truncQ 2,
consensus chimera removal). Paired-end datasets whose truncated reads cannot
span the amplicon are denoised from the forward reads only, and Ion Torrent and
454 reads are processed as single-end (`modules/denoise`).

Datasets with synthetic quality scores (a constant Phred score along the read)
are denoised with Deblur after primer removal, read merging and trimming to a
common length (`--denoiser deblur`, `modules/deblur`).

## 5. Taxonomic assignment

One naive Bayes classifier is trained per primer pair on SILVA 138.2 SSU NR99,
restricted to the amplified fragment, and reused for every dataset with the
same primers (`modules/classifier_build`, [`reference_database.md`](reference_database.md)).
ASVs are classified with `classify-sklearn` (confidence 0.7); only features
assigned to Bacteria are kept before collapsing to genus (`modules/taxonomy`).

## 6. *Lactobacillus* species groups

In reproductive-tract datasets, ASVs assigned to Lactobacillaceae are classified
against a GTDB r95 training set labelled with four species groups (*L.
crispatus*, *L. iners*, *L. gasseri* and *L. jensenii* groups), following
Lebeer et al. (2023), with dada2 `assignTaxonomy` (minimum bootstrap 50,
reverse-complement matching) and a consensus of 10 runs ([`../species/`](../species/)).

## 7. Steps outside this pipeline

These steps are described in the Methods but are carried out when the Atlas
tables are built, in the web platform's data-loading code:

- removal of mitochondrial and chloroplast reads (SILVA `f__Mitochondria`, `o__Chloroplast`) before relative abundances are computed;
- the within-site filter of rare genera (mean relative abundance below 0.01% and below 1% in every sample), always keeping *Lactobacillus*;
- removal of samples without reads assigned to any genus;
- the contaminant watchlist and the per-study contamination flags.
