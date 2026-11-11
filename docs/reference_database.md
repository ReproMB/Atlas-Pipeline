# Reference database

The pipeline takes SILVA as two QIIME 2 artifacts (`--silva_seqs`,
`--silva_tax`) and does not download them itself. For release 2026.11 they were
built with RESCRIPt from SILVA 138.2 SSURef NR99, with rank propagation and
without species labels, and the RNA sequences were reverse-transcribed to DNA.
No culling, dereplication or length filtering was applied:

```bash
qiime rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --p-rank-propagation \
    --p-no-include-species-labels \
    --o-silva-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.2-ssu-nr99-tax.qza

qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.2-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.2-ssu-nr99-seqs.qza
```

Because ranks are propagated, an ASV without a genus-level assignment can carry
the name of a higher rank in the genus position after `taxa collapse --p-level 6`.
The Atlas keeps these labels as they are.

## Classifiers

From these artifacts the pipeline builds one naive Bayes classifier per primer
pair (`subworkflows/classifier_preparation.nf`):

1. `feature-classifier extract-reads` with the primer pair, 80% identity, both
   read orientations and region-specific length limits
   (`modules/classifier_build/main.nf`);
2. `feature-classifier fit-classifier-naive-bayes` with alpha 0.001.

Each classifier is written as `classifier_<region>_<hash>.qza`, where the hash
is derived from the primer sequences (`bin/clean_primers.py`), and is reused for
every unit with the same primer pair.

## *Lactobacillus* training set

The species-group step uses `trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz`,
built with `species/build_database.py` following the method of Lebeer et al.
(2023), *Nat Microbiol* 8:2183–2195 (doi:10.1038/s41564-023-01500-0): the GTDB r95 16S representatives are reduced to family Lactobacillaceae
(325 sequences) and then to genus *Lactobacillus* sensu stricto (38 sequences),
whose species are relabelled with four groups (*L. crispatus*, *L. iners*,
*L. gasseri* and *L. jensenii* groups; 23 sequences keep the genus label). The
species-to-group mapping is in `species/database/changes.csv`. The training set
is not included in this repository; `build_database.py` downloads GTDB r95 and
builds it.
