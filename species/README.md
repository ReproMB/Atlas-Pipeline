# *Lactobacillus* species groups

SILVA resolves *Lactobacillus* only to genus. This step assigns the ASVs that
SILVA places in Lactobacillaceae to species groups, following the
amplicon-taxonomy approach of Lebeer et al. (2023), *Nat Microbiol* 8:2183–2195
(doi:10.1038/s41564-023-01500-0). It runs after the pipeline, on the output
folder of each processing unit.

1. `subclassify_lactobacillus.py` classifies the Lactobacillaceae ASVs with
   dada2 `assignTaxonomy` (`assign_taxonomy_gtdb.R`; minBoot 50, `tryRC = TRUE`).
   ASVs below the bootstrap threshold stay as *Lactobacillus* without a group.
2. `consensus.py` combines 10 runs: an ASV keeps its group only if it is the
   same in all of them.
3. `export_resolved_table.py` writes `exports/<unit>_resolved_genus_table.tsv`,
   the genus table with *Lactobacillus* split into groups. Every other genus is
   unchanged, and the group counts add up to the original *Lactobacillus* count.

`run_species_step.sh` runs the three in sequence:

```bash
species/run_species_step.sh results/<unit> <unit> trainingset_Lactobacillaceae_GTDB_05-RS95.fasta.gz
```

The training set is described in [`../docs/reference_database.md`](../docs/reference_database.md).
The step needs Python with biom-format and pandas, and Rscript with dada2
(both in the `qiime2-amplicon-2025.7` conda environment).
