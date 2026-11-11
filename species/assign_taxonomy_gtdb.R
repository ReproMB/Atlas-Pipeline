#!/usr/bin/env Rscript
# assign_taxonomy_gtdb.R
# Classify Lactobacillaceae ASVs against the GTDB-derived custom training set
# using dada2::assignTaxonomy() — the exact method from Lebeer et al. 2023.
#
# Usage (called from subclassify_lactobacillus.py via subprocess):
#   Rscript assign_taxonomy_gtdb.R <input_fasta> <db_fasta_gz> <output_tsv> <min_boot> <threads>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    cat("Usage: Rscript assign_taxonomy_gtdb.R <input_fasta> <db_fasta_gz> <output_tsv> <min_boot> <threads>\n",
        file = stderr())
    quit(status = 1)
}

input_fasta <- args[1]
db_path     <- args[2]
output_tsv  <- args[3]
min_boot    <- as.integer(args[4])
threads     <- as.integer(args[5])

if (!requireNamespace("dada2", quietly = TRUE)) {
    cat("ERROR: dada2 package not found. Activate the qiime2 conda environment.\n",
        file = stderr())
    quit(status = 1)
}

library(dada2)

cat(sprintf("[R] Input FASTA : %s\n", input_fasta))
cat(sprintf("[R] Database    : %s\n", db_path))
cat(sprintf("[R] min_boot    : %d\n", min_boot))
cat(sprintf("[R] threads     : %d\n", threads))

# Read sequences — getSequences() returns a named character vector:
#   names  = FASTA header (Feature ID)
#   values = DNA sequence string
seqs <- getSequences(input_fasta)
cat(sprintf("[R] Sequences loaded: %d\n", length(seqs)))

if (length(seqs) == 0) {
    cat("[R] No sequences found — writing empty output.\n")
    empty <- data.frame(ASV_ID = character(0))
    write.table(empty, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    quit(status = 0)
}

# assignTaxonomy returns a matrix where rownames = raw sequences (not ASV IDs).
# outputBootstraps=TRUE returns a list: $tax (matrix) and $boot (matrix).
# tryRC=TRUE is required: GTDB r95 SSU representative sequences for several key
# vaginal Lactobacillus species (iners, gasseri, jensenii) are stored on the
# reverse-complement strand, so forward-only matching misses them entirely.
result <- assignTaxonomy(
    seqs          = seqs,
    refFasta      = db_path,
    minBoot       = min_boot,
    multithread   = threads,
    outputBootstraps = TRUE,
    tryRC         = TRUE,
    verbose       = TRUE
)

tax_mat  <- result$tax
boot_mat <- result$boot

# Map rownames (sequences) back to Feature IDs using the named vector.
seq_to_id <- setNames(names(seqs), seqs)
feature_ids <- seq_to_id[rownames(tax_mat)]

tax_df  <- as.data.frame(tax_mat,  stringsAsFactors = FALSE)
boot_df <- as.data.frame(boot_mat, stringsAsFactors = FALSE)
colnames(boot_df) <- paste0("boot_", colnames(boot_df))

output_df <- cbind(
    data.frame(ASV_ID = feature_ids, stringsAsFactors = FALSE),
    tax_df,
    boot_df
)
rownames(output_df) <- NULL

write.table(output_df, output_tsv, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat(sprintf("[R] Written: %s (%d rows)\n", output_tsv, nrow(output_df)))

# Summary per genus group
genus_counts <- table(output_df$Genus)
cat("[R] Classification summary (Genus level):\n")
for (g in sort(names(genus_counts))) {
    cat(sprintf("    %-45s  %d ASVs\n", g, genus_counts[[g]]))
}

unclassified <- sum(is.na(output_df$Genus) | output_df$Genus == "")
cat(sprintf("[R] Unclassified at genus level: %d ASVs\n", unclassified))
