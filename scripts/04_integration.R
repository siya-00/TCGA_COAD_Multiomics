#!/usr/bin/env Rscript
# Integrate all data types

source("config.R")

cat("Integrating data...\n")

# Read all data
clin <- read_tsv(CLINICAL_CLEAN, show_col_types = FALSE)
mut <- read_tsv(MUTATION_BINARY, show_col_types = FALSE)
expr <- read_tsv(EXPRESSION_DE, show_col_types = FALSE)

# Find common samples
common <- Reduce(intersect, list(clin$sample_id, mut$sample_id, expr$sample_id))
cat("Final common samples:", length(common), "\n")

# Filter to common samples
clin_final <- clin %>% filter(sample_id %in% common)
mut_final <- mut %>% filter(sample_id %in% common)
expr_final <- expr %>% filter(sample_id %in% common)

# Save final datasets
write_tsv(clin_final, CLINICAL_FINAL)
write_tsv(mut_final, MUTATION_FINAL)
write_tsv(expr_final, EXPRESSION_FINAL)

cat("Saved intersected clinical, mutation, and expression files.\n")