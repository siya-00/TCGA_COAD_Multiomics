#!/usr/bin/env Rscript
# Process mutation data

source("config.R")

cat("Processing mutation data...\n")

# Read clinical and mutation data
clin <- read_tsv(CLINICAL_CLEAN, show_col_types = FALSE)
mut <- read_tsv(MUTATION_RAW, comment = "#", show_col_types = FALSE)

cat("Clinical samples:", nrow(clin), "\n")
cat("Raw mutation rows:", nrow(mut), "\n")

# Filter mutations
nonsyn_classes <- c(
  "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
  "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation"
)

mut_filt <- mut %>%
  filter(
    is.na(GDC_FILTER),
    Variant_Classification %in% nonsyn_classes,
    Mutation_Status == "Somatic",
    is.na(MAX_AF) | MAX_AF < 0.01,
    is.na(gnomAD_AF) | gnomAD_AF < 0.01,
    is.na(`1000G_AF`) | `1000G_AF` < 0.01,
    !is.na(t_alt_count) & t_alt_count >= 20
  )

# Match sample IDs
try_lengths <- 12:16
overlaps <- map_int(try_lengths, function(L) {
  length(intersect(substr(mut_filt$Tumor_Sample_Barcode, 1, L), clin$sample_id))
})
best_L <- try_lengths[which.max(overlaps)]
cat("Best substring length:", best_L, "\n")

mut_filt <- mut_filt %>%
  mutate(sample_id = substr(Tumor_Sample_Barcode, 1, best_L))

# Build binary matrix
mut_bin <- mut_filt %>%
  distinct(sample_id, Hugo_Symbol) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = value, values_fill = 0)

# Align with clinical samples and filter rare genes
mut_bin <- clin %>%
  select(sample_id) %>%
  left_join(mut_bin, by = "sample_id") %>%
  mutate(across(-sample_id, ~replace_na(.x, 0)))

min_cnt <- ceiling(MIN_MUTATION_FREQ * nrow(mut_bin))
gene_counts <- sort(colSums(mut_bin[,-1], na.rm = TRUE), decreasing = TRUE)
keep_genes <- names(gene_counts)[gene_counts >= min_cnt]

cat("Keeping genes mutated in >=", min_cnt, "samples (", MIN_MUTATION_FREQ*100, "%).\n")
cat("Original gene count:", length(gene_counts), "-> kept:", length(keep_genes), "\n")

mut_bin_final <- mut_bin %>% select(sample_id, all_of(keep_genes))

# Save
write_tsv(mut_bin_final, MUTATION_BINARY)

cat("Saved mutation matrix:", MUTATION_BINARY, "\n")
cat("Dimensions:", dim(mut_bin_final), "\n")
cat("Total mutated entries:", sum(mut_bin_final[,-1]), "\n")