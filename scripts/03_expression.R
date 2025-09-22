#!/usr/bin/env Rscript
# Process expression data

source("config.R")

cat("Processing expression data...\n")

# Read data
clin <- read_tsv(CLINICAL_CLEAN, show_col_types = FALSE)
counts <- read_tsv(EXPRESSION_RAW, show_col_types = FALSE)

# Prepare count matrix
gene_col <- names(counts)[1]
counts_mat <- counts %>% column_to_rownames(var = gene_col) %>% as.matrix()

# Match sample IDs
try_lengths <- 12:16
overlaps <- map_int(try_lengths, function(L) {
  length(intersect(substr(colnames(counts_mat), 1, L), clin$sample_id))
})
best_L <- try_lengths[which.max(overlaps)]
colnames(counts_mat) <- substr(colnames(counts_mat), 1, best_L)

samples_keep <- intersect(clin$sample_id, colnames(counts_mat))
clin_sub <- clin %>% filter(sample_id %in% samples_keep) %>% arrange(match(sample_id, samples_keep))
counts_mat <- counts_mat[, clin_sub$sample_id]

cat("Samples used in DESeq2:", ncol(counts_mat), "\n")

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),
  colData = clin_sub %>% column_to_rownames("sample_id"),
  design = ~ side
)

# Filter low counts
keep <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(dds))
dds <- dds[keep,]

# Run DESeq2
dds$side <- relevel(factor(dds$side), ref = "left")
dds <- DESeq(dds)

# Results
res <- results(dds, contrast = c("side","right","left"), alpha = 0.05)
res_df <- as.data.frame(res) %>% rownames_to_column("gene") %>% arrange(padj)
write_tsv(res_df, file.path(RESULTS_DIR, "DE_side_results.tsv"))

# Select DE genes
de_genes <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) >= 1)
cat("DE genes selected:", nrow(de_genes), "\n")

# VST transform and reduce
vsd <- vst(dds, blind = FALSE)
expr_vsd <- assay(vsd)

# Map Entrez IDs to symbols
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(expr_vsd),
  keytype = "ENTREZID",
  columns = c("SYMBOL")
) %>% dplyr::distinct(ENTREZID, .keep_all = TRUE)

rownames(expr_vsd) <- id_map$SYMBOL[match(rownames(expr_vsd), id_map$ENTREZID)]

# Ensure de_genes uses SYMBOL
de_genes <- de_genes %>%
  mutate(symbol = id_map$SYMBOL[match(gene, id_map$ENTREZID)]) %>%
  filter(!is.na(symbol))

expr_reduced <- expr_vsd[rownames(expr_vsd) %in% de_genes$symbol, , drop = FALSE]
expr_reduced_df <- as.data.frame(t(expr_reduced)) %>% rownames_to_column("sample_id")

write_tsv(expr_reduced_df, EXPRESSION_DE)
cat("Saved reduced expression matrix:", EXPRESSION_DE, "\n")