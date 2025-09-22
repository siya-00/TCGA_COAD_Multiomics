#!/usr/bin/env Rscript
# Quality control and exploratory analysis

source("config.R")

cat("Running quality control...\n")

# Load data
clin <- read_tsv(CLINICAL_FINAL, show_col_types = FALSE)
mut <- read_tsv(MUTATION_FINAL, show_col_types = FALSE)
expr <- read_tsv(EXPRESSION_FINAL, show_col_types = FALSE)

# Dimensions
qc_summary <- tibble(
  dataset = c("Clinical", "Mutation", "Expression"),
  n_samples = c(nrow(clin), nrow(mut), nrow(expr)),
  n_features = c(ncol(clin), ncol(mut)-1, ncol(expr)-1)
)
write_tsv(qc_summary, file.path(RESULTS_DIR, "qc", "qc_summary_dimensions.tsv"))

# Survival analysis
fit <- survfit(Surv(time_months, event) ~ side, data = clin)
cox <- coxph(Surv(time_months, event) ~ side, data = clin)
cox_sum <- broom::tidy(cox, exponentiate = TRUE, conf.int = TRUE)
write_tsv(cox_sum, file.path(RESULTS_DIR, "qc", "qc_survival_cox.tsv"))

p1 <- ggsurvplot(
  fit, data = clin, pval = TRUE, risk.table = TRUE,
  palette = c("firebrick2", "steelblue3"),
  legend.title = "Tumor Side",
  legend.labs = c("Left", "Right"),
  xlab = "Time (months)", ylab = "Overall survival probability",
  title = "Overall Survival by Tumor Side"
)
ggsave(file.path(RESULTS_DIR, "qc", "qc_survival_curve.png"), p1$plot, width = 8, height = 6, dpi = 300)

# Mutation burden
mut_counts <- rowSums(mut[,-1])
clin$mut_burden <- mut_counts

p2 <- ggplot(clin, aes(x = side, y = mut_burden, fill = side)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, fill = "white") +
  scale_fill_manual(values = c("firebrick2", "steelblue3")) +
  theme_classic() +
  labs(title = "Mutation burden by tumor side", x = "Tumor side", y = "Mutated genes per sample")
ggsave(file.path(RESULTS_DIR, "qc", "qc_mutation_burden.png"), p2, width = 7, height = 5, dpi = 300)

top_genes <- sort(colSums(mut[,-1]), decreasing = TRUE)[1:20]
write_tsv(enframe(top_genes, name = "gene", value = "n_samples"), 
          file.path(RESULTS_DIR, "qc", "qc_top_mutated_genes.tsv"))

# Expression PCA
expr_matrix <- expr[,-1] %>% as.matrix()
expr_pca <- prcomp(expr_matrix, scale. = TRUE)
var_exp <- round(100 * (expr_pca$sdev^2 / sum(expr_pca$sdev^2))[1:2], 1)

pca_df <- data.frame(
  PC1 = expr_pca$x[,1], PC2 = expr_pca$x[,2], side = clin$side
)
p3 <- ggplot(pca_df, aes(PC1, PC2, color = side)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("firebrick2", "steelblue3")) +
  theme_classic() +
  labs(
    title = "PCA of Expression (DE genes)",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)"),
    color = "Tumor side"
  )
ggsave(file.path(RESULTS_DIR, "qc", "qc_expression_pca.png"), p3, width = 7, height = 6, dpi = 300)

cat("QC complete. Saved plots + summary tables.\n")