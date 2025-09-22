# Configuration file for TCGA COAD Multi-omics Analysis

# Paths
DATA_DIR <- "data"
RAW_DIR <- file.path(DATA_DIR, "raw")
PROCESSED_DIR <- file.path(DATA_DIR, "processed")
RESULTS_DIR <- "results"

# Create directories if they don't exist
dir.create(DATA_DIR, showWarnings = FALSE)
dir.create(RAW_DIR, showWarnings = FALSE)
dir.create(PROCESSED_DIR, showWarnings = FALSE)
dir.create(RESULTS_DIR, showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "qc"), showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "classification"), showWarnings = FALSE)
dir.create(file.path(RESULTS_DIR, "survival"), showWarnings = FALSE)

# File names
CLINICAL_RAW <- file.path(RAW_DIR, "coad_tcga_gdc_clinical_data.tsv")
MUTATION_RAW <- file.path(RAW_DIR, "data_mutations.txt")
EXPRESSION_RAW <- file.path(RAW_DIR, "data_mrna_seq_read_counts.txt")

CLINICAL_CLEAN <- file.path(PROCESSED_DIR, "clinical_clean.tsv")
MUTATION_BINARY <- file.path(PROCESSED_DIR, "x1_mutation_binary.tsv")
EXPRESSION_DE <- file.path(PROCESSED_DIR, "x2_expr_DE_reduced_vsd.tsv")

CLINICAL_FINAL <- file.path(PROCESSED_DIR, "clinical_final.tsv")
MUTATION_FINAL <- file.path(PROCESSED_DIR, "x1_mutation_final.tsv")
EXPRESSION_FINAL <- file.path(PROCESSED_DIR, "x2_expr_final.tsv")

# Analysis parameters
SEED <- 2025
N_CORES <- parallel::detectCores() - 1
N_FOLDS_OUTER <- 5
N_FOLDS_INNER <- 5
N_REPEATS_STABILITY <- 50
TOP_EXPR_GENES <- 500
TOP_MUT_GENES <- 100
MIN_MUTATION_FREQ <- 0.05

# Load required packages
required_packages <- c(
  "tidyverse", "survival", "survminer", "glmnet", "caret", 
  "pROC", "ranger", "foreach", "doParallel", "DESeq2", 
  "org.Hs.eg.db", "AnnotationDbi", "broom", "reshape2"
)

# Install missing packages
install_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
  # Load all packages
  suppressPackageStartupMessages(invisible(lapply(packages, library, character.only = TRUE)))
}

install_missing(required_packages)

select <- dplyr::select  # Make sure dplyr select is used by default

# Set seed and parallel backend
set.seed(SEED)
registerDoParallel(cores = N_CORES)