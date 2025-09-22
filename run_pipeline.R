#!/usr/bin/env Rscript
# Master script to run the entire pipeline

source("config.R")

cat("Starting TCGA COAD Multi-omics Analysis Pipeline\n")
cat("===============================================\n")

# Run each step
cat("\n1. Processing clinical data...\n")
source("scripts/01_clinical.R")

cat("\n2. Processing mutation data...\n")
source("scripts/02_mutation.R")

cat("\n3. Processing expression data...\n")
source("scripts/03_expression.R")

cat("\n4. Integrating data...\n")
source("scripts/04_integration.R")

cat("\n5. Running quality control...\n")
source("scripts/05_qc.R")

cat("\n6. Running classification modeling...\n")
source("scripts/06_classification.R")

cat("\n7. Running survival analysis...\n")
source("scripts/07_survival.R")

cat("\nPipeline completed successfully!\n")
cat("Results saved to:", RESULTS_DIR, "\n")