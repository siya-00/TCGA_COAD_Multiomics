# Package requirements for TCGA COAD Multi-omics Analysis

# CRAN packages
cran_packages <- c(
  "tidyverse",      # Data manipulation and visualization
  "survival",       # Survival analysis
  "survminer",      # Survival plots
  "glmnet",         # LASSO regression
  "caret",          # Machine learning framework
  "pROC",           # ROC curve analysis
  "ranger",         # Random forest implementation
  "foreach",        # Parallel processing
  "doParallel",     # Parallel backend
  "broom",          # Tidy model outputs
  "reshape2"        # Data reshaping
)

# Bioconductor packages
bioc_packages <- c(
  "DESeq2",         # Differential expression analysis
  "org.Hs.eg.db",   # Human gene annotation
  "AnnotationDbi"   # Database interface
)

# Install CRAN packages if needed
install_if_missing <- function(packages, repos = "https://cloud.r-project.org/") {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE, repos = repos)
  }
}

# Install Bioconductor packages if needed
install_bioc_if_missing <- function(packages) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    BiocManager::install(new_packages)
  }
}

# Install packages
install_if_missing(cran_packages)
install_bioc_if_missing(bioc_packages)

# Load all packages
invisible(lapply(c(cran_packages, bioc_packages), library, character.only = TRUE))

cat("All required packages are installed and loaded.\n")