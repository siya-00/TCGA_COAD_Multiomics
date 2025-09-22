# TCGA COAD Multi-omics Analysis

A reproducible pipeline for analyzing multi-omics data from TCGA COAD (Colon Adenocarcinoma) to predict tumor side and investigate survival associations.

## Overview

This pipeline integrates clinical, mutation, and gene expression data to:
1. Predict tumor side (left vs. right) using machine learning models
2. Investigate survival differences based on tumor side
3. Identify molecular features that interact with tumor side to influence survival

## Data Requirements

Place the following files in the `data/raw/` directory:
- `coad_tcga_gdc_clinical_data.tsv` - Clinical data
- `data_mutations.txt` - Mutation data
- `data_mrna_seq_read_counts.txt` - Gene expression counts

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/TCGA_COAD_Multiomics.git
cd TCGA_COAD_Multiomics