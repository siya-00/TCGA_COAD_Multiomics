#!/usr/bin/env Rscript
# Process clinical data

source("config.R")

cat("Processing clinical data...\n")

# Read raw clinical data
clin <- read_tsv(CLINICAL_RAW, comment = "#", show_col_types = FALSE)
cat("Raw clinical rows:", nrow(clin), "\n")

# Process clinical data
clin_clean <- clin %>%
  filter(`Sample Type` == "Primary Tumor") %>%
  drop_na(`Sample ID`, `Patient ID`) %>%
  mutate(
    time_months = suppressWarnings(as.numeric(`Overall Survival (Months)`)),
    event = case_when(
      str_detect(`Overall Survival Status`, "DECEASED|Dead|1") ~ 1,
      str_detect(`Overall Survival Status`, "LIVING|Alive|0") ~ 0,
      TRUE ~ NA_real_
    ),
    side = case_when(
      `ICD-10 Classification` %in% c("C18.0","C18.2","C18.3","C18.4") ~ "right",
      `ICD-10 Classification` %in% c("C18.5","C18.6","C18.7") ~ "left",
      TRUE ~ NA_character_
    ),
    patient_id = substr(`Patient ID`, 1, 12),
    sample_id  = `Sample ID`
  ) %>%
  drop_na(time_months, event, side, patient_id, sample_id)

# Save cleaned clinical data
write_tsv(clin_clean, CLINICAL_CLEAN)

cat("Saved clinical_clean.tsv with", nrow(clin_clean), "samples\n")
cat("Side distribution:\n")
print(table(clin_clean$side))
cat("Event distribution:\n")
print(table(clin_clean$event))