#!/usr/bin/env Rscript
# Survival analysis

source("config.R")

cat("Running survival analysis...\n")

# Load data
clin_file <- read_tsv(CLINICAL_FINAL, show_col_types = FALSE)
mut_file <- read_tsv(MUTATION_FINAL, show_col_types = FALSE)
expr_file <- read_tsv(EXPRESSION_FINAL, show_col_types = FALSE)


nfold_inner <- 5
nrepeats_stability <- 100   # increase or decrease as compute allows
nthreads <- max(1, parallel::detectCores() - 1)
registerDoParallel(cores = nthreads)

# Feature reduction thresholds (practical):
top_expr_genes <- 500   # top variable genes to keep for penalized Cox
top_mut_genes  <- 100   # top mutated genes by frequency

# ------------------------------------------------------------------------------

cat("Loading files...\n")
clin_raw <- read_tsv(clin_file, show_col_types = FALSE)
mut_raw  <- read_tsv(mut_file, show_col_types = FALSE)
expr_raw <- read_tsv(expr_file, show_col_types = FALSE)

# Basic existence checks
stopifnot("sample_id" %in% names(clin_raw))
stopifnot("sample_id" %in% names(mut_raw))
stopifnot("sample_id" %in% names(expr_raw))

# Some helpful prints
cat("Clinical rows:", nrow(clin_raw), "Mut rows:", nrow(mut_raw), "Expr rows:", nrow(expr_raw), "\n")

# Copy to working variables
clin <- clin_raw
mut  <- mut_raw
expr <- expr_raw

# ---------- 0) Harmonize clinical, check time column ----------
required_cols <- c("sample_id", "time_months", "event", "side")
if (!all(required_cols %in% names(clin))) {
  cat("WARNING: clinical file missing required columns. Trying to detect alternatives...\n")
}

# Print summary of time_months (raw)
cat("\nInitial time_months summary:\n")
print(summary(clin$time_months))
cat("Count non-positive times (<=0):", sum(!is.na(clin$time_months) & clin$time_months <= 0), "\n")
cat("Count NA times:", sum(is.na(clin$time_months)), "\n")

# Attempt to recompute times from date columns if time_months has many non-positive or NA values.
# We will look for possible date columns: any column with 'diagnos', 'death', 'last' in name (case-insensitive).
date_cols <- names(clin)[str_detect(names(clin), regex("diagnos|death|last|contact|date", ignore_case = TRUE))]
cat("Potential date-like columns in clinical:", paste(date_cols, collapse = ", "), "\n")

# helper to parse dates: detect if column looks like YYYY-MM-DD or similar
is_date_like <- function(x) {
  if (!is.character(x) && !is.factor(x)) return(FALSE)
  s <- as.character(na.omit(x))[1]
  if (is.na(s)) return(FALSE)
  return(str_detect(s, "\\d{4}-\\d{2}-\\d{2}") | str_detect(s, "/"))
}

# Try recompute if we can find diagnosis and death/last-contact columns
diag_col <- NULL; event_date_col <- NULL
possible_diag <- names(clin)[str_detect(names(clin), regex("diagnos", ignore_case = TRUE))]
possible_death <- names(clin)[str_detect(names(clin), regex("death|last.*contact|last contact|last_contact|lastcontact", ignore_case = TRUE))]

# choose columns that are date-like
if (length(possible_diag) > 0) {
  for (c in possible_diag) if (is_date_like(clin[[c]])) { diag_col <- c; break }
}
if (length(possible_death) > 0) {
  for (c in possible_death) if (is_date_like(clin[[c]])) { event_date_col <- c; break }
}

if (!is.null(diag_col) && !is.null(event_date_col)) {
  cat("Found date columns for recomputing time:", diag_col, "and", event_date_col, "\n")
  # try to coerce to Date
  clin <- clin %>%
    mutate(
      dx_date = as.Date(.data[[diag_col]]),
      ev_date = as.Date(.data[[event_date_col]])
    )
  # compute days and convert to months
  if (all(!is.na(clin$dx_date))) {
    clin <- clin %>%
      mutate(time_days_calc = as.numeric(difftime(ev_date, dx_date, units = "days")),
             time_months_calc = time_days_calc / 30.44)
    # if recomputed values look reasonable, use them for NA or non-positive originals
    bad_before <- which(is.na(clin$time_months) | clin$time_months <= 0)
    if (length(bad_before) > 0) {
      use_idx <- which(!is.na(clin$time_months_calc) & clin$time_months_calc > 0)
      n_replace <- sum(use_idx %in% bad_before)
      if (n_replace > 0) {
        clin$time_months[use_idx] <- clin$time_months_calc[use_idx]
        cat("Replaced", n_replace, "bad time_months values with computed time_months_calc.\n")
      }
    }
  } else {
    cat("Recomputed diagnosis dates contained NA -> skipping automated recompute.\n")
  }
}

# After recompute attempt, re-check times
cat("Post-recompute time_months summary:\n")
print(summary(clin$time_months))
cat("Count non-positive (<=0):", sum(!is.na(clin$time_months) & clin$time_months <= 0), "\n")
cat("Count NA:", sum(is.na(clin$time_months)), "\n")

# If still any non-positive or NA times, drop those rows (but log them)
bad_idx <- which(is.na(clin$time_months) | clin$time_months <= 0)
if (length(bad_idx) > 0) {
  cat("Dropping", length(bad_idx), "samples with NA or non-positive survival time. Writing list to disk.\n")
  write_tsv(clin[bad_idx, c("sample_id", "time_months")], file.path(out_dir, "dropped_samples_bad_time.tsv"))
  keep_ids <- setdiff(clin$sample_id, clin$sample_id[bad_idx])
  clin <- clin %>% filter(sample_id %in% keep_ids)
  mut  <- mut  %>% filter(sample_id %in% keep_ids)
  expr <- expr %>% filter(sample_id %in% keep_ids)
  cat("Remaining samples:", nrow(clin), "\n")
} else {
  cat("All survival times positive and non-NA. Proceeding.\n")
}

# final basic clinical checks
stopifnot(all(c("sample_id","time_months","event","side") %in% names(clin)))
clin$event <- as.integer(clin$event)   # ensure numeric 0/1
cat("Events (1) count:", sum(clin$event == 1, na.rm = TRUE), "Censored:", sum(clin$event == 0, na.rm = TRUE), "\n")

# ---------- 1) Keep only samples present in all datasets ----------
common_ids <- Reduce(intersect, list(clin$sample_id, mut$sample_id, expr$sample_id))
cat("Samples present in all 3 data types:", length(common_ids), "\n")
clin <- clin %>% filter(sample_id %in% common_ids) %>% arrange(match(sample_id, common_ids))
mut  <- mut  %>% filter(sample_id %in% common_ids)  %>% arrange(match(sample_id, common_ids))
expr <- expr %>% filter(sample_id %in% common_ids) %>% arrange(match(sample_id, common_ids))

# safety
stopifnot(all(clin$sample_id == mut$sample_id), all(clin$sample_id == expr$sample_id))

# ---------- 2) Kaplan-Meier by side ----------
cat("Fitting Kaplan-Meier by side...\n")
surv_obj <- Surv(time = clin$time_months, event = clin$event)
fit_km <- survfit(surv_obj ~ side, data = clin)
p <- ggsurvplot(fit_km, data = clin, pval = TRUE, risk.table = TRUE,
                legend.title = "Side", palette = c("#1f77b4","#d62728"))
ggsave(file.path(out_dir, "KM_by_side.png"), plot = p$plot, width = 7, height = 6)

# ---------- 3) Baseline Cox (side + basic covariates if present) ----------
# pick covariates if present
covars <- c()
if ("Diagnosis Age" %in% names(clin)) {
  clin <- clin %>% mutate(Age = as.numeric(`Diagnosis Age`))
  covars <- c(covars, "Age")
} 
if ("AJCC Pathologic Stage" %in% names(clin)) {
  clin <- clin %>% mutate(Stage = factor(`AJCC Pathologic Stage`))
  covars <- c(covars, "Stage")
}
if ("Sex" %in% names(clin)) {
  clin <- clin %>% mutate(Sex = factor(Sex))
  covars <- c(covars, "Sex")
}
# build formula
base_formula <- as.formula(paste0("Surv(time_months, event) ~ side", if (length(covars)>0) paste0(" + ", paste(covars, collapse = " + ")) else ""))
cat("Baseline Cox formula:", deparse(base_formula), "\n")
cox_base <- tryCatch(coxph(base_formula, data = clin), error = function(e) { e })
if (inherits(cox_base, "error")) {
  cat("Baseline Cox error:", cox_base$message, "\n")
} else {
  cox_tab <- broom::tidy(cox_base, exponentiate = TRUE, conf.int = TRUE)
  write_tsv(cox_tab, file.path(out_dir, "cox_baseline_table.tsv"))
  # save cox zph results for top terms
  zph <- tryCatch(cox.zph(cox_base), error = function(e) NULL)
  if (!is.null(zph)) {
    capture.output(zph, file = file.path(out_dir, "cox_baseline_zph.txt"))
  }
}

# ---------- 4) Feature-wise Cox with side interaction (limited candidate set) ----------
cat("Preparing candidate feature lists for univariate interaction screening...\n")
# Mutation candidates: pick driver genes if present, else top mutated genes by frequency
drivers <- c("TP53","KRAS","BRAF","APC","PIK3CA","SMAD4","FBXW7","NRAS","ATM")
mut_genes_present <- setdiff(colnames(mut), "sample_id")
driver_genes <- intersect(drivers, mut_genes_present)
if (length(driver_genes) == 0) {
  # choose top N mutated by frequency
  mut_freq <- colSums(mut[ , mut_genes_present, drop = FALSE], na.rm = TRUE)
  driver_genes <- names(sort(mut_freq, decreasing = TRUE))[1:min(top_mut_genes, length(mut_freq))]
}
cat("Mutation candidate genes (n):", length(driver_genes), "\n")

# Expression candidates: top variable genes
expr_mat <- expr %>% column_to_rownames("sample_id")
expr_var <- apply(expr_mat, 2, var, na.rm = TRUE)
top_expr <- names(sort(expr_var, decreasing = TRUE))[1:min(top_expr_genes, length(expr_var))]
cat("Expression candidates (n):", length(top_expr), "\n")

# Run featurewise Cox with side:feat interaction
feature_results <- list()

# Mutations: binary variable
for (g in driver_genes) {
  cat("Testing mutation feature:", g, "\n")
  df <- clin %>% mutate(feat = as.numeric(mut[[g]]))
  # skip if feature constant
  if (length(unique(df$feat)) <= 1) {
    cat(" - skipped (constant)\n"); next
  }
  # fit model with interaction
  f <- tryCatch(coxph(Surv(time_months, event) ~ side * feat + ., data = df %>% select(time_months, event, side, feat, all_of(covars))), error = function(e) e)
  if (inherits(f,"error")) {
    cat(" - coxph error:", f$message, "\n")
    next
  }
  tb <- broom::tidy(f, exponentiate = TRUE, conf.int = TRUE) %>% mutate(feature = g, type = "mutation")
  feature_results[[paste0("mut_", g)]] <- tb
}

# Expression features: continuous
for (g in top_expr) {
  cat("Testing expression feature:", g, "\n")
  df <- clin %>% mutate(feat = as.numeric(expr_mat[clin$sample_id, g]))
  if (all(is.na(df$feat))) { cat(" - skipped (NA)\n"); next }
  # standardize feature for interpretability
  df$feat <- scale(df$feat)[,1]
  f <- tryCatch(coxph(Surv(time_months, event) ~ side * feat + ., data = df %>% select(time_months, event, side, feat, all_of(covars))), error = function(e) e)
  if (inherits(f,"error")) { cat(" - coxph error:", f$message, "\n"); next }
  tb <- broom::tidy(f, exponentiate = TRUE, conf.int = TRUE) %>% mutate(feature = g, type = "expression")
  feature_results[[paste0("expr_", g)]] <- tb
}

feature_tab <- bind_rows(feature_results)
write_tsv(feature_tab, file.path(out_dir, "featurewise_cox_with_interaction.tsv"))
cat("Featurewise Cox written to:", file.path(out_dir, "featurewise_cox_with_interaction.tsv"), "\n")

# ---------- 5) Penalized multivariable Cox (LASSO) on reduced feature set ----------
# We will reduce features to:
#  - expression: top_expr (top variable genes)
#  - mutation: top X genes by frequency (top_mut_genes)
cat("Preparing features for penalized Cox...\n")
# Clean expression matrix: select top_expr and scale
expr_keep <- expr_mat[clin$sample_id, top_expr, drop = FALSE]
expr_keep_scaled <- scale(expr_keep)

# Mutation: choose top by frequency
mut_mat <- mut %>% column_to_rownames("sample_id")
mut_freq_all <- colSums(mut_mat, na.rm = TRUE)
mut_top <- names(sort(mut_freq_all, decreasing = TRUE))[1:min(top_mut_genes, length(mut_freq_all))]
mut_keep <- mut_mat[clin$sample_id, mut_top, drop = FALSE]

cat("Final dims expr:", dim(expr_keep_scaled), "mut:", dim(mut_keep), "\n")

# Clean columns: remove NA cols and constant cols
clean_mat <- function(M, name="mat") {
  # Remove NA-containing columns
  na_cols <- which(apply(M, 2, function(x) any(is.na(x))))
  if (length(na_cols) > 0) {
    cat(sprintf(" - removing %d NA columns from %s\n", length(na_cols), name))
    M <- M[, -na_cols, drop = FALSE]
  }
  const_cols <- which(apply(M, 2, function(x) length(unique(x[!is.na(x)])) <= 1))
  if (length(const_cols) > 0) {
    cat(sprintf(" - removing %d constant columns from %s\n", length(const_cols), name))
    M <- M[, -const_cols, drop = FALSE]
  }
  return(M)
}
expr_clean <- clean_mat(expr_keep_scaled, "expr")
mut_clean  <- clean_mat(as.matrix(mut_keep), "mut")

# Combine
X_all <- cbind(expr_clean, mut_clean)
cat("Combined feature matrix dims:", dim(X_all), "\n")
# Sanity checks
if (any(is.na(X_all))) stop("NA present in X_all — handle before glmnet.")
if (min(clin$time_months, na.rm = TRUE) <= 0) stop("Non-positive survival times remain — aborting.")

# If features >> events, warn and optionally reduce further via univariate Cox filter
n_events <- sum(clin$event == 1, na.rm = TRUE)
n_features <- ncol(X_all)
cat("Events:", n_events, "Features:", n_features, "\n")
if (n_features > 5 * n_events) {
  cat("WARNING: number of features is large relative to events. Applying univariate Cox filter (p<0.1) to reduce features.\n")
  pvals <- rep(NA_real_, n_features)
  names(pvals) <- colnames(X_all)
  for (i in seq_len(n_features)) {
    feat <- colnames(X_all)[i]
    df <- data.frame(time_months = clin$time_months, event = clin$event, feat = X_all[,i])
    f <- tryCatch(coxph(Surv(time_months, event) ~ feat + ., data = df), error = function(e) NULL)
    if (is.null(f)) { pvals[feat] <- NA; next }
    tb <- tryCatch(broom::tidy(f), error = function(e) NULL)
    if (is.null(tb)) pvals[feat] <- NA else {
      # take p-value for feat (first row after intercept)
      pvals[feat] <- tb$p.value[ which(tb$term == "feat") ][1]
    }
  }
  keep_small <- names(pvals)[!is.na(pvals) & pvals < 0.1]
  if (length(keep_small) < 10) {
    # ensure at least some features (take top 50 by smallest p)
    keep_small <- names(sort(pvals, na.last = NA))[1:min(50, length(na.omit(pvals)))]
  }
  cat("Reducing features from", n_features, "to", length(keep_small), "by univariate p<0.1\n")
  X_all <- X_all[, keep_small, drop = FALSE]
  n_features <- ncol(X_all)
}

# Prepare Surv object and fit cv.glmnet for Cox
surv_y <- Surv(time = clin$time_months, event = clin$event)
set.seed(2025)
cat("Running cv.glmnet (Cox) with", nfold_inner, "folds. This can take time...\n")
cvfit <- cv.glmnet(x = as.matrix(X_all), y = surv_y, family = "cox", nfolds = nfold_inner, parallel = TRUE)
cat("cv.glmnet completed. lambda.min =", cvfit$lambda.min, "\n")

fit_coxnet <- glmnet(x = as.matrix(X_all), y = surv_y, family = "cox", lambda = cvfit$lambda.min)
coefs <- coef(fit_coxnet)
sel_idx <- which(as.numeric(coefs) != 0)
selected_features <- rownames(coefs)[sel_idx]
cat("Selected features (non-zero):", length(selected_features), "\n")
writeLines(selected_features, file.path(out_dir, "penalized_cox_selected_features.txt"))
write_tsv(tibble(feature = selected_features), file.path(out_dir, "penalized_cox_selected_features.tsv"))

# Save coefficients table (sparse)
coef_df <- tibble(feature = rownames(coefs), coef = as.numeric(coefs)) %>% filter(coef != 0)
write_tsv(coef_df, file.path(out_dir, "penalized_cox_coefficients.tsv"))

# ---------- 6) Stability selection (subsampling) for LASSO -->
cat("Running LASSO stability selection with", nrepeats_stability, "subsamples (parallel)...\n")
set.seed(2025)
stability_res <- foreach(rep = 1:nrepeats_stability, .combine = bind_rows, .packages = c("glmnet")) %dopar% {
  # subsample 70% of samples
  samp <- sample(seq_len(nrow(X_all)), size = floor(0.7 * nrow(X_all)), replace = FALSE)
  Xs <- X_all[samp, , drop = FALSE]
  ys <- Surv(time = clin$time_months[samp], event = clin$event[samp])
  # small check: remove constant cols in this subsample
  const_cols <- which(apply(Xs, 2, function(z) length(unique(z[!is.na(z)])) <= 1))
  if (length(const_cols) > 0) Xs <- Xs[, -const_cols, drop = FALSE]
  # cv and fit
  cvg <- tryCatch(cv.glmnet(x = as.matrix(Xs), y = ys, family = "cox", nfolds = min(5, floor(nrow(Xs)/2))), error = function(e) NULL)
  if (is.null(cvg)) return(tibble(rep = rep, feature = character()))
  lam <- cvg$lambda.min
  fit <- glmnet(x = as.matrix(Xs), y = ys, family = "cox", lambda = lam)
  cfs <- as.matrix(coef(fit))
  sel <- rownames(cfs)[which(cfs[,1] != 0)]
  if (length(sel) == 0) return(tibble(rep = rep, feature = character()))
  tib <- tibble(rep = rep, feature = sel)
  tib
}

stability_summary <- stability_res %>%
  group_by(feature) %>%
  summarise(freq = n(), frac = n() / nrepeats_stability) %>%
  arrange(desc(freq))

write_tsv(stability_summary, file.path(out_dir, "lasso_stability_summary.tsv"))
cat("Wrote stability summary to:", file.path(out_dir, "lasso_stability_summary.tsv"), "\n")

# ---------- 7) Diagnostics & proportional hazards check for top selected features ----------
if (length(selected_features) > 0) {
  top_for_check <- selected_features[1:min(20, length(selected_features))]
  diag_tab <- list()
  for (f in top_for_check) {
    # build a Cox model with feature + side + covariates
    dft <- data.frame(time_months = clin$time_months, event = clin$event, feat = X_all[, f], side = clin$side)
    if (length(covars) > 0) {
      # attempt to bind covariates from clin
      for (cv in covars) dft[[cv]] <- clin[[cv]]
    }
    mm <- tryCatch(coxph(Surv(time_months, event) ~ feat + side + ., data = dft), error = function(e) NULL)
    if (!is.null(mm)) {
      zph <- tryCatch(cox.zph(mm), error = function(e) NULL)
      diag_tab[[f]] <- list(cox = mm, zph = zph)
    } else {
      diag_tab[[f]] <- list(cox = NULL, zph = NULL)
    }
  }
  # Save brief zph outputs
  zph_out <- lapply(names(diag_tab), function(nm) {
    z <- diag_tab[[nm]]$zph
    if (is.null(z)) return(tibble(feature = nm, term = NA, chisq = NA, p = NA))
    dfz <- as_tibble(z$table, rownames = "term")
    dfz <- dfz %>% mutate(feature = nm) %>% select(feature, everything())
    dfz
  }) %>% bind_rows()
  write_tsv(zph_out, file.path(out_dir, "proportional_hazards_zph_top_features.tsv"))
}

# ---------- Save session info and small report ----------
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_stepC.txt"))
report_lines <- c(
  paste0("Step C completed at: ", Sys.time()),
  paste0("Samples used: ", nrow(clin)),
  paste0("Events (deaths): ", sum(clin$event == 1, na.rm = TRUE)),
  paste0("Final feature count for penalized Cox: ", ncol(X_all)),
  paste0("Selected penalized-cox features (n): ", length(selected_features))
)
writeLines(report_lines, file.path(out_dir, "summary_stepC.txt"))

cat("Step C pipeline finished. Outputs in:", out_dir, "\n")

