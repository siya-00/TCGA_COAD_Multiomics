#!/usr/bin/env Rscript
# Classification modeling

source("config.R")

cat("Running classification modeling...\n")

# Load data
clin <- read_tsv(CLINICAL_FINAL, show_col_types = FALSE)
mut <- read_tsv(MUTATION_FINAL, show_col_types = FALSE)
expr <- read_tsv(EXPRESSION_FINAL, show_col_types = FALSE)

# Prepare data
y <- factor(clin$side)
X_expr <- expr %>% column_to_rownames("sample_id") %>% as.matrix()
X_mut <- mut %>% column_to_rownames("sample_id") %>% as.matrix()

# Clinical covariates
clin_cov <- clin %>%
  select(sample_id, Age = `Diagnosis Age`, Stage = `AJCC Pathologic Stage`, Sex = Sex) %>%
  mutate(
    Age = as.numeric(Age),
    Stage = factor(Stage),
    Sex = factor(Sex)
  )

# Create stratified folds
set.seed(SEED)
folds <- createFolds(y, k = N_FOLDS_OUTER, list = TRUE, returnTrain = FALSE)

# Store results
results_summary <- tibble()
roc_list <- list()
lasso_feature_list <- list()
rf_importance_list <- list()

# Outer CV loop
for (fold_i in seq_along(folds)) {
  cat("\n====================\nOuter fold", fold_i, "of", length(folds), "\n====================\n")
  
  test_idx <- folds[[fold_i]]
  train_idx <- setdiff(seq_along(y), test_idx)
  
  # Prepare training and test data
  sample_ids_train <- clin$sample_id[train_idx]
  sample_ids_test <- clin$sample_id[test_idx]
  
  Xexpr_train <- X_expr[sample_ids_train, , drop = FALSE]
  Xexpr_test <- X_expr[sample_ids_test, , drop = FALSE]
  Xmut_train <- X_mut[sample_ids_train, , drop = FALSE]
  Xmut_test <- X_mut[sample_ids_test, , drop = FALSE]
  
  clin_train <- clin_cov %>% filter(sample_id %in% sample_ids_train) %>% arrange(match(sample_id, sample_ids_train))
  clin_test <- clin_cov %>% filter(sample_id %in% sample_ids_test) %>% arrange(match(sample_id, sample_ids_test))
  
  # Scale expression data using training stats
  expr_means <- colMeans(Xexpr_train, na.rm = TRUE)
  expr_sds <- apply(Xexpr_train, 2, sd, na.rm = TRUE)
  expr_sds[expr_sds == 0] <- 1
  
  Xexpr_train_s <- sweep(Xexpr_train, 2, expr_means, "-")
  Xexpr_train_s <- sweep(Xexpr_train_s, 2, expr_sds, "/")
  Xexpr_test_s <- sweep(Xexpr_test, 2, expr_means, "-")
  Xexpr_test_s <- sweep(Xexpr_test_s, 2, expr_sds, "/")
  
  # Prepare response variables
  y_train_num <- ifelse(y[train_idx] == "right", 1, 0)
  y_test_num <- ifelse(y[test_idx] == "right", 1, 0)
  
  # Expression-only LASSO
  cv_glm <- cv.glmnet(Xexpr_train_s, y_train_num, family = "binomial",
                      nfolds = N_FOLDS_INNER, type.measure = "auc", parallel = TRUE)
  fit_lasso <- glmnet(Xexpr_train_s, y_train_num, family = "binomial",
                      lambda = cv_glm$lambda.min, standardize = FALSE)
  preds_prob <- predict(fit_lasso, newx = Xexpr_test_s, type = "response")[,1]
  roc_obj <- roc(y_test_num, preds_prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  results_summary <- bind_rows(
    results_summary,
    tibble(model = "expr_lasso", fold = fold_i, auc = auc_val,
           n_selected = sum(coef(fit_lasso)[-1,1] != 0))
  )
  roc_list[[paste0("expr_lasso_fold", fold_i)]] <- roc_obj
  
  # Save selected features
  coefs <- as.matrix(coef(fit_lasso))
  sel_features <- rownames(coefs)[which(coefs[,1] != 0)]
  sel_features <- sel_features[sel_features != "(Intercept)"]
  lasso_feature_list[[paste0("expr_fold", fold_i)]] <- tibble(fold = fold_i, feature = sel_features)
  
  # Mutation-only LASSO
  keep_mut_cols <- which(apply(Xmut_train, 2, function(col) length(unique(col)) > 1))
  Xmut_train2 <- if (length(keep_mut_cols) > 0) Xmut_train[, keep_mut_cols, drop = FALSE] else matrix(nrow = nrow(Xmut_train), ncol = 0)
  Xmut_test2 <- if (length(keep_mut_cols) > 0) Xmut_test[, keep_mut_cols, drop = FALSE] else matrix(nrow = nrow(Xmut_test), ncol = 0)
  
  if (ncol(Xmut_train2) > 0) {
    cv_glm_mut <- cv.glmnet(Xmut_train2, y_train_num, family = "binomial",
                            nfolds = N_FOLDS_INNER, type.measure = "auc", parallel = TRUE)
    fit_lasso_mut <- glmnet(Xmut_train2, y_train_num, family = "binomial",
                            lambda = cv_glm_mut$lambda.min, standardize = FALSE)
    preds_prob_mut <- predict(fit_lasso_mut, newx = Xmut_test2, type = "response")[,1]
    roc_obj_mut <- roc(y_test_num, preds_prob_mut, quiet = TRUE)
    auc_val_mut <- as.numeric(auc(roc_obj_mut))
    
    results_summary <- bind_rows(
      results_summary,
      tibble(model = "mut_lasso", fold = fold_i, auc = auc_val_mut,
             n_selected = sum(coef(fit_lasso_mut)[-1,1] != 0))
    )
    roc_list[[paste0("mut_lasso_fold", fold_i)]] <- roc_obj_mut
  } else {
    cat("Skipping mutation LASSO - no variable mutation columns\n")
  }
  
  # Combined LASSO (expression + mutation + clinical)
  # Prepare clinical features
  clin_train_cov <- clin_train %>% select(-sample_id)
  clin_test_cov <- clin_test %>% select(-sample_id)
  
  if (ncol(clin_train_cov) > 0) {
    dv <- dummyVars(~ ., data = clin_train_cov, fullRank = FALSE)
    mm_train <- predict(dv, newdata = clin_train_cov)
    mm_test <- predict(dv, newdata = clin_test_cov)
    mm_train[is.na(mm_train)] <- 0
    mm_test[is.na(mm_test)] <- 0
  } else {
    mm_train <- matrix(nrow = nrow(Xexpr_train_s), ncol = 0)
    mm_test <- matrix(nrow = nrow(Xexpr_test_s), ncol = 0)
  }
  
  # Combine all features
  Xtrain_comb <- cbind(Xexpr_train_s, Xmut_train2, mm_train)
  Xtest_comb <- cbind(Xexpr_test_s, Xmut_test2, mm_test)
  colnames(Xtrain_comb) <- make.unique(colnames(Xtrain_comb))
  colnames(Xtest_comb) <- make.unique(colnames(Xtest_comb))
  
  # Combined LASSO
  cv_glm_comb <- cv.glmnet(Xtrain_comb, y_train_num, family = "binomial",
                           nfolds = N_FOLDS_INNER, type.measure = "auc", parallel = TRUE)
  fit_lasso_comb <- glmnet(Xtrain_comb, y_train_num, family = "binomial",
                           lambda = cv_glm_comb$lambda.min, standardize = FALSE)
  preds_prob_comb <- predict(fit_lasso_comb, newx = Xtest_comb, type = "response")[,1]
  roc_obj_comb <- roc(y_test_num, preds_prob_comb, quiet = TRUE)
  auc_val_comb <- as.numeric(auc(roc_obj_comb))
  
  results_summary <- bind_rows(
    results_summary,
    tibble(model = "comb_lasso", fold = fold_i, auc = auc_val_comb,
           n_selected = sum(coef(fit_lasso_comb)[-1,1] != 0))
  )
  roc_list[[paste0("comb_lasso_fold", fold_i)]] <- roc_obj_comb
  
  # Random Forest - FIXED: Extract importance directly from ranger model
  df_train_rf <- as.data.frame(Xtrain_comb)
  df_train_rf$side <- factor(ifelse(y[train_idx] == "right", "right", "left"))
  df_test_rf <- as.data.frame(Xtest_comb)
  df_test_rf$side <- factor(ifelse(y[test_idx] == "right", "right", "left"))
  
  trControl <- trainControl(method = "cv", number = N_FOLDS_INNER,
                            classProbs = TRUE, summaryFunction = twoClassSummary,
                            allowParallel = TRUE)
  
  tuneGrid <- expand.grid(
    mtry = c(25, 50, 100, max(1, floor(ncol(df_train_rf)/3))),
    splitrule = "gini",
    min.node.size = 1
  )
  tuneGrid <- tuneGrid[tuneGrid$mtry <= ncol(df_train_rf)-1, , drop = FALSE]
  
  set.seed(SEED + fold_i)
  rf_fit <- train(side ~ ., data = df_train_rf, method = "ranger",
                  trControl = trControl, tuneGrid = tuneGrid,
                  metric = "ROC", importance = "permutation")
  
  rf_preds_prob <- predict(rf_fit, newdata = df_test_rf, type = "prob")[, "right"]
  rf_roc <- roc(ifelse(y[test_idx] == "right", 1, 0), rf_preds_prob, quiet = TRUE)
  rf_auc <- as.numeric(auc(rf_roc))
  
  results_summary <- bind_rows(
    results_summary,
    tibble(model = "rf", fold = fold_i, auc = rf_auc, n_selected = NA_real_)
  )
  roc_list[[paste0("rf_fold", fold_i)]] <- rf_roc
  
  # Extract variable importance - FIXED: Extract directly from ranger model
  rf_imp <- tryCatch({
    # Get the final ranger model
    ranger_model <- rf_fit$finalModel
    
    # Check if importance exists and extract it
    if (!is.null(ranger_model$variable.importance)) {
      importance_values <- ranger_model$variable.importance
      
      # Create tibble with importance values
      tibble(
        feature = names(importance_values),
        importance = as.numeric(importance_values),
        fold = fold_i
      )
    } else {
      # Try to manually compute importance if not available
      # This is a fallback method
      cat("Using fallback method for variable importance extraction\n")
      
      # Get the trained model and predict OOB
      pred_train <- predict(ranger_model, data = df_train_rf[, -ncol(df_train_rf)])$predictions
      baseline_accuracy <- mean(pred_train == df_train_rf$side)
      
      # Compute importance by permutation
      importance_values <- sapply(colnames(df_train_rf[, -ncol(df_train_rf)]), function(feat) {
        df_permuted <- df_train_rf
        df_permuted[[feat]] <- sample(df_permuted[[feat]])
        pred_permuted <- predict(ranger_model, data = df_permuted[, -ncol(df_permuted)])$predictions
        permuted_accuracy <- mean(pred_permuted == df_permuted$side)
        baseline_accuracy - permuted_accuracy
      })
      
      tibble(
        feature = names(importance_values),
        importance = as.numeric(importance_values),
        fold = fold_i
      )
    }
  }, error = function(e) {
    cat("Error extracting RF importance for fold", fold_i, ":", e$message, "\n")
    tibble(feature = character(), importance = numeric(), fold = integer())
  })
  
  rf_importance_list[[fold_i]] <- rf_imp
}

# Aggregate results
res_df <- results_summary %>%
  group_by(model) %>%
  summarise(
    n_folds = n(),
    mean_auc = mean(auc, na.rm = TRUE),
    sd_auc = sd(auc, na.rm = TRUE),
    median_selected = median(n_selected, na.rm = TRUE)
  ) %>% arrange(desc(mean_auc))

write_tsv(res_df, file.path(RESULTS_DIR, "classification", "classification_results_summary.tsv"))

# Save feature importance
lasso_features_all <- bind_rows(lasso_feature_list)
lasso_freq <- lasso_features_all %>%
  group_by(feature) %>%
  summarise(freq = n()) %>%
  arrange(desc(freq))
write_tsv(lasso_freq, file.path(RESULTS_DIR, "classification", "lasso_feature_frequency.tsv"))

rf_imp_df <- bind_rows(rf_importance_list)
if (nrow(rf_imp_df) > 0) {
  rf_imp_summary <- rf_imp_df %>%
    group_by(feature) %>%
    summarise(
      mean_imp = mean(importance, na.rm = TRUE),
      median_imp = median(importance, na.rm = TRUE)
    ) %>%
    arrange(desc(mean_imp))
  write_tsv(rf_imp_summary, file.path(RESULTS_DIR, "classification", "rf_feature_importance.tsv"))
} else {
  cat("No RF importance data to save\n")
}

# Plot ROC curves
plot_rocs <- function(prefix, title) {
  rocs <- roc_list[grepl(prefix, names(roc_list))]
  if (length(rocs) == 0) return(NULL)
  
  png(file.path(RESULTS_DIR, "classification", paste0(prefix, "_roc_curves.png")), 
      width = 900, height = 700, res = 150)
  plot(NULL, xlim = c(1, 0), ylim = c(0, 1), xlab = "Specificity", ylab = "Sensitivity",
       main = paste(title, "- ROC curves across folds"))
  
  aucs <- c()
  for (i in seq_along(rocs)) {
    r <- rocs[[i]]
    aucs <- c(aucs, auc(r))
    lines(1 - r$specificities, r$sensitivities, col = rgb(0, 0, 1, alpha = 0.3))
  }
  
  text(0.6, 0.2, paste0("Mean AUC = ", round(mean(aucs, na.rm = TRUE), 3),
                        " (sd = ", round(sd(aucs, na.rm = TRUE), 3), ")"))
  dev.off()
}

plot_rocs("expr_lasso", "Expression LASSO")
plot_rocs("mut_lasso", "Mutation LASSO")
plot_rocs("comb_lasso", "Combined LASSO")
plot_rocs("rf", "Random Forest")

# Stability selection for LASSO (repeated subsampling)
cat("Running LASSO stability selection with", N_REPEATS_STABILITY, "subsamples (this may take time)...\n")
set.seed(SEED)
stability_counts <- foreach(rep = 1:N_REPEATS_STABILITY, .combine = rbind, .packages = c("glmnet")) %dopar% {
  # subsample 70% of samples
  samp <- sample(seq_len(nrow(X_expr)), size = floor(0.7 * nrow(X_expr)), replace = FALSE)
  Xsub <- X_expr[samp, , drop = FALSE]
  ysub <- ifelse(y[samp] == "right", 1, 0)
  
  # scale Xsub
  mns <- colMeans(Xsub, na.rm = TRUE)
  sds <- apply(Xsub, 2, sd, na.rm = TRUE)
  sds[sds == 0] <- 1
  
  Xsub_s <- sweep(Xsub, 2, mns, "-")
  Xsub_s <- sweep(Xsub_s, 2, sds, "/")
  
  cvg <- cv.glmnet(x = as.matrix(Xsub_s), y = ysub, family = "binomial", nfolds = N_FOLDS_INNER)
  lam <- cvg$lambda.min
  fit <- glmnet(x = as.matrix(Xsub_s), y = ysub, family = "binomial", lambda = lam)
  coefs <- as.matrix(coef(fit))
  sel <- rownames(coefs)[which(coefs[,1] != 0)]
  sel <- sel[sel != "(Intercept)"]
  
  tibble(rep = rep, feature = sel)
}

# Compute frequencies
stability_df <- stability_counts %>%
  group_by(feature) %>%
  summarise(freq = n()) %>%
  arrange(desc(freq))

write_tsv(stability_df, file.path(RESULTS_DIR, "classification", "lasso_stability_subsampling.tsv"))

# Save session info
writeLines(capture.output(sessionInfo()), file.path(RESULTS_DIR, "classification", "session_info.txt"))

cat("Classification modeling complete. Results saved to:", file.path(RESULTS_DIR, "classification"), "\n")