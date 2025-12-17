#######################################################################################
## ProtAge: Sex-Stratified XGBoost / LightGBM / CatBoost Stacked Ensemble Model
##
## Pipeline assumption:
## - Step 1: Feature selection / QC 
## - Step 2: KNN imputation (k=10)
## - Step 3: Outlier detection (Isolation Forest) and removal
##
## Input  : results/protage_step3_outliers/protage_step3_cleaned_no_outliers.rds
## Output : results/protage_step4_stacked_model/
#######################################################################################

# ---- 0) Setup ----
required_pkgs <- c("dplyr","caret","xgboost","lightgbm","glmnet","Metrics","ggplot2")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# CatBoost is optional (manual install on some systems)
HAS_CATBOOST <- requireNamespace("catboost", quietly = TRUE)
if (HAS_CATBOOST) suppressPackageStartupMessages(library(catboost))
if (!HAS_CATBOOST) message("NOTE: 'catboost' not installed. CatBoost will be skipped.\n",
                           "Install instructions: https://catboost.ai/docs/installation/r-installation.html")

set.seed(20250101)

# ---- USER SETTINGS ----
IN_DIR    <- "results/protage_step3_outliers"
INPUT_RDS <- file.path(IN_DIR, "protage_step3_cleaned_no_outliers.rds")

OUT_DIR <- "results/protage_step4_stacked_model"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Column names (edit to match your dataset)
ID_COL  <- "f.eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# Holdout + CV settings
HOLDOUT_FRAC <- 0.30     # 70/30 split within each sex
N_FOLDS      <- 5        # OOF folds within TRAIN only (for meta-learner)

# Base model hyperparameters 
XGB_PARAMS <- list(
  objective = "reg:squarederror",
  eta = 0.03,
  max_depth = 6,
  subsample = 0.9,
  colsample_bytree = 0.9,
  eval_metric = "rmse"
)
XGB_NROUNDS <- 800

LGB_PARAMS <- list(
  objective = "regression",
  metric = "rmse",
  learning_rate = 0.03,
  num_leaves = 31,
  feature_fraction = 0.9,
  bagging_fraction = 0.9,
  bagging_freq = 1
)
LGB_NROUNDS <- 3000
LGB_EARLY_STOP <- 50

CAT_PARAMS <- list(
  loss_function = "RMSE",
  iterations = 2000,
  learning_rate = 0.03,
  depth = 6,
  od_type = "Iter",
  od_wait = 50,
  verbose = FALSE
)

# Elastic Net meta-learner
META_ALPHA <- 0.5

# ---- 1) Load data (already cleaned/imputed/no-outliers) ----
data_in <- readRDS(INPUT_RDS)

required_cols <- c(ID_COL, SEX_COL, AGE_COL)
missing_cols <- setdiff(required_cols, names(data_in))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "),
       "\nEdit ID_COL/SEX_COL/AGE_COL to match your data.")
}

data_ml <- data_in %>%
  dplyr::filter(!is.na(.data[[AGE_COL]]), !is.na(.data[[SEX_COL]])) %>%
  dplyr::mutate(
    !!SEX_COL := as.factor(.data[[SEX_COL]]),
    !!AGE_COL := as.numeric(.data[[AGE_COL]]),
    !!ID_COL  := as.character(.data[[ID_COL]])
  )

feature_cols <- setdiff(names(data_ml), required_cols)
feature_cols <- feature_cols[sapply(data_ml[, feature_cols, drop = FALSE], is.numeric)]
if (length(feature_cols) < 2) stop("Need at least 2 numeric protein features for modeling.")

# Safe names (avoid special-char issues in ML libraries)
safe_names <- make.names(feature_cols, unique = TRUE)
names_map <- setNames(safe_names, feature_cols)
names(data_ml)[match(feature_cols, names(data_ml))] <- safe_names
feature_cols <- safe_names

# ---- 2) Metrics function ----
calc_metrics <- function(truth, pred){
  valid_idx <- !is.na(pred) & !is.na(truth)
  truth <- truth[valid_idx]; pred <- pred[valid_idx]
  if(length(truth) < 2) return(data.frame(R_Squared=NA, Correlation=NA, RMSE=NA, MAE=NA))
  ss_res <- sum((truth - pred)^2)
  ss_tot <- sum((truth - mean(truth))^2)
  data.frame(
    R_Squared   = 1 - ss_res/ss_tot,
    Correlation = cor(truth, pred),
    RMSE        = Metrics::rmse(truth, pred),
    MAE         = Metrics::mae(truth, pred)
  )
}

# ---- 3) Base model helpers ----
fit_xgb <- function(train_df, valid_df, feats, label_col){
  dtrain <- xgboost::xgb.DMatrix(as.matrix(train_df[, feats, drop = FALSE]), label = train_df[[label_col]])
  dvalid <- xgboost::xgb.DMatrix(as.matrix(valid_df[, feats, drop = FALSE]), label = valid_df[[label_col]])
  xgboost::xgb.train(
    params = XGB_PARAMS,
    data = dtrain,
    nrounds = XGB_NROUNDS,
    watchlist = list(valid = dvalid),
    verbose = 0
  )
}
pred_xgb <- function(model, df, feats){
  dmat <- xgboost::xgb.DMatrix(as.matrix(df[, feats, drop = FALSE]))
  as.numeric(predict(model, dmat))
}

fit_lgb <- function(train_df, valid_df, feats, label_col){
  dtrain <- lightgbm::lgb.Dataset(as.matrix(train_df[, feats, drop = FALSE]), label = train_df[[label_col]])
  dvalid <- lightgbm::lgb.Dataset(as.matrix(valid_df[, feats, drop = FALSE]), label = valid_df[[label_col]])
  lightgbm::lgb.train(
    params = LGB_PARAMS,
    data = dtrain,
    nrounds = LGB_NROUNDS,
    valids = list(valid = dvalid),
    early_stopping_rounds = LGB_EARLY_STOP,
    verbose = -1
  )
}
pred_lgb <- function(model, df, feats){
  as.numeric(predict(model, as.matrix(df[, feats, drop = FALSE])))
}

fit_cat <- function(train_df, valid_df, feats, label_col){
  if (!HAS_CATBOOST) return(NULL)
  pool_train <- catboost::catboost.load_pool(as.matrix(train_df[, feats, drop = FALSE]), label = train_df[[label_col]])
  pool_valid <- catboost::catboost.load_pool(as.matrix(valid_df[, feats, drop = FALSE]), label = valid_df[[label_col]])
  catboost::catboost.train(pool_train, pool_valid, params = CAT_PARAMS)
}
pred_cat <- function(model, df, feats){
  if (is.null(model)) return(rep(NA_real_, nrow(df)))
  pool <- catboost::catboost.load_pool(as.matrix(df[, feats, drop = FALSE]))
  as.numeric(catboost::catboost.predict(model, pool))
}

# ---- 4) Sex-stratified 70/30 split + OOF stacking within TRAIN ----
sex_levels <- levels(data_ml[[SEX_COL]])

pred_all <- data_ml %>%
  dplyr::select(all_of(c(ID_COL, SEX_COL, AGE_COL))) %>%
  dplyr::mutate(
    ProtAge_XGB      = NA_real_,
    ProtAge_LGB      = NA_real_,
    ProtAge_CatBoost = NA_real_,
    ProtAge_Stack    = NA_real_,  # OOF in train + final in test
    Split            = NA_character_
  )

metrics_list <- list()
models_by_sex <- list()

for (s in sex_levels) {

  df_s <- data_ml %>% dplyr::filter(.data[[SEX_COL]] == s)
  if (nrow(df_s) < 50) {
    message("Skipping sex='", s, "' due to small sample size (n=", nrow(df_s), ").")
    next
  }

  # Holdout split within sex
  idx_train <- caret::createDataPartition(df_s[[AGE_COL]], p = 1 - HOLDOUT_FRAC, list = FALSE)
  train_s <- df_s[idx_train, , drop = FALSE]
  test_s  <- df_s[-idx_train, , drop = FALSE]

  pred_all$Split[pred_all[[SEX_COL]] == s & pred_all[[ID_COL]] %in% train_s[[ID_COL]]] <- "train"
  pred_all$Split[pred_all[[SEX_COL]] == s & pred_all[[ID_COL]] %in% test_s[[ID_COL]]]  <- "test"

  # Folds for OOF predictions in TRAIN only
  folds <- caret::createFolds(train_s[[AGE_COL]], k = N_FOLDS, list = TRUE, returnTrain = FALSE)

  oof_xgb <- rep(NA_real_, nrow(train_s))
  oof_lgb <- rep(NA_real_, nrow(train_s))
  oof_cat <- rep(NA_real_, nrow(train_s))

  # ---- OOF loop ----
  for (k in seq_along(folds)) {
    val_idx <- folds[[k]]
    tr_idx  <- setdiff(seq_len(nrow(train_s)), val_idx)

    tr <- train_s[tr_idx, , drop = FALSE]
    va <- train_s[val_idx, , drop = FALSE]

    mod_xgb <- fit_xgb(tr, va, feature_cols, AGE_COL)
    oof_xgb[val_idx] <- pred_xgb(mod_xgb, va, feature_cols)

    mod_lgb <- fit_lgb(tr, va, feature_cols, AGE_COL)
    oof_lgb[val_idx] <- pred_lgb(mod_lgb, va, feature_cols)

    if (HAS_CATBOOST) {
      mod_cat <- fit_cat(tr, va, feature_cols, AGE_COL)
      oof_cat[val_idx] <- pred_cat(mod_cat, va, feature_cols)
    }
  }

  # Meta-learner training data (OOF)
  meta_train <- data.frame(
    XGB = oof_xgb,
    LGB = oof_lgb,
    Cat = if (HAS_CATBOOST) oof_cat else NULL,
    age = train_s[[AGE_COL]]
  )
  meta_train <- meta_train[complete.cases(meta_train), , drop = FALSE]

  if (nrow(meta_train) < 50) {
    message("Skipping stacking for sex='", s, "' due to insufficient complete OOF rows.")
    next
  }

  x_meta <- as.matrix(meta_train[, setdiff(names(meta_train), "age"), drop = FALSE])
  y_meta <- meta_train$age

  cv <- glmnet::cv.glmnet(x_meta, y_meta, alpha = META_ALPHA, nfolds = N_FOLDS)
  meta_mod <- glmnet::glmnet(x_meta, y_meta, alpha = META_ALPHA, lambda = cv$lambda.min)

  # Fit final base models on full TRAIN, then predict TRAIN+TEST
  mod_xgb_full <- fit_xgb(train_s, test_s, feature_cols, AGE_COL)
  mod_lgb_full <- fit_lgb(train_s, test_s, feature_cols, AGE_COL)
  mod_cat_full <- if (HAS_CATBOOST) fit_cat(train_s, test_s, feature_cols, AGE_COL) else NULL

  base_train <- data.frame(
    XGB = pred_xgb(mod_xgb_full, train_s, feature_cols),
    LGB = pred_lgb(mod_lgb_full, train_s, feature_cols),
    Cat = if (HAS_CATBOOST) pred_cat(mod_cat_full, train_s, feature_cols) else NULL
  )
  base_test <- data.frame(
    XGB = pred_xgb(mod_xgb_full, test_s, feature_cols),
    LGB = pred_lgb(mod_lgb_full, test_s, feature_cols),
    Cat = if (HAS_CATBOOST) pred_cat(mod_cat_full, test_s, feature_cols) else NULL
  )

  # TRAIN stacked: use OOF base predictions (leakage-safe)
  stack_train <- rep(NA_real_, nrow(train_s))
  ok_rows <- which(complete.cases(data.frame(XGB=oof_xgb, LGB=oof_lgb, Cat=if(HAS_CATBOOST) oof_cat else NULL)))
  stack_train[ok_rows] <- as.numeric(predict(
    meta_mod,
    newx = as.matrix(data.frame(
      XGB=oof_xgb[ok_rows],
      LGB=oof_lgb[ok_rows],
      Cat=if(HAS_CATBOOST) oof_cat[ok_rows] else NULL
    )),
    s = cv$lambda.min
  ))

  # TEST stacked: use final base predictions + meta learner
  stack_test <- as.numeric(predict(meta_mod, newx = as.matrix(base_test), s = cv$lambda.min))

  # Write predictions back (map by ID)
  idx_train_global <- which(pred_all[[SEX_COL]] == s & pred_all[[ID_COL]] %in% train_s[[ID_COL]])
  idx_test_global  <- which(pred_all[[SEX_COL]] == s & pred_all[[ID_COL]] %in% test_s[[ID_COL]])

  # Base preds
  pred_all$ProtAge_XGB[idx_train_global] <- base_train$XGB
  pred_all$ProtAge_LGB[idx_train_global] <- base_train$LGB
  if (HAS_CATBOOST) pred_all$ProtAge_CatBoost[idx_train_global] <- base_train$Cat

  pred_all$ProtAge_XGB[idx_test_global] <- base_test$XGB
  pred_all$ProtAge_LGB[idx_test_global] <- base_test$LGB
  if (HAS_CATBOOST) pred_all$ProtAge_CatBoost[idx_test_global] <- base_test$Cat

  # Stacked preds
  train_stack_map <- data.frame(id=train_s[[ID_COL]], stack=stack_train, stringsAsFactors = FALSE)
  pred_all$ProtAge_Stack[idx_train_global] <- train_stack_map$stack[match(pred_all[[ID_COL]][idx_train_global], train_stack_map$id)]
  pred_all$ProtAge_Stack[idx_test_global]  <- stack_test

  # Evaluate on TEST
  metrics_list[[length(metrics_list)+1]] <- calc_metrics(test_s[[AGE_COL]], base_test$XGB) %>% mutate(Sex=s, Model="XGBoost")
  metrics_list[[length(metrics_list)+1]] <- calc_metrics(test_s[[AGE_COL]], base_test$LGB) %>% mutate(Sex=s, Model="LightGBM")
  if (HAS_CATBOOST) metrics_list[[length(metrics_list)+1]] <- calc_metrics(test_s[[AGE_COL]], base_test$Cat) %>% mutate(Sex=s, Model="CatBoost")
  metrics_list[[length(metrics_list)+1]] <- calc_metrics(test_s[[AGE_COL]], stack_test) %>% mutate(Sex=s, Model="Stacked_ElasticNet(alpha=0.5)")

  # Save models for deployment
  models_by_sex[[as.character(s)]] <- list(
    sex_level = s,
    feature_cols = feature_cols,
    name_map = names_map,
    holdout_frac = HOLDOUT_FRAC,
    n_folds = N_FOLDS,
    base_models = list(xgb = mod_xgb_full, lgb = mod_lgb_full, cat = mod_cat_full),
    meta_model = meta_mod,
    meta_lambda = cv$lambda.min,
    meta_alpha = META_ALPHA
  )

  message("Finished sex='", s, "': train n=", nrow(train_s), ", test n=", nrow(test_s))
}

# ---- 5) Save metrics ----
metrics_tbl <- dplyr::bind_rows(metrics_list)
write.csv(metrics_tbl, file.path(OUT_DIR, "protage_step4_test_metrics_by_sex.csv"), row.names = FALSE)

test_rows <- pred_all$Split == "test"
overall_metrics <- calc_metrics(pred_all[[AGE_COL]][test_rows], pred_all$ProtAge_Stack[test_rows]) %>%
  mutate(Sex="Overall", Model="Stacked_ElasticNet(alpha=0.5)")
write.csv(overall_metrics, file.path(OUT_DIR, "protage_step4_test_metrics_overall.csv"), row.names = FALSE)

# ---- 6) Plot (TEST only; stacked) ----
p <- ggplot2::ggplot(pred_all[test_rows, ], ggplot2::aes(x=.data[[AGE_COL]], y=ProtAge_Stack)) +
  ggplot2::geom_point(size=0.6, alpha=0.7) +
  ggplot2::geom_abline(intercept=0, slope=1, linewidth=1) +
  ggplot2::labs(
    title = paste0("Chronological Age vs ProtAge (Stacked, TEST only; 70/30 holdout)\n",
                   "Correlation=", round(overall_metrics$Correlation, 3),
                   " | RMSE=", round(overall_metrics$RMSE, 3),
                   " | MAE=", round(overall_metrics$MAE, 3)),
    x = "Chronological Age",
    y = "ProtAge (Predicted Age)"
  ) +
  ggplot2::theme_classic(base_size = 14)

ggplot2::ggsave(
  filename = file.path(OUT_DIR, "protage_step4_age_vs_protage_test.png"),
  plot = p, width = 7, height = 5, dpi = 300
)

# ---- 7) Save predictions + models ----
write.csv(pred_all, file.path(OUT_DIR, "protage_step4_all_predictions.csv"), row.names = FALSE)
saveRDS(models_by_sex, file.path(OUT_DIR, "protage_step4_models_by_sex.rds"))

cat("\nDONE ✅ ProtAge stacked modeling complete.\n")
cat("Input :", INPUT_RDS, "\n")
cat("Predictions:", file.path(OUT_DIR, "protage_step4_all_predictions.csv"), "\n")
cat("Models     :", file.path(OUT_DIR, "protage_step4_models_by_sex.rds"), "\n")
cat("Metrics    :", file.path(OUT_DIR, "protage_step4_test_metrics_by_sex.csv"), "\n\n")
