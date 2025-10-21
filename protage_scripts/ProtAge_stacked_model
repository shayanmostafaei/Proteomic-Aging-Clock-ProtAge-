#######################################################################################
## ProtAge: Sex-Stratified XGBoost / LightGBM / CatBoost Stacked Ensemble Model
## Author: Shayan Mostafaei  
## GitHub: https://github.com/shayanmostafaei/Proteomic-Aging-Clock-ProtAge-
#######################################################################################

# ---- 0) Setup ----
# Install and load required packages
required_pkgs <- c(
  "dplyr", "caret", "glmnet", "xgboost", "lightgbm", "matrixStats",
  "data.table", "Metrics", "purrr", "catboost", "ggplot2", "tidyr", "tibble"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    if (pkg == "catboost") {
      message("CatBoost requires manual installation: https://catboost.ai/docs/installation/r-installation.html")
    } else {
      install.packages(pkg)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

set.seed(123)

# ---- 1) Data preparation ----
# Assumption: `proteins_data` exists with columns: f.eid, sex, age, protein features
# Replace placeholder with your actual dataset
# Example: proteins_data <- read.csv("proteins_cleaned.csv")
proteins_data <- proteins_preprocessed # <- Replace with your actual data object

data_for_ml <- proteins_data %>%
  dplyr::select(-f.eid) %>% 
  dplyr::filter(!is.na(age)) %>%
  dplyr::mutate(
    sex = as.factor(sex),
    age = as.numeric(age)
  )

names(data_for_ml) <- make.names(names(data_for_ml))
data_final <- data_for_ml

# Train/Test Split 70/30
train_index <- sample(seq_len(nrow(data_final)), size = 0.7 * nrow(data_final))
train_data <- data_final[train_index, ]
test_data  <- data_final[-train_index, ]

# ---- 2) Metrics function ----
calc_metrics <- function(truth, pred) {
  valid_idx <- !is.na(pred) & !is.na(truth)
  truth <- truth[valid_idx]; pred <- pred[valid_idx]
  
  if (length(truth) < 2) return(data.frame(R_Squared = NA, Correlation = NA, RMSE = NA, MAE = NA))
  
  ss_res <- sum((truth - pred)^2)
  ss_tot <- sum((truth - mean(truth))^2)
  data.frame(
    R_Squared = 1 - ss_res/ss_tot,
    Correlation = cor(truth, pred),
    RMSE = Metrics::rmse(truth, pred),
    MAE  = Metrics::mae(truth, pred)
  )
}

# ---- 3) Model training functions ----
# XGBoost
fit_predict_xgb <- function(train_df, test_df, full_df = NULL, nrounds = 500, eta = 0.03, max_depth = 6) {
  feats <- setdiff(names(train_df), c("age", "sex"))
  dtrain <- xgboost::xgb.DMatrix(as.matrix(train_df[, feats]), label = train_df$age)
  dtest  <- xgboost::xgb.DMatrix(as.matrix(test_df[, feats]), label = test_df$age)
  param <- list(objective = "reg:squarederror", eta = eta, max_depth = max_depth, eval_metric = "rmse")
  model <- xgboost::xgb.train(param, dtrain, nrounds = nrounds, watchlist = list(train=dtrain, eval=dtest),
                              early_stopping_rounds = 25, verbose = 0)
  list(
    model = model,
    pred_test = predict(model, dtest),
    pred_full = if (!is.null(full_df)) predict(model, xgboost::xgb.DMatrix(as.matrix(full_df[, feats]))) else NULL
  )
}

# LightGBM
fit_predict_lgb <- function(train_df, test_df, full_df = NULL, nrounds = 1000, learning_rate = 0.03, num_leaves = 31) {
  feats <- setdiff(names(train_df), c("age", "sex"))
  dtrain <- lightgbm::lgb.Dataset(as.matrix(train_df[, feats]), label = train_df$age)
  val <- list(valid = lightgbm::lgb.Dataset(as.matrix(test_df[, feats]), label = test_df$age))
  params <- list(objective = "regression", metric = "rmse", learning_rate = learning_rate, num_leaves = num_leaves)
  model <- lightgbm::lgb.train(params, dtrain, nrounds = nrounds, valids = val,
                               early_stopping_rounds = 25, verbose = -1)
  list(
    model = model,
    pred_test = predict(model, as.matrix(test_df[, feats])),
    pred_full = if (!is.null(full_df)) predict(model, as.matrix(full_df[, feats])) else NULL
  )
}

# CatBoost
fit_predict_cat <- function(train_df, test_df, full_df = NULL, iterations = 1000, learning_rate = 0.03, depth = 6) {
  if (!requireNamespace("catboost", quietly = TRUE)) {
    warning("CatBoost unavailable, skipping.")
    return(list(model=NULL, pred_test=rep(NA, nrow(test_df)), pred_full=rep(NA, if(!is.null(full_df)) nrow(full_df) else 0)))
  }
  feats <- setdiff(names(train_df), c("age", "sex"))
  pool_train <- catboost::catboost.load_pool(as.matrix(train_df[, feats]), label=train_df$age)
  pool_test  <- catboost::catboost.load_pool(as.matrix(test_df[, feats]), label=test_df$age)
  params <- list(loss_function="RMSE", iterations=iterations, learning_rate=learning_rate, depth=depth,
                 od_type="Iter", od_wait=25, verbose=FALSE)
  model <- catboost::catboost.train(pool_train, pool_test, params=params)
  list(
    model = model,
    pred_test = catboost::catboost.predict(model, pool_test),
    pred_full = if (!is.null(full_df)) catboost::catboost.predict(model, catboost::catboost.load_pool(as.matrix(full_df[, feats]))) else NULL
  )
}

# ---- 4) Run models by sex ----
run_models_for_sex <- function(sex_label, train_df, test_df, full_df) {
  train_s <- dplyr::filter(train_df, sex == sex_label)
  test_s  <- dplyr::filter(test_df, sex == sex_label)
  full_s  <- dplyr::filter(full_df, sex == sex_label)
  if (nrow(train_s) < 10) return(NULL)
  
  message("Running models for sex: ", sex_label)
  xgb_res <- fit_predict_xgb(train_s, test_s, full_s)
  lgb_res <- fit_predict_lgb(train_s, test_s, full_s)
  cat_res <- fit_predict_cat(train_s, test_s, full_s)
  
  list(
    xgb_test_metrics = calc_metrics(test_s$age, xgb_res$pred_test),
    lgb_test_metrics = calc_metrics(test_s$age, lgb_res$pred_test),
    cat_test_metrics = calc_metrics(test_s$age, cat_res$pred_test),
    xgb_full_pred = xgb_res$pred_full,
    lgb_full_pred = lgb_res$pred_full,
    cat_full_pred = cat_res$pred_full,
    models = list(xgb=xgb_res$model, lgb=lgb_res$model, cat=cat_res$model),
    counts = c(n_train=nrow(train_s), n_test=nrow(test_s), n_full=nrow(full_s))
  )
}

# ---- 5) Run all models ----
sex_levels <- levels(data_final$sex)
all_sex_results <- list()
protage_preds <- data_final %>%
  dplyr::select(sex, age) %>%
  dplyr::mutate(ProtAge_XGB=NA, ProtAge_LGB=NA, ProtAge_CAT=NA, ProtAge_Stack=NA)

for (s in sex_levels) {
  res <- run_models_for_sex(s, train_data, test_data, data_final)
  all_sex_results[[s]] <- res
  if (!is.null(res)) {
    idx <- which(data_final$sex == s)
    protage_preds$ProtAge_XGB[idx] <- res$xgb_full_pred
    protage_preds$ProtAge_LGB[idx] <- res$lgb_full_pred
    protage_preds$ProtAge_CAT[idx] <- res$cat_full_pred
  }
}

# ---- 6) Stacked Elastic Net Model ----
stack_results <- list()
train_preds <- protage_preds[train_index, ]
test_preds  <- protage_preds[-train_index, ]
base_feats <- c("ProtAge_XGB","ProtAge_LGB","ProtAge_CAT")

for (s in sex_levels) {
  train_s <- train_preds %>% dplyr::filter(sex==s)
  test_s  <- test_preds %>% dplyr::filter(sex==s)
  valid_idx <- rowSums(is.na(train_s[, base_feats]))==0
  train_s <- train_s[valid_idx, ]
  if (nrow(train_s) < 5) next
  x_train <- as.matrix(train_s[, base_feats]); y_train <- train_s$age
  x_test  <- as.matrix(test_s[, base_feats])
  cv <- glmnet::cv.glmnet(x_train, y_train, alpha=0.5, nfolds=5)
  meta_mod <- glmnet::glmnet(x_train, y_train, alpha=0.5, lambda=cv$lambda.min)
  protage_preds$ProtAge_Stack[protage_preds$sex==s] <- predict(meta_mod, newx=as.matrix(protage_preds[protage_preds$sex==s, base_feats]), s=cv$lambda.min)
}

# ---- 7) Visualize final stacked predictions ----
p <- ggplot(protage_preds, aes(x=age, y=ProtAge_Stack)) +
  geom_point(size=0.5, alpha=0.8, color="#0072B2", position=position_jitter(width=0.3, height=0)) +
  geom_abline(intercept=0, slope=1, color="black", linewidth=1.2) +
  labs(title="Chronological Age vs ProtAge (Stacked Model)", x="Chronological Age", y="ProtAge (Predicted Age)") +
  theme_classic(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold", size=16),
        axis.title=element_text(face="bold"),
        axis.text=element_text(color="black")) +
  scale_x_continuous(limits=c(38,72), breaks=seq(40,70,5)) +
  scale_y_continuous(limits=c(38,72), breaks=seq(40,70,5))
print(p)

# ---- 8) Save predictions ----
write.csv(protage_preds, "protage_all_predictions_stacked.csv", row.names=FALSE)
