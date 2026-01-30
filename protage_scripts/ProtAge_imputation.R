# ============================================================
# ProtAge_imputation.R
# Missing value imputation for Olink NPX using KNN (k = 9)
# Input  : results/protage_step1_feature_selection/protage_step1_npx_ready.rds
# Output : results/protage_step2_imputation/protage_step2_imputed.rds
# Used by: ProtAge_outlier_detection.R
# ============================================================

suppressPackageStartupMessages({
  library(VIM)     # kNN()
  library(dplyr)
})

# --------------------------
# USER SETTINGS 
# --------------------------

set.seed(20250101)

IN_DIR     <- "results/protage_step1_feature_selection"
INPUT_RDS  <- file.path(IN_DIR, "protage_step1_npx_ready.rds")

OUT_DIR <- "results/protage_step2_imputation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata columns to keep unchanged (edit to match your dataset)
ID_COL  <- "f.eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# KNN settings
K_IMPUTE <- 9

# --------------------------
# LOAD DATA (AFTER STEP 1)
# --------------------------

protein_data <- readRDS(INPUT_RDS)

required_cols <- c(ID_COL, SEX_COL, AGE_COL)
missing_cols <- setdiff(required_cols, names(protein_data))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "),
       "\nEdit ID_COL/SEX_COL/AGE_COL to match your data.")
}

meta_cols <- required_cols
feature_cols <- setdiff(names(protein_data), meta_cols)

# Only impute numeric protein columns
feature_cols <- feature_cols[sapply(protein_data[, feature_cols, drop = FALSE], is.numeric)]
if (length(feature_cols) < 1) stop("No numeric protein features found to impute.")

meta_df <- protein_data[, meta_cols, drop = FALSE]
X <- protein_data[, feature_cols, drop = FALSE]

# --------------------------
# IMPUTATION REPORT (BEFORE)
# --------------------------

miss_before <- sapply(X, function(v) mean(is.na(v))) * 100
imputation_summary_before <- data.frame(
  feature = names(miss_before),
  pct_missing_before = as.numeric(miss_before),
  stringsAsFactors = FALSE
)

# --------------------------
# KNN IMPUTATION (k = 9)
# NOTE: VIM::kNN appends indicator columns by default.
# --------------------------

X_imp <- VIM::kNN(
  X,
  k = K_IMPUTE,
  imp_var = FALSE
)

# --------------------------
# IMPUTATION REPORT (AFTER)
# --------------------------

miss_after <- sapply(X_imp, function(v) mean(is.na(v))) * 100
imputation_summary_after <- data.frame(
  feature = names(miss_after),
  pct_missing_after = as.numeric(miss_after),
  stringsAsFactors = FALSE
)

imputation_report <- imputation_summary_before %>%
  left_join(imputation_summary_after, by = "feature") %>%
  mutate(
    k = K_IMPUTE,
    n_rows = nrow(X),
    n_features = ncol(X)
  )

write.csv(
  imputation_report,
  file.path(OUT_DIR, "protage_step2_imputation_report.csv"),
  row.names = FALSE
)

# --------------------------
# FINAL OUTPUT DATAFRAME
# --------------------------

protein_data_imputed <- bind_cols(meta_df, X_imp)

# Save for next step
saveRDS(protein_data_imputed, file.path(OUT_DIR, "protage_step2_imputed.rds"))

cat("\nDONE ✅ ProtAge KNN imputation completed.\n")
cat("Input :", INPUT_RDS, "\n")
cat("Output:", file.path(OUT_DIR, "protage_step2_imputed.rds"), "\n")
cat("Report:", file.path(OUT_DIR, "protage_step2_imputation_report.csv"), "\n")
cat("k =", K_IMPUTE, "| Proteins imputed:", ncol(X_imp), "| Rows:", nrow(X_imp), "\n\n")
