# ============================================================
# ProtAge_outlier_detection.R
# Outlier detection using Isolation Forest (iForest)
# Input  : results/protage_step2_imputation/protage_step2_imputed.rds
# Output : results/protage_step3_outliers/protage_step3_cleaned_no_outliers.rds
# Used by: ProtAge_stacked_model.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# --------------------------
# USER SETTINGS 
# --------------------------

set.seed(20250101)

IN_DIR     <- "results/protage_step2_imputation"
INPUT_RDS  <- file.path(IN_DIR, "protage_step2_imputed.rds")

OUT_DIR <- "results/protage_step3_outliers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata columns to keep unchanged (edit to match your dataset)
ID_COL  <- "f.eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# Isolation Forest settings
# quantile cutoff for outliers: 0.95 means flag top 5% highest anomaly scores
QUANTILE_CUTOFF <- 0.95

# Model settings (reasonable defaults)
N_TREES   <- 500
SAMPLE_SIZE <- 256  # typical default; will be capped at nrow(X) automatically

# Prefer 'isotree' (widely used + stable). 
HAS_ISOTREE <- requireNamespace("isotree", quietly = TRUE)
HAS_ISOFOREST_PKG <- requireNamespace("isolationForest", quietly = TRUE)

if (!HAS_ISOTREE && !HAS_ISOFOREST_PKG) {
  stop("Please install either 'isotree' (recommended) or 'isolationForest' package.")
}

if (HAS_ISOTREE) suppressPackageStartupMessages(library(isotree))

# --------------------------
# LOAD DATA (AFTER IMPUTATION)
# --------------------------

protein_data_imputed <- readRDS(INPUT_RDS)

required_cols <- c(ID_COL, SEX_COL, AGE_COL)
missing_cols <- setdiff(required_cols, names(protein_data_imputed))
if (length(missing_cols) > 0) {
  stop("Missing required columns in input: ", paste(missing_cols, collapse = ", "),
       "\nEdit ID_COL/SEX_COL/AGE_COL to match your data.")
}

meta_cols <- required_cols
feature_cols <- setdiff(names(protein_data_imputed), meta_cols)

# Only numeric protein columns go into iForest
feature_cols <- feature_cols[sapply(protein_data_imputed[, feature_cols, drop = FALSE], is.numeric)]
if (length(feature_cols) < 2) stop("Need at least 2 numeric protein features for iForest outlier detection.")

X <- protein_data_imputed[, feature_cols, drop = FALSE]

# Safety: ensure no NA remain (imputation should have removed all)
if (anyNA(X)) stop("NA values detected after imputation. Fix imputation before outlier detection.")

# --------------------------
# FIT ISOLATION FOREST + SCORE
# --------------------------

if (HAS_ISOTREE) {
  # isotree returns anomaly scores: higher = more anomalous
  sample_size_used <- min(SAMPLE_SIZE, nrow(X))
  iso_model <- isotree::isolation.forest(
    data = as.matrix(X),
    ntrees = N_TREES,
    sample_size = sample_size_used
  )
  anomaly_score <- isotree::predict.isolation_forest(iso_model, as.matrix(X), type = "score")
  model_used <- "isotree::isolation.forest"
} else {
  # fallback to isolationForest package (R6 style, like your original)
  suppressPackageStartupMessages(library(isolationForest))
  iso <- isolationForest::isolationForest$new()
  iso$fit(X)
  scores <- iso$predict(X)
  anomaly_score <- scores$anomaly_score
  iso_model <- iso
  model_used <- "isolationForest::isolationForest$new()"
}

# --------------------------
# OUTLIER THRESHOLD + FLAG
# --------------------------

outlier_threshold <- as.numeric(stats::quantile(anomaly_score, probs = QUANTILE_CUTOFF, na.rm = TRUE))
is_outlier <- anomaly_score > outlier_threshold
outliers_idx <- which(is_outlier)

# Exclude outliers
cleaned_protein_data <- protein_data_imputed[!is_outlier, , drop = FALSE]

# --------------------------
# SAVE OUTPUTS + REPORT
# --------------------------

outlier_report <- protein_data_imputed %>%
  select(all_of(meta_cols)) %>%
  mutate(
    anomaly_score = as.numeric(anomaly_score),
    quantile_cutoff = QUANTILE_CUTOFF,
    threshold = outlier_threshold,
    is_outlier = is_outlier,
    model = model_used,
    n_trees = N_TREES,
    sample_size = if (HAS_ISOTREE) min(SAMPLE_SIZE, nrow(X)) else NA
  )

write.csv(outlier_report, file.path(OUT_DIR, "protage_step3_outlier_report.csv"), row.names = FALSE)

saveRDS(cleaned_protein_data, file.path(OUT_DIR, "protage_step3_cleaned_no_outliers.rds"))

# Save artifacts for reproducibility
artifacts <- list(
  model_used = model_used,
  quantile_cutoff = QUANTILE_CUTOFF,
  threshold = outlier_threshold,
  n_trees = N_TREES,
  sample_size = if (HAS_ISOTREE) min(SAMPLE_SIZE, nrow(X)) else NA,
  has_isotree = HAS_ISOTREE
)
saveRDS(artifacts, file.path(OUT_DIR, "protage_step3_outlier_artifacts.rds"))

cat("\nDONE ✅ ProtAge outlier detection (Isolation Forest) completed.\n")
cat("Input :", INPUT_RDS, "\n")
cat("Model :", model_used, "\n")
cat("Quantile cutoff :", QUANTILE_CUTOFF, " (flagged top ", round((1-QUANTILE_CUTOFF)*100, 1), "%)\n", sep = "")
cat("Outliers detected:", length(outliers_idx), "of", nrow(protein_data_imputed), "\n")
cat("Output:", file.path(OUT_DIR, "protage_step3_cleaned_no_outliers.rds"), "\n")
cat("Report:", file.path(OUT_DIR, "protage_step3_outlier_report.csv"), "\n\n")
