# ============================================================
# ProtAge_feature_selection.R
# Feature selection + basic preprocessing for Olink NPX (log2)
#
# IMPORTANT:
# - Olink NPX values are already on a log2 scale
#
# Output is used by: ProtAge_imputation.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(caret)   # findCorrelation()
})

# --------------------------
# USER SETTINGS 
# --------------------------

set.seed(20250101)

# Input: an .rds containing a data.frame named "protein_data"
# (Recommended for a clean, reproducible pipeline)
INPUT_RDS <- "data/protage_input_raw.rds"

# Output folder
OUT_DIR <- "results/protage_step1_feature_selection"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Metadata columns to keep unchanged (edit to match your dataset)
ID_COL  <- "f.eid"
SEX_COL <- "sex"
AGE_COL <- "age"

# Correlation cutoff for optional redundancy removal
COR_CUTOFF <- 0.99

# --------------------------
# LOAD DATA
# --------------------------

obj <- readRDS(INPUT_RDS)

# Accept either a data.frame directly OR a list containing protein_data
if (is.data.frame(obj)) {
  protein_data <- obj
} else if (!is.null(obj$protein_data) && is.data.frame(obj$protein_data)) {
  protein_data <- obj$protein_data
} else {
  stop("INPUT_RDS must contain a data.frame or a list with $protein_data (data.frame).")
}

# --------------------------
# SPLIT METADATA vs FEATURES
# --------------------------

meta_cols <- intersect(c(ID_COL, SEX_COL, AGE_COL), names(protein_data))
meta_df   <- if (length(meta_cols) > 0) protein_data[, meta_cols, drop = FALSE] else NULL

# Protein features = numeric columns excluding metadata
feature_cols <- setdiff(names(protein_data), meta_cols)
feature_cols <- feature_cols[sapply(protein_data[, feature_cols, drop = FALSE], is.numeric)]

if (length(feature_cols) < 2) {
  stop("Not enough numeric protein features found. Check your input and metadata column names.")
}

X <- protein_data[, feature_cols, drop = FALSE]

# --------------------------
# QC: DROP ALL-NA / CONSTANT FEATURES
# --------------------------

all_na <- names(X)[sapply(X, function(v) all(is.na(v)))]
X <- X %>% select(-any_of(all_na))

constant <- names(X)[sapply(X, function(v) {
  vv <- v[!is.na(v)]
  length(unique(vv)) < 2
})]
X <- X %>% select(-any_of(constant))

feature_cols2 <- names(X)

# Optional: Near-zero variance filter (recommended for stability)
nzv_idx <- caret::nearZeroVar(X)
drop_nzv <- if (length(nzv_idx) > 0) names(X)[nzv_idx] else character(0)
if (length(drop_nzv) > 0) X <- X %>% select(-any_of(drop_nzv))

# --------------------------
# REMOVE HIGHLY CORRELATED PROTEINS (|r| > 0.99)
# --------------------------

to_remove_corr <- character(0)
if (ncol(X) >= 2) {
  cor_matrix <- cor(X, use = "pairwise.complete.obs")
  to_remove_corr <- caret::findCorrelation(cor_matrix, cutoff = COR_CUTOFF, names = TRUE, exact = TRUE)
  if (length(to_remove_corr) > 0) {
    X <- X %>% select(-any_of(to_remove_corr))
  }
}

kept_features <- names(X)

# --------------------------
# IMPORTANT: NO TRANSFORMATION
# Olink NPX is already log2; use directly  
# --------------------------

protein_data_step1 <- if (!is.null(meta_df)) {
  bind_cols(meta_df, X)
} else {
  X
}

# --------------------------
# SAVE OUTPUTS FOR NEXT SCRIPTS
# --------------------------

saveRDS(protein_data_step1, file.path(OUT_DIR, "protage_step1_npx_ready.rds"))

artifacts <- list(
  meta_cols = meta_cols,
  removed_all_na = all_na,
  removed_constant = constant,
  removed_nzv = drop_nzv,
  removed_high_corr = to_remove_corr,
  kept_features = kept_features,
  cor_cutoff = COR_CUTOFF,
  note = "Olink NPX is already log2; no Box-Cox/log transform applied."
)

saveRDS(artifacts, file.path(OUT_DIR, "protage_step1_artifacts.rds"))

summary_tbl <- data.frame(
  step = c("start_numeric_features", "removed_all_na", "removed_constant", "removed_nzv", "removed_high_corr", "kept_features"),
  n = c(length(feature_cols), length(all_na), length(constant), length(drop_nzv), length(to_remove_corr), length(kept_features))
)
write.csv(summary_tbl, file.path(OUT_DIR, "protage_step1_summary.csv"), row.names = FALSE)

cat("\nDONE ✅ ProtAge feature selection completed (NPX used directly; no transformation).\n")
cat("Input :", INPUT_RDS, "\n")
cat("Output:", file.path(OUT_DIR, "protage_step1_npx_ready.rds"), "\n")
cat("Artifacts:", file.path(OUT_DIR, "protage_step1_artifacts.rds"), "\n")
cat("Kept proteins:", length(kept_features), "\n\n")
