# ProtAge Feature Selection & Preprocessing
library(dplyr)
library(MASS)

# Assume 'protein_data' is a NPX matrix (rows=samples, columns=proteins)
# Step 1: Remove highly correlated proteins (optional for very high-dimensional data)
cor_matrix <- cor(protein_data, use = "pairwise.complete.obs")
high_corr <- which(abs(cor_matrix) > 0.99, arr.ind = TRUE)
high_corr_pairs <- as.data.frame(high_corr)
high_corr_pairs <- high_corr_pairs[high_corr_pairs$row < high_corr_pairs$col, ]
to_remove <- unique(colnames(protein_data)[high_corr_pairs$col])
protein_data <- protein_data %>% select(-one_of(to_remove))

# Step 2: Box-Cox transform
min_vals <- sapply(protein_data, min, na.rm = TRUE)
if(any(min_vals <= 0)) {
  protein_data <- protein_data + abs(min(min_vals)) + 1
}
# Apply Box-Cox or log transformation for normalization as appropriate for the data
# Example: protein_data <- as.data.frame(lapply(protein_data, function(x) MASS::boxcox(x~1)$y))

str(protein_data)
