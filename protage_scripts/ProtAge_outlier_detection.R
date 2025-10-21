# ProtAge Outlier Detection using Machine Learning (iForest) 
library(isolationForest)
set.seed(123)
# protein_data_imputed is the matrix after imputation
iso <- isolationForest$new()
iso$fit(protein_data_imputed)
scores <- iso$predict(protein_data_imputed)
quantile_cutoff <- 0.95
outlier_threshold <- quantile(scores$anomaly_score, quantile_cutoff)
outliers <- which(scores$anomaly_score > outlier_threshold)
# Exclude outliers
cleaned_protein_data <- protein_data_imputed[-outliers, ]
