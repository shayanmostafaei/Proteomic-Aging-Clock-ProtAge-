# ProtAge Missing Value Imputation (KNN)
library(VIM)
set.seed(123)
id_col <- protein_data[,1]
data_to_impute <- protein_data[,-1]
protein_data_imputed <- kNN(data_to_impute, k = 10)
protein_data_imputed <- cbind(id_col, protein_data_imputed) 
