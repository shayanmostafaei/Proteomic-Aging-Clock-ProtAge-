# ProtAge: Proteomic Aging Clock

This repository contains the code and workflow for constructing the Proteomic Aging Clock ("ProtAge") using Olink NPX proteomic data from the UK Biobank.

## Overview

- **Normalization**: Box-Cox transformation to address skewed distributions.
- **Missing Value Imputation**: k-Nearest Neighbors (k=10).
- **Outlier Detection**: Machine learning methods (Isolation Forest). 
- **Modeling**: Stacked ensemble (XGBoost, Random Forest, Neural Network) with Elastic Net Regression as a meta-learner.
- **Sex-wise Modeling**: Separate models for Men and Women, addressing sex-specific proteomic aging rates.
- **Subgroup Modeling**: Separate models for age (Under 50, 50-59, 60 and older) and sex groups.

  ## Files

- `ProtAge_feature_selection.R`: Feature selection and preprocessing.
- `ProtAge_imputation.R`: KNN imputation for missing values.
- `ProtAge_outlier_detection.R`: Outlier detection and exclusion.
- `ProtAge_stacked_model.R`: Main model training and evaluation.
-  ProtAge_sexwise_models.R`:Dedicated script for training and evaluating ProtAge separately for Women and Men
- `ProtAge_groupwise_models.R`: Age/sex subgroup modeling.
- `README.md`: This file.
- `LICENSE`: License

## How to Use

1. Prepare your Olink NPX proteomic data in R.
2. Run the scripts in order as listed above.
3. Follow the instructions and comments in each script for detailed steps and configuration.

## Reference

If you use this code, please cite:

> Mostafaei S, et al. (2025) "Precision Prediction of Alzheimer's Disease and Related Dementias Using Integrative Multi-Omics Aging Clocks and Genetic Data" [Manuscript].  

## License

This project is licensed under the MIT License

## Contact

For questions or contributions, please contact: • Dr. Shayan Mostafaei (shayan.mostafaei@ki.se) • Dr. Sara Hägg (sara.hagg@ki.se)
