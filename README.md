# ProtAge: Proteomic Aging Clock (UK Biobank, Olink)

This repository provides a **leakage-safe**, **sex-stratified** workflow to construct the Proteomic Aging Clock (**ProtAge**) from **Olink NPX proteomic data** (already **log2-scaled**) using a **stacked ensemble** (XGBoost, LightGBM, CatBoost) with an **Elastic Net meta-learner (α = 0.5)**.

ProtAge produces predicted “proteomic age” estimates that can be used as a biological aging (BA) marker in downstream analyses (e.g., ADRD risk prediction). To avoid information leakage, the clock should be trained only in the training split and then applied to the held-out test split.

---

## Overview of the ProtAge pipeline

**Input data**
- Olink **NPX (log2 scale)** proteomic matrix (columns = proteins; rows = participants)
- Required metadata: **participant ID**, **sex**, **chronological age**

**Preprocessing & modeling**
- **Missing value imputation**: **KNN imputation (k = 10)**
- **Outlier detection**: **Isolation Forest (iForest)** to identify and exclude outliers in high-dimensional proteomic space
- **Modeling**: **Stacked ensemble** combining:
  - Base learners: **XGBoost**, **LightGBM**, **CatBoost**
  - Meta-learner: **Elastic Net regression (α = 0.5)** predicting chronological age
- **Sex-wise modeling**: Separate models for **Women** and **Men**
- **Leakage avoidance**: Train/validation operations are performed within the training split; the finalized model is used to generate predictions for the held-out test split.

---

## Key results (from the manuscript)

- **Clock-development dataset performance (UK Biobank; N = 50,482):**
  - **Correlation with chronological age:** *r ≈ 0.97*
  - **RMSE:** *≈ 1.95 years*

- **Cross-validated performance (training split), stacked model by sex:**
  - **Women:** correlation = **0.9715**, RMSE = **1.9299**, MAE = **1.3084**
  - **Men:** correlation = **0.9748**, RMSE = **1.8734**, MAE = **1.2247**

Stacking consistently outperformed all single base learners (XGBoost, LightGBM, CatBoost) across sexes.

---

## Repository structure

### Scripts (run in order)
1. `ProtAge_feature_selection.R`  
   Loads NPX data, performs basic QC / feature handling (as needed), prepares the modeling matrix.

2. `ProtAge_imputation.R`  
   Performs **KNN imputation (k = 10)** for missing NPX values.

3. `ProtAge_outlier_detection.R`  
   Detects and removes outliers using **Isolation Forest (iForest)**.

4. `ProtAge_stacked_model.R`  
   Trains **sex-stratified stacked ensembles** and generates ProtAge predictions in a **70/30 train–test split** framework.

### Outputs
Each step writes a clean `.rds` artifact and a small report into:
- `results/protage_step1_feature_selection/`
- `results/protage_step2_imputation/`
- `results/protage_step3_outliers/`
- `results/protage_step4_stacked_model/`

---

## How to run

1. Prepare your proteomic dataset in R as a dataframe with:
   - ID column (e.g., `f.eid`)
   - `sex`
   - `age`
   - proteomic NPX features (numeric columns)

2. Run scripts **in order**:
   - `ProtAge_feature_selection.R`
   - `ProtAge_imputation.R`
   - `ProtAge_outlier_detection.R`
   - `ProtAge_stacked_model.R`

3. The final ProtAge estimates will be written to:
   - `results/protage_step4_stacked_model/protage_step4_all_predictions.csv`

---

## Dependencies

Core packages used:
- `dplyr`, `caret`, `glmnet`, `Metrics`, `ggplot2`
- Base learners: `xgboost`, `lightgbm`
- Optional: `catboost` *(may require manual installation)*

If CatBoost is unavailable, the pipeline can be configured to run with XGBoost + LightGBM only.

---

## Citation

If you use this code, please cite:

Mostafaei S, et al. (2025).  
**“Precision Prediction of Alzheimer's Disease and Related Dementias Using Integrative Multi-Omics Aging Clocks and Genetic Data.”**  
Manuscript.

---

## License

This project is licensed under the **MIT License**.

---

## Contact

- Dr. Shayan Mostafaei — shayan.mostafaei@ki.se  
- Dr. Sara Hägg — sara.hagg@ki.se
