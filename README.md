# Project Overview

This repo holds R scripts for preparing data, training models (IWRF, XGBoost, horseshoe logistic), calibrating with a two-stage GP (BNE), and making plots.

## Files

- **Data_processing.R**: convert STATLOG dataset from .dat to .csv.
- **IWRF.R**: select features with Inf-FS, tune & run IWRF. 
- **IWRF_XGB.R**: blend IWRF and XGBoost outputs.
- **horseshoe.R**: fit logistic-horseshoe in Stan.
- **BNE.R**: take IWRF OOFs + PCs, fit two-stage GP in Stan.
- **Plots.R**: read saved results and generate plots in the appendices.