# Topic
Comparison of VAR Estimation Methods for Multivariate Time Series Forecasting

# Group
Pei-Ling Lin (plin073@ucr.edu), Nancy Lopez (nlope089@ucr.edu)

# Summary
This project investigates methods for estimating **multivariate VAR(p) models** for **intensive longitudinal data (ILD)** with both shared and individual-specific effects. We implement two approaches: (1) **Multi-VAR**, which jointly estimates coefficients across all subjects, and (2) **LASSO VAR**, applied individually to each subject. Model performance is evaluated using metrics such as **estimation error**, **sensitivity and specificity**, and **RMSFE** to assess forecasting accuracy and variable selection.

# File Descriptions

## Requirements
**00_requirements.R** – Lists and loads all R packages needed to run the project.

## Data Simulation
**01_DataSimulation.R** – Contains functions for simulating multivariate VAR(p) time series data, including both shared and individual-specific effects.

## Methods

### 1. Multi-VAR
**12_Method1_aux.R** – Auxiliary function `multivar_B()` reshapes lists of estimated coefficients into 3D arrays `[d, d, K]`.

*Note: 12_Method1 was not needed since we used the `multivar()` package for Multi-VAR.*

### 2. LASSO VAR
**21_Method2.R** – Fits LASSO-based VAR models for individual subjects using `glmnet`.

**22_Method2_aux.R** – Helper functions:

* `gen_design_mat()`: Constructs the VAR(p) design matrix.
* `extract_B_hat()`: Extracts estimated coefficients from LASSO fits.
* `reshape_to_array()`: Converts flat matrices into 3D coefficient arrays.

## Evaluation
**3_Evaluation.R** – Contains functions to compute evaluation metrics:

* `diff_func()`: Computes squared estimation error.
* `sens_spec_func()`: Computes sensitivity and specificity.
* `calc_metrics()`: Calculates average or subject-wise metrics.
* `compute_rmsfe_p()` / `get_rmsfe_p()`: Computes RMSFE for single or multiple subjects.

## Analysis
**41_Examples.Rmd (41_Examples.pdf)** – Generates simulations for both methods and evaluates the results.

**42_Results_Plots.Rmd (42_Results_Plots.pdf)** – Generates plots summarizing the results.
