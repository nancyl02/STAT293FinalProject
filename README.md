# Topic
Comparison of VAR Estimation Methods for Multivariate Time Series Forecasting

# Group
Pei-Ling Lin (plin073@ucr.edu), Nancy Lopez (nlope089@ucr.edu)

# Summary
This project investigates methods for estimating **multivariate VAR(p) models** for **intensive longitudinal data (ILD)** with both shared and individual-specific effects. We implement two approaches: (1) **Multi-VAR**, which jointly estimates coefficients across all subjects, and (2) **LASSO VAR**, applied individually to each subject. Model performance is evaluated using metrics such as **estimation error**, **sensitivity and specificity**, and **RMSFE** to assess forecasting accuracy and variable selection.

# File Descriptions

## RCode 

### Requirements
**00_requirements.R** – Lists and loads all R packages needed to run the project.

### Data Simulation
**01_DataSimulation.R** – Contains functions for simulating multivariate VAR(p) time series data, including both shared and individual-specific effects.

### Methods

#### 1. Multi-VAR
**12_Method1_aux.R** – Auxiliary function `multivar_B()` reshapes lists of estimated coefficients into 3D arrays `[d, d, K]`.

*Note: 12_Method1 was not needed since we used the `multivar()` package for Multi-VAR.*

#### 2. LASSO VAR
**21_Method2.R** – Fits LASSO-based VAR models for individual subjects using `glmnet`.

**22_Method2_aux.R** – Helper functions:

* `gen_design_mat()`: Constructs the VAR(p) design matrix.
* `extract_B_hat()`: Extracts estimated coefficients from LASSO fits.
* `reshape_to_array()`: Converts flat matrices into 3D coefficient arrays.

### Evaluation
**3_Evaluation.R** – Contains functions to compute evaluation metrics:

* `diff_func()`: Computes squared estimation error.
* `sens_spec_func()`: Computes sensitivity and specificity.
* `calc_metrics()`: Calculates average or subject-wise metrics.
* `compute_rmsfe_p()` / `get_rmsfe_p()`: Computes RMSFE for single or multiple subjects.

### Analysis
**41_Examples.Rmd (41_Examples.pdf)** – Generates simulations for both methods and evaluates the results.

**42_Results Folder** - Contains all simulation outputs for the project, including `.pdf`, `.Rmd`, and `.RData` files, provided for reference.

**42_Results_Plots.Rmd (42_Results_Plots.pdf)** – Generates plots summarizing the results.

## Report

**FinalReport.pdf** - Final compiled report summarizing the background information, methodology, implementation and simulations, results, and conclusions.

**FinalReport.zip** - Source files used to generate `FinalReport.pdf` (e.g., LaTeX/RMarkdown, figures, and supporting code).

## Slides

**Slides.pdf** - Final presentation slides summarizing the key background, methodology, simulations, and results.

**Slides.zip** - Beamer (LaTeX) source files used to generate `Slides.pdf`.

# Reproducing the Results

To reproduce the results in this project, just follow these steps:

1. **Set up your environment**  
   - Open RStudio and set your working directory to the project folder.  
   - Run `00_requirements.R` to install and load all necessary packages.

2. **Run example analyses**  
   - Open `41_Examples.Rmd` and adjust parameters directly in *Step 2*.  
   - This will automatically simulate data, fit both Multi-VAR and LASSO VAR models, and compute all evaluation metrics.

3. **Generate plots**  
   - Open `42_Results_Plots.Rmd` and update the results dataset as needed.  
   - This produces figures summarizing model performance across different settings.



