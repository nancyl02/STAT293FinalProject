# Topic
Comparison of VAR Estimation Methods for Multivariate Time Series Forecasting

# Group
Pei-Ling Lin (plin073@ucr.edu), Nancy Lopez (nlope089@ucr.edu)

# Summary
This project investigates methods for estimating multivariate VAR(p) models for intensive longitudinal data (ILD) with both shared and individual-specific effects. We implement two approaches: (1) Multi-VAR, which jointly estimates coefficients across all subjects, and (2) LASSO VAR, applied individually to each subject. Model performance is evaluated using metrics such as estimation error, sensitivity and specificity, and RMSFE to assess forecasting accuracy and variable selection.

# File Descriptions

## 00_requirements.R
Lists and loads all R packages needed to run the project.

## 01_DataSimulation.R
Contains functions for simulating multivariate VAR(p) time series data, including both shared and individual-specific effects.

## 12_Method1_aux.R


## 21_Method2.R

## 22_Method2_aux.R

## 3_Evaluation.R

## 41_Examples.Rmd (41_Examples.pdf)

## 42_Results_Plots.Rmd (42_Results_Plots.pdf)
