# Network-meta-analysis
Here's a concise GitHub repository description for your network meta-analysis project: Network Meta-Analysis of Breast Cancer Treatments A comprehensive R implementation for Bayesian network meta-analysis (NMA) comparing breast cancer treatment regimens. This repository contains reproducible code

# Network Meta-Analysis of Breast Cancer Treatments

A comprehensive R implementation for Bayesian network meta-analysis (NMA) comparing breast cancer treatment regimens. This repository contains reproducible code for evidence synthesis using the gemtc package.

## Key Features

- **Bayesian NMA**with random-effects consistency model
- **Treatment ranking**with SUCRA values and probability plots
- **Model diagnostics**including convergence assessment and heterogeneity testing
- **High-quality visualizations**: network plots, forest plots, ranking charts
- **Inconsistency assessment**using node-splitting and unrelated mean effects models

## Analysis Pipeline

1. Data preparation and treatment standardization
2. Network construction and visualization
3. Bayesian MCMC simulation (3 chains, 100K iterations)
4. Model validation and diagnostic checks
5. Treatment effects estimation and ranking
6. Publication-ready figure generation

## Technical Stack

- R 4.3.2 + gemtc, netmeta, ggplot2 packages
- Bayesian framework with JAGS backend
- High-resolution TIFF outputs for academic publishing

Ideal for researchers conducting systematic reviews and comparative effectiveness research in oncology.
