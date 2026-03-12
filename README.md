# A Guide for Exploring Pleiotropic Associations in Genome Wide Association Studies using Summary Statistics

**Authors:** Christina Y. Feng, Pierre Emmanuel Sugier, Nan Zou, and Benoit Liquet

## Overview

This repository provides data objects and reproducible R code accompanying the tutorial.

The materials support the examples presented in the paper and focus on pleiotropic analyses using GWAS summary statistics for breast cancer and thyroid cancer.

No individual level data are shared.

## Contents

- GWAS summary statistics, including regression coefficients, standard errors, Z statistics, and covariance matrices
- Breast cancer and thyroid cancer example datasets
- R scripts implementing the pleiotropy methods used in the tutorial

The methods illustrated include ASSET, PLACO, GPA, CPBayes, and GCPBayes, covering both variant level and gene level pleiotropic inference.

## Software

All analyses were conducted in R using packages including `tidyverse`, `broom`, `MASS`, `ASSET`, `PLACO`, `GPA`, `CPBayes`, `BhGLM`, and `GCPBayes`. Package versions and dependencies are documented within the scripts.

## Scope

This repository is intended for educational and methodological purposes. It demonstrates how pleiotropic association methods can be applied using GWAS summary statistics, with an emphasis on cross phenotype inference and correlated effect estimates. The results are provided for illustration and should not be interpreted as novel genetic discoveries.

## References

- Bhattacharjee S, Qi G, Chatterjee N, Wheeler W (2025). *ASSET: An R package for subset based association analysis of heterogeneous traits and subtypes*. Bioconductor. R package version 2.28.0. doi:10.18129/B9.bioc.ASSET.

- Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between type 2 diabetes and prostate cancer. *PLoS Genetics*, 16(12), e1009218. doi:10.1371/journal.pgen.1009218.

- Chung D, Yang C, Li C, Gelernter J, Zhao H (2014). GPA: A statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data. *PLoS Genetics*, 10, e1004787. doi:10.1371/journal.pgen.1004787.

- Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018). A Bayesian meta analysis method for studying cross phenotype genetic associations using summary statistics. *PLoS Genetics*, 14, e1007139. doi:10.1371/journal.pgen.1007139.

- Baghfalaki T, Sugier PE, Truong T, Pettitt AN, Mengersen K, Liquet B (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. *Statistics in Medicine*, 40(6), 1498 to 1518.
