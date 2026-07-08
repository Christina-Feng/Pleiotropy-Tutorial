# A Guide for Exploring Pleiotropic Associations in Genome-Wide Association Studies using Summary Statistics

**Authors:** Christina Y. Feng, Pierre Emmanuel Sugier, Nan Zou, and Benoit Liquet

---

## Overview

This repository provides data objects and a reproducible vignette-style R workflow accompanying the tutorial *A Guide for Exploring Pleiotropic Associations in Genome-Wide Association Studies using Summary Statistics*. The materials support the breast cancer and thyroid cancer examples presented in the paper and demonstrate how to apply variant-level and gene-level pleiotropic association methods using GWAS summary statistics.

No individual-level data are provided.

---

## Repository Structure

```text
Pleiotropy-Tutorial/
├── README.md
├── Pleiotropy-Tutorial.Rproj
├── data/
│   ├── GWAS_sumstats.RData
│   └── BCAT1_sumstats.RData
├── R/
│   ├── 00_setup.R
│   ├── 01_prepare_data.R
│   ├── 02_asset.R
│   ├── 03_placo.R
│   ├── 04_gpa.R
│   ├── 05_cpbayes.R
│   ├── 06_gcpbayes.R
│   └── 07_run_all.R
├── vignettes/
│   └── pleiotropy_workflow.Rmd
└── results/
    └── .gitkeep
```

---

## Reproducible Workflow

To improve transparency and ease of use, the example analyses have been organized into a sequential, end-to-end workflow:

1. **`00_setup.R`**: Loads required libraries, verifies package installations, and configures folders.
2. **`01_prepare_data.R`**: Imports raw GWAS summary statistics and formats inputs into a unified structured file.
3. **`02_asset.R`**: Applies the subset-based meta-analysis method ASSET.
4. **`03_placo.R`**: Conducts composite null testing with PLACO.
5. **`04_gpa.R`**: Runs GPA for pleiotropy-informed and annotation-informed prioritisation.
6. **`05_cpbayes.R`**: Executes the Bayesian cross-phenotype model CPBayes.
7. **`06_gcpbayes.R`**: Applies the Bayesian group-level model GCPBayes to gene-level summary statistics.
8. **`07_run_all.R`**: Companion script that runs the entire analysis pipeline in sequence.

### Running the Workflow

To run the full workflow in one step from the root directory of the repository, use:

```r
source("R/07_run_all.R")
```

Alternatively, follow the step-by-step annotated vignette and compile it to an HTML report using:

```r
rmarkdown::render("vignettes/pleiotropy_workflow.Rmd")
```

All output results will be saved as `.rds` files inside the `results/` folder.

---

## Software & Dependencies

All analyses are conducted in R. We recommend opening the project via `Pleiotropy-Tutorial.Rproj` in RStudio to set the project working directory automatically.

The workflow uses the following R packages:

* `tidyverse`
* `here`
* `MASS`
* `devtools`
* `ASSET`
* `GPA`
* `CPBayes`
* `GCPBayes`

Most packages can be installed from CRAN or Bioconductor. The PLACO functions are sourced directly from the authors' GitHub repository using `devtools::source_url()` within the PLACO analysis script.

### Package Installation

Install the required CRAN packages using:

```r
install.packages(c(
  "tidyverse",
  "here",
  "MASS",
  "devtools",
  "CPBayes",
  "GCPBayes"
))
```

Install the Bioconductor packages using:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("ASSET")
BiocManager::install("GPA")
```

If `GPA` is not available through Bioconductor for your R version, it can alternatively be installed from GitHub using:

```r
devtools::install_github("dongjunchung/GPA")
```

The PLACO functions are sourced directly from GitHub in `R/03_placo.R`:

```r
devtools::source_url(
  "https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE"
)
```

After installing the required packages, the setup script can be used to check that the workflow runs from the project root:

```r
source("R/00_setup.R")
```

---

## References

1. Bhattacharjee S, Qi G, Chatterjee N, Wheeler W (2025). *ASSET: An R package for subset based association analysis of heterogeneous traits and subtypes*. Bioconductor. R package version 2.28.0. doi:10.18129/B9.bioc.ASSET.

2. Ray D, Chatterjee N (2020). A powerful method for pleiotropic analysis under composite null hypothesis identifies novel shared loci between type 2 diabetes and prostate cancer. *PLoS Genetics*, 16(12), e1009218. doi:10.1371/journal.pgen.1009218.

3. Chung D, Yang C, Li C, Gelernter J, Zhao H (2014). GPA: A statistical approach to prioritizing GWAS results by integrating pleiotropy information and annotation data. *PLoS Genetics*, 10, e1004787. doi:10.1371/journal.pgen.1004787.

4. Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018). A Bayesian meta-analysis method for studying cross-phenotype genetic associations using summary statistics. *PLoS Genetics*, 14, e1007139. doi:10.1371/journal.pgen.1007139.

5. Baghfalaki T, Sugier PE, Truong T, Pettitt AN, Mengersen K, Liquet B (2021). Bayesian meta-analysis models for cross-cancer genomic investigation of pleiotropic effects using group structure. *Statistics in Medicine*, 40(6), 1498--1518.
