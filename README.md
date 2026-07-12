# A Guide for Exploring Pleiotropic Associations in Genome-Wide Association Studies using Summary Statistics

**Authors:** Christina Y. Feng, Pierre Emmanuel Sugier, Nan Zou, and Benoit Liquet

---

## Overview

This repository provides data objects and a reproducible vignette-style R workflow accompanying the tutorial *A Guide for Exploring Pleiotropic Associations in Genome-Wide Association Studies using Summary Statistics*. The materials support the breast cancer and thyroid cancer examples presented in the paper and demonstrate how to apply variant-level and gene-level pleiotropic association methods using GWAS summary statistics.

No individual-level data are provided.

---

## Reproducible Workflow

To improve transparency and usability, the example analyses have been organised into a sequential, end-to-end workflow. 

> [!NOTE]
> You can go through this tutorial in two ways:
> - **Run the R scripts directly:** Fork or download the repository, open `Pleiotropy-Tutorial.Rproj` and run the raw `.R` scripts from the `R/` directory.
> - **Follow the step-by-step Markdown files:** Browse the detailed walkthroughs, code explanations, and method descriptions in the Markdown files linked below or located in the `vignettes/` directory.

Each step is documented in a dedicated Markdown file:

1. [**`00_setup.md`**](vignettes/00_setup.md): Loads required libraries, verifies package installations, and configures folders.
2. [**`01_prepare_data.md`**](vignettes/01_prepare_data.md): Imports raw GWAS summary statistics and formats inputs into a unified structured file.
3. [**`02_asset.md`**](vignettes/02_asset.md): Applies the subset-based meta-analysis method ASSET.
4. [**`03_placo.md`**](vignettes/03_placo.md): Conducts composite null testing with PLACO.
5. [**`04_gpa.md`**](vignettes/04_gpa.md): Runs GPA for pleiotropy-informed and annotation-informed prioritisation.
6. [**`05_cpbayes.md`**](vignettes/05_cpbayes.md): Executes the Bayesian cross-phenotype model CPBayes.
7. [**`06_gcpbayes.md`**](vignettes/06_gcpbayes.md): Applies the Bayesian group-level model GCPBayes to gene-level summary statistics.
8. [**`07_run_all.md`**](vignettes/07_run_all.md): Describes how to run the complete analysis pipeline in sequence using `R/07_run_all.R`.

### Running the Workflow

To run the full workflow in one step and compile a complete HTML report from the root directory of the repository, use:

```r
source("R/07_run_all.R")
```

This will automatically render the R Markdown report into the `results/` directory as `pleiotropy_workflow.html`. If the `rmarkdown` package is not installed, it will fall back to executing the individual R scripts sequentially.

Each step can be followed individually in its respective Markdown file listed above. All output results are saved as `.rds` files inside the `results/` folder.

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
* `rmarkdown`
* `knitr`

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
  "GCPBayes",
  "rmarkdown",
  "knitr"
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
