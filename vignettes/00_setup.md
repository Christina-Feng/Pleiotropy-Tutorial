# Step 00: Environment Setup

This script sets up the workspace for the cross-phenotype analysis workflow used to explore potential pleiotropic signals. It checks that the required R packages are installed and creates the `results/` output directory if it does not already exist.

## R Code (`R/00_setup.R`)

```R
############################################################
# 00_setup.R
# Load required packages and define paths
############################################################

required_packages <- c(
  "tidyverse",
  "here",
  "MASS",
  "ASSET",
  "GPA",
  "CPBayes",
  "GCPBayes"
)

# Note: Bioconductor and GitHub packages should be installed manually.
# For PLACO, it is sourced directly from GitHub, so no package load is needed.
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Package not installed: ", pkg)
  }
}

data_dir <- here::here("data")
results_dir <- here::here("results")

if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
```

### Results / Console Output

This script initialises the environment. It creates the `results/` directory if it does not exist and does not print any messages unless a package is missing:
```text
# (Runs silently if all packages are installed successfully)
```

## Description & Key Tasks

1. **Verify Packages**: Checks if the required CRAN, Bioconductor, and GitHub packages are installed.
2. **Directory Initialisation**: Automatically initialises a `results/` folder in the root project directory (using the `here` package to construct pathing robustly).

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/00_setup.R")
```
