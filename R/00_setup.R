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
  "GCPBayes",
  "rmarkdown",
  "knitr"
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
