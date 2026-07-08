############################################################
# Companion script for:
# A Guide for Exploring Pleiotropic Associations in GWAS
# using Summary Statistics
#
# This script reproduces the main workflow described in the
# tutorial, from loading summary statistics to applying
# pleiotropic association methods.
############################################################

rm(list = ls())

source("R/00_setup.R")
source("R/01_prepare_data.R")
source("R/02_asset.R")
source("R/03_placo.R")
source("R/04_gpa.R")
source("R/05_cpbayes.R")
source("R/06_gcpbayes.R")

cat("Workflow completed successfully.\n")
