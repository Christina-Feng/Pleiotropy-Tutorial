###############
#### ASSET ####
###############

# To install ASSET package, enter:
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("ASSET")

# Load R package and prepared summary statistics
library(ASSET)
library(tidyverse)
library(here)

dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))

# Function h.traits to perform meta-analysis 

# First Element - SNPs (character)
snps <- dat$snps

# Second Element - Traits or Studies (character)
traits.lab <- dat$traits_lab

# Third Element - beta (matrix)
beta.hat <- dat$beta_hat

# Fourth Element - sigma (matrix)
sigma.hat <- dat$sigma_hat

# Fifth & Sixth Elements - Case and Control Vectors 
ncase <- dat$N11
ncntl <- dat$N00

# Seventh Element Cor Matrix (assume independent studies)
k <- length(traits.lab)
cor <- list(
  N11 = diag(ncase, nrow = k, ncol = k),
  N00 = diag(ncntl, nrow = k, ncol = k),
  N10 = matrix(0, nrow = k, ncol = k)
)

# Show two sided results using Discrete Local Maximum approach (DLM)
res <- h.traits(
  snp.vars = snps, traits.lab = traits.lab, beta.hat = beta.hat, sigma.hat = sigma.hat,
  ncase = ncase, ncntl = ncntl,
  cor = cor, cor.numr = FALSE,
  search = NULL, side = 2, meta = TRUE,
  zmax.args = NULL, meth.pval = "DLM"
)

# Show results - print first three rows of output
cat("--- Meta Summary ---\n")
print(head(h.summary(res)$Meta,3))
cat("--- Subset 1-sided ---\n")
print(head(h.summary(res)$Subset.1sided,3))
cat("--- Subset 2-sided ---\n")
print(head(h.summary(res)$Subset.2sided,3))

# Save results to the results folder
asset_summary <- list(
  meta = h.summary(res)$Meta,
  subset_1sided = h.summary(res)$Subset.1sided,
  subset_2sided = h.summary(res)$Subset.2sided
)
saveRDS(asset_summary, here::here("results", "asset_results.rds"))

cat("ASSET analysis complete. Results saved to results/asset_results.rds\n")
