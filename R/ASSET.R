###############
#### ASSET ####
###############

# To install ASSET package, enter:
# if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("ASSET")

# Load R package and summary statistics
library(ASSET)
library(tidyverse)
library(here)
load("GWAS_sumstats.RData")

# Function h.traits to perform meta-analysis 

# First Element - SNPs (character)
snps <- annot$snp %>% as.character()

# Second Element - Traits or Studies (character)
traits.lab <- c("Breast","Thyroid")

# Third Element - beta (matrix)
beta.hat <- cbind(B_sumstats[2],T_sumstats[2]) %>% as.matrix() 

# Fourth Element - sigma (matrix)
sigma.hat <- cbind(B_sumstats[3],T_sumstats[3]) %>% as.matrix() 

# Fifth & Sixth Elements - Case and Control Vectors 
# N11 <- c(B_case_total, T_case_total)
# N00 <- c(B_cntl_total, T_cntl_total)
ncase <- N11
ncntl <- N00

# Seventh Element Cor Matrix (assume indept studies)
k <- 2
cor <- list(
  N11 = diag(ncase, nrow = k, ncol = k),
  N00 = diag(ncntl, nrow = k, ncol = k),
  N10 = matrix(0, nrow = k, ncol = k)
)

# Show two sided results using Discrete Local Maximum approach (DLM)
res <- h.traits(
  snps, traits.lab, beta.hat, sigma.hat,
  ncase, ncntl,
  cor = cor, cor.numr = FALSE,
  search = NULL, side = 2, meta = TRUE,
  zmax.args = NULL, meth.pval = "DLM"
)

# Show results - print first three rows of output
head(h.summary(res)$Meta,3)
head(h.summary(res)$Subset.1sided,3)
head(h.summary(res)$Subset.2sided,3)
