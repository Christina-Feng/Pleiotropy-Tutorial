# Step 05: CPBayes Analysis

This script applies **CPBayes** (Bayesian cross-phenotype association analysis) to study multi-trait genetic associations using summary statistics. It models independent traits, estimates posterior probabilities of association (PPA), and identifies which specific traits contribute to the pleiotropic signal.

## Step-by-Step Analysis & Results

### 1. Load CPBayes & Structured GWAS Statistics
Load packages, prepared data, and target SNP stats:

```R
library(CPBayes)
library(here)
library(tidyverse)

# Load summary statistics
dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))

# Extract illustrative variant rs12914272
rs12914272 <- rbind(dat$B_sumstats %>% filter(SNPs == "rs12914272"), 
                    dat$T_sumstats %>% filter(SNPs == "rs12914272"))

Beta_Hat <- rs12914272$estimate
SE <- rs12914272$std.error
Traits <- c("Breast","Thyroid")
SNP <- "rs12914272"
```

### 2. Compute Analytical locFDR and log10_BF
Perform single-variant pre-filtering calculation:

```R
result <- analytic_locFDR_BF_uncor(Beta_Hat, SE)

cat("--- Analytic locFDR and log10_BF ---\n")
print(result)
```

**Output:**
```text
--- Analytic locFDR and log10_BF ---
$locFDR
[1] 0.009236311

$log10_BF
[1] 1.55335
```

### 3. Fit CPBayes MCMC Model on Single Variant
Estimate posterior metrics via MCMC simulations for `rs12914272` (5,500 iterations, 500 burn-in):

```R
# Fit MCMC model
result_CP <- cpbayes_uncor(Beta_Hat, SE, 
                           Phenotypes = Traits, Variant = SNP, 
                           MCMCiter = 5500, Burnin = 500)

# Summary of posterior probability of associations (PPAj)
PleioSumm <- post_summaries(result_CP, level = 0.05) 

cat("--- Posterior Summaries: Important Traits ---\n")
print(PleioSumm$important_traits)
cat("--- Posterior Summaries: Beta ---\n")
print(PleioSumm$poste_summary_beta)
cat("--- Posterior Summaries: OR ---\n")
print(PleioSumm$poste_summary_OR)
```

**Output:**
```text
--- Posterior Summaries: Important Traits ---
   traits   PPAj direction
1  Breast 0.9970  negative
2 Thyroid 0.3658  positive
--- Posterior Summaries: Beta ---
   traits  poste_mean poste_median   poste_se         lCl        uCl
1  Breast -0.23334679 -0.233277756 0.06167295 -0.35045995 -0.1124438
2 Thyroid  0.03911667  0.009368279 0.06087601 -0.01710717  0.1905788
--- Posterior Summaries: OR ---
   traits poste_mean poste_median       lCl       uCl
1  Breast  0.7933895    0.7919336 0.7043640 0.8936476
2 Thyroid  1.0418724    1.0094123 0.9830383 1.2099498
```

### 4. Run CPBayes MCMC SNP-by-SNP for Filtered Targets
Fast-filter variants (locFDR $< 0.05$ & log10(BF) $> 1$) and perform full MCMC runs on target variants:

```R
# Extract matrices
all_snps <- dat$snps
beta_mat <- dat$beta_hat
se_mat <- dat$sigma_hat

# Compute fast analytical locFDR and log10_BF for all 3,707 SNPs
bfs <- sapply(1:length(all_snps), function(i) {
  res <- tryCatch(analytic_locFDR_BF_uncor(beta_mat[i,], se_mat[i,]), error = function(e) NULL)
  if(!is.null(res)) c(res$locFDR, res$log10_BF) else c(NA, NA)
})
bfs <- t(bfs)

# Filter SNPs passing threshold and matching illustrative targets
target_snps <- c("rs1713411", "rs10812800", "rs12914272", "rs1482057", "rs878156")
sig_indices <- which(bfs[,1] < 0.05 & bfs[,2] > 1 & all_snps %in% target_snps)

# Run full MCMC runs on matching targets
for (i in sig_indices) {
  snp_name <- all_snps[i]
  post <- cpbayes_uncor(beta_mat[i,], se_mat[i,], 
                        Phenotypes = Traits, Variant = snp_name, 
                        MCMCiter = 5500, Burnin = 500)
  post_summ <- post_summaries(post)
  
  cat("\nSNP:", snp_name, " (locFDR =", round(bfs[i,1],4), ", log10_BF =", round(bfs[i,2],4), ")\n")
  print(post_summ$important_traits)
}
```

**Output:**
```text
SNP: rs10812800  (locFDR = 0.0247 , log10_BF = 1.1197 )
   traits   PPAj direction
1  Breast 0.4462  negative
2 Thyroid 0.9840  positive

SNP: rs1713411  (locFDR = 0.0285 , log10_BF = 1.0554 )
   traits   PPAj direction
1  Breast 0.2468  negative
2 Thyroid 0.9070  positive

SNP: rs878156  (locFDR = 0.0013 , log10_BF = 2.4122 )
   traits   PPAj direction
1  Breast 0.2166  negative
2 Thyroid 0.9714  positive

SNP: rs12914272  (locFDR = 0.0092 , log10_BF = 1.5534 )
   traits   PPAj direction
1  Breast 0.9970  negative
2 Thyroid 0.3658  positive

SNP: rs1482057  (locFDR = 0.0051 , log10_BF = 1.8093 )
   traits   PPAj direction
1  Breast 1.0000  negative
2 Thyroid 0.2766  positive
```

### 5. Save Results
Save MCMC details and posteriors to `results/cpbayes_results.rds`:

```R
saveRDS(
  list(
    analytic = result,
    cpbayes_model = result_CP,
    post_summary = PleioSumm
  ),
  here::here("results", "cpbayes_results.rds")
)

cat("CPBayes analysis complete. Results saved to results/cpbayes_results.rds\n")
```

**Output:**
```text
CPBayes analysis complete. Results saved to results/cpbayes_results.rds
```

## Description & Key Tasks

1. **Analytical Filter**: Employs `analytic_locFDR_BF_uncor()` to rapidly compute analytical local FDR (locFDR) and Bayes factors (BF) for all 3,707 SNPs. This serves as a fast pre-filtering step.
2. **Bayesian Model Fitting (`cpbayes_uncor`)**: Fits the Bayesian model using MCMC sampling (5,500 iterations, 500 burn-in) for the variants of interest to calculate PPAs and effect sizes.
3. **Posterior Summaries**: Extracts and prints trait-specific posterior probabilities of association (PPAj), posterior mean/SD of beta coefficients, and odds ratios.
4. **Save Results**: Outputs analytical and MCMC-based posterior statistics to `results/cpbayes_results.rds`.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/05_cpbayes.R")
```
