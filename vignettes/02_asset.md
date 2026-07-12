# Step 02: ASSET Analysis

This script applies the subset-based meta-analysis method **ASSET** (association analysis of heterogeneous traits/subtypes) to explore pleiotropic association patterns across breast and thyroid cancers.

## Step-by-Step Analysis & Results

### 1. Load ASSET & Structured GWAS Statistics
First, load the required libraries and the structured summary statistics:

```R
library(ASSET)
library(tidyverse)
library(here)

dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))
```

### 2. Format Variables & Configure Correlation Structure
Extract components required by ASSET and define the diagonal correlation lists under the assumption of independent studies:

```R
# Extract inputs
snps <- dat$snps
traits.lab <- dat$traits_lab
beta.hat <- dat$beta_hat
sigma.hat <- dat$sigma_hat
ncase <- dat$N11
ncntl <- dat$N00

# Set up correlation list
k <- length(traits.lab)
cor <- list(
  N11 = diag(ncase, nrow = k, ncol = k),
  N00 = diag(ncntl, nrow = k, ncol = k),
  N10 = matrix(0, nrow = k, ncol = k)
)
```

### 3. Run Subset Meta-Analysis
Perform subset search using the Discrete Local Maximum approach (DLM) to adjust for multi-trait testing:

```R
res <- h.traits(
  snp.vars = snps, traits.lab = traits.lab, beta.hat = beta.hat, sigma.hat = sigma.hat,
  ncase = ncase, ncntl = ncntl,
  cor = cor, cor.numr = FALSE,
  search = NULL, side = 2, meta = TRUE,
  zmax.args = NULL, meth.pval = "DLM"
)
```

### 4. Summarise Meta-Analysis Results
Print the first 3 rows of each meta-analysis summary section:

* **Meta Summary (Overall association):**
```R
cat("--- Meta Summary ---\n")
print(head(h.summary(res)$Meta, 3))
```

**Output:**
```text
--- Meta Summary ---
         SNP      Pvalue    OR CI.low CI.high
1  rs7531583 0.005022145 0.868  0.864   0.873
2 rs12044597 0.836495764 1.009  1.005   1.012
3   rs697679 0.805849716 1.012  1.007   1.017
```

* **Subset 1-sided Summary (Directional subset search):**
```R
cat("--- Subset 1-sided ---\n")
print(head(h.summary(res)$Subset.1sided, 3))
```

**Output:**
```text
--- Subset 1-sided ---
         SNP    Pvalue    OR CI.low  CI.high          Pheno
1  rs7531583 0.0130482 0.868  0.777    0.971 Breast,Thyroid
2 rs12044597 1.0000000 1.028  0.000      Inf         Breast
3   rs697679 0.9929280 1.042  0.000 8885.827        Thyroid
```

* **Subset 2-sided Summary (Heterogeneous subset search):**
```R
cat("--- Subset 2-sided ---\n")
print(head(h.summary(res)$Subset.2sided, 3))
```

**Output:**
```text
--- Subset 2-sided ---
         SNP     Pvalue  Pvalue.1   Pvalue.2  OR.1 CI.low.1 CI.high.1  OR.2  CI.low.2 CI.high.2 Pheno.1        Pheno.2
1  rs7531583 0.01209296 1.0000000 0.01209296    NA       NA        NA 0.868     0.778     0.970         Breast,Thyroid
2 rs12044597 0.87571060 0.6371509 0.85514351 1.028    0.916     1.153 0.989     0.881     1.111  Breast        Thyroid
3   rs697679 0.79912068 0.5493094 0.79634187 1.042    0.911     1.191 0.982     0.857     1.125 Thyroid         Breast
```

### 5. Save Results
Save the meta-analysis outputs to `results/asset_results.rds`:

```R
asset_summary <- list(
  meta = h.summary(res)$Meta,
  subset_1sided = h.summary(res)$Subset.1sided,
  subset_2sided = h.summary(res)$Subset.2sided
)
saveRDS(asset_summary, here::here("results", "asset_results.rds"))

cat("ASSET analysis complete. Results saved to results/asset_results.rds\n")
```

**Output:**
```text
ASSET analysis complete. Results saved to results/asset_results.rds
```

## Description & Key Tasks

1. **Format Input Variables**: Re-structures the prepared summary data into vectors and matrices expected by ASSET (SNPs, traits, beta, standard errors, case/control sample sizes).
2. **Correlation Matrix Setup**: Constructs the correlation structure. In this case, studies are assumed independent, so diagonal matrix formulations are used.
3. **Subset Search (`h.traits`)**: Applies the `h.traits()` function to search all possible subsets of traits for association, utilising the **Discrete Local Maximum (DLM)** approach to adjust for multi-trait testing.
4. **Save Results**: Outputs the meta-analysis summary tables to `results/asset_results.rds`.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/02_asset.R")
```
