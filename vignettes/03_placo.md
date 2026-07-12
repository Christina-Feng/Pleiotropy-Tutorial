# Step 03: PLACO Analysis

This script applies **PLACO** (Pleiotropic Association Test under Composite Null Hypothesis) to test for pleiotropic effects between breast cancer and thyroid cancer using summary statistics.

## Step-by-Step Analysis & Results

### 1. Load Libraries and Source PLACO from GitHub
Load required dependencies, source the official PLACO R script from GitHub, and read the prepared summary statistics:

```R
library(devtools)
library(here)
library(tidyverse)
library(MASS)

# Source PLACO utility functions directly from GitHub
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")

# Load prepared summary statistics
dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))
```

### 2. Check Z-score Pearson Correlation
Isolate Z-scores and p-values to investigate the correlation structure of the traits:

```R
Z.matrix <- dat$z_matrix
P.matrix <- dat$p_matrix
CorZ <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)

cat("--- Z-score Correlation (CorZ) ---\n")
print(CorZ)
```

**Output:**
```text
--- Z-score Correlation (CorZ) ---
            Breast    Thyroid
Breast  1.00000000 0.02076901
Thyroid 0.02076901 1.00000000
```
*(Since CorZ is $< 0.5$, the traits are uncorrelated, allowing standard `placo` composite null testing.)*

### 3. Estimate Variance Parameters & Perform Association Tests
Estimate the variance matrix under the composite null and test each variant:

```R
# Step 1: Estimate the variance parameters
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)

# Step 2: Pleiotropic association test and its asymptotic p-value
PLACO_res <- sapply(1:nrow(Z.matrix), 
                    function(i) placo(Z=Z.matrix[i,], VarZ=VarZ))
res <- t(PLACO_res) 
```

### 4. Output Statistics for Specific Variant (rs7531583)
Examine statistics for the first variant in the dataset as an example:

```R
cat("--- Output of specific variant (rs7531583, index 1) ---\n")
cat("T.placo:", res[1,]$T.placo, "\n")
cat("p.placo:", res[1,]$p.placo, "\n")
```

**Output:**
```text
--- Output of specific variant (rs7531583, index 1) ---
T.placo: 3.914022 
p.placo: 0.0067312 
```

### 5. Display Summary of Illustrative Target SNPs
Merge test statistics with biological annotations and display statistics for selected target SNPs reported in the manuscript:

```R
PLACO_res_df <- res %>% as.data.frame() %>% 
  mutate(SNPs = dat$annot$snp) %>% 
  left_join(dat$B_sumstats, by = join_by(SNPs == SNPs)) 

target_snps <- c("rs7531583", "rs12049503", "rs6676419")

cat("--- PLACO results for selected SNPs ---\n")
print(
  PLACO_res_df %>%
    filter(SNPs %in% target_snps) %>%
    arrange(match(SNPs, target_snps)) %>%
    dplyr::select(SNPs, T.placo, p.placo, chr, pos, minor, major, gene, pathway)
)
```

**Output:**
```text
--- PLACO results for selected SNPs ---
        SNPs   T.placo     p.placo chr      pos minor major   gene            pathway
1  rs7531583  3.914022   0.0067312   1  1706160     A     G   NADK F_tobacco_hsa00760
2 rs12049503 -3.675333 0.008784198   1 46021556     A     C AKR1A1             F_xeno
3  rs6676419  2.380087  0.03836835   1 66052209     A     G   LEPR          F_obesity
```

### 6. Save Results
Save final estimates and data frame containing all statistics to `results/placo_results.rds`:

```R
saveRDS(
  list(
    CorZ = CorZ,
    VarZ = VarZ,
    results_df = PLACO_res_df
  ),
  here::here("results", "placo_results.rds")
)

cat("PLACO analysis complete. Results saved to results/placo_results.rds\n")
```

**Output:**
```text
PLACO analysis complete. Results saved to results/placo_results.rds
```

## Description & Key Tasks

1. **Source Code**: Fetches the core PLACO functions directly from the official GitHub repository.
2. **Correlation Estimation**: Computes the Pearson correlation coefficient between the Z-scores of the two traits using variants that are not strongly associated with either trait.
3. **Variance Parameter Estimation**: Estimates variance and covariance parameters ($\text{Var}(Z)$ matrix) under the null hypothesis.
4. **Pleiotropic Association Test**: Computes the PLACO test statistic ($T_{\text{placo}}$) and corresponding p-value ($p_{\text{placo}}$) for each variant.
5. **Save Results**: Outputs the final data frame containing statistics, variant annotations, and genomic coordinates to `results/placo_results.rds`.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/03_placo.R")
```
