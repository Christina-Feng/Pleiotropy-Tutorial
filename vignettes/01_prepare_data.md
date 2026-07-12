# Step 01: Data Preparation

This script loads the raw GWAS summary statistics, aligns and extracts the relevant matrices (estimates, standard errors, Z-statistics, and p-values) for Breast and Thyroid cancer, and packages them into a structured dataset.

## Step-by-Step Analysis & Results

### 1. Load Libraries and Import Summary Statistics
First, load the `tidyverse` library and import the raw summary statistics:

```R
library(tidyverse)

# Load raw summary statistics
load(here::here("data", "GWAS_sumstats.RData"))

traits_lab <- c("Breast", "Thyroid")
snps <- annot$snp %>% as.character()
```

### 2. Align Statistics and Extract Matrices
Extract effect size estimates ($\hat{\beta}$), standard errors ($\hat{\sigma}$), Z-scores, and p-values for both breast and thyroid cancer:

```R
# Extract beta and sigma matrices
beta_hat <- cbind(B_sumstats$estimate, T_sumstats$estimate) %>% as.matrix()
sigma_hat <- cbind(B_sumstats$std.error, T_sumstats$std.error) %>% as.matrix()

# Extract Z-statistic and p-value matrices
z_matrix <- cbind(
  Breast = B_sumstats$statistic,
  Thyroid = T_sumstats$statistic
) %>% as.matrix()

p_matrix <- cbind(
  Breast = B_sumstats$p.value,
  Thyroid = T_sumstats$p.value
) %>% as.matrix()
```

### 3. Preview Aligned Data
Print the first 3 rows of the aligned estimates and Z-scores to verify the layout:

```R
cat("--- First 3 rows of Aligned Estimates (beta_hat) ---\n")
print(head(beta_hat, 3))
```

**Output:**
```text
--- First 3 rows of Aligned Estimates (beta_hat) ---
            [,1]        [,2]
[1,] -0.15333205 -0.12903026
[2,]  0.02767578 -0.01080159
[3,] -0.01792820  0.04093252
```

```R
cat("--- First 3 rows of Aligned Z-scores (z_matrix) ---\n")
print(head(z_matrix, 3))
```

**Output:**
```text
--- First 3 rows of Aligned Z-scores (z_matrix) ---
         Breast    Thyroid
[1,] -2.1443557 -1.8252670
[2,]  0.4716861 -0.1825597
[3,] -0.2580843  0.5987953
```

### 4. Save Structured Dataset
Save the final aligned list to `results/prepared_summary_statistics.rds` for subsequent steps:

```R
saveRDS(
  list(
    snps = snps,
    traits_lab = traits_lab,
    beta_hat = beta_hat,
    sigma_hat = sigma_hat,
    z_matrix = z_matrix,
    p_matrix = p_matrix,
    B_sumstats = B_sumstats,
    T_sumstats = T_sumstats,
    annot = annot,
    N11 = N11,
    N00 = N00
  ),
  here::here("results", "prepared_summary_statistics.rds")
)

cat("Data preparation complete. Saved to results/prepared_summary_statistics.rds\n")
```

**Output:**
```text
Data preparation complete. Saved to results/prepared_summary_statistics.rds
```

## Description & Key Tasks

1. **Load Raw Data**: Imports the summary statistics data (`data/GWAS_sumstats.RData`).
2. **Align Statistics**: Extracts effect size estimates ($\hat{\beta}$), standard errors ($\hat{\sigma}$), Z-scores, and p-values for both breast cancer and thyroid cancer datasets.
3. **Save Output**: Saves a structured list containing all these matrices and annotations as `results/prepared_summary_statistics.rds`, which will serve as the input for subsequent analyses.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/01_prepare_data.R")
```
