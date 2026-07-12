# Step 06: GCPBayes Analysis

This script implements **GCPBayes** (Bayesian group-level pleiotropy association analysis) to analyse gene-level pleiotropic associations using multi-variant summary statistics. It fits three types of spike-and-slab models: Dirac Spike (DS), Continuous Spike (CS), and Hierarchical Spike (HS).

## Step-by-Step Analysis & Results

### 1. Load GCPBayes & BCAT1 Summary Statistics
Load packages and multi-variant dataset for gene *BCAT1* (containing 50 SNPs across 2 studies):

```R
library(GCPBayes)
library(here)
library(tidyverse)

load(here::here("data", "BCAT1_sumstats.RData"))

K <- 2 # Two Studies
m <- length(snpnames) %>% as.numeric() # 50 SNPs
```

### 2. Fit Dirac Spike (DS) Model
Fit the DS model (MCMC over 20,000 iterations) and examine the fitted criteria (PPA, BF, and local FDR):

```R
set.seed(321)
RES_DS <- DS(Betah, Sigmah, kappa0 = 0.5, sigma20 = 1, 
             m = m, K = K, niter = 20000, burnin = 10000, nthin = 2, nchains = 1, 
             a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename 
             )

cat("--- DS Model Criteria ---\n")
print(RES_DS$Criteria)
```

**Output:**
```text
--- DS Model Criteria ---
$`Name of Gene`
[1] "BCAT1"

$PPA
$PPA[[1]]
[1] 0.9998

$PPA[[2]]
[1] 1

$log10BF
[1] 9.193099

$lBFDR
[1] 2.136877e-10

$theta
[1] 0.9999647
```

### 3. Fit Continuous Spike (CS) Model
Initialise prior configurations and fit the CS model:

```R
# Initialise zinit for the CS prior
set.seed(321)
p_vals <- matrix(0, K, m) 
for (k in 1:K) {
  p_vals[k, ] <- 2*pnorm(-abs(Betah[[k]] / sqrt(diag(Sigmah[[k]]))))
}

zinit <- rep(0, K)
for (j in 1:K) {
  index <- 1:m
  PVALUE <- p.adjust(p_vals[j, ]) 
  SIGNALS <- index[PVALUE < 0.05] 
  modelf1 <- rep(0, m)
  modelf1[SIGNALS] <- 1 
  if (max(modelf1) == 1) (zinit[j] <- 1)
}

# Fit CS model
RES_CS <- CS(Betah, Sigmah,
             kappa0 = 0.5, tau20 = 1, zeta0 = zinit,
             m = m, K = K, niter = 20000, burnin = 10000, nthin = 2, nchains = 1, 
             a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, 
             sigma2 = 10^-3, snpnames = snpnames, genename = genename
             )

cat("--- CS Model Criteria ---\n")
print(RES_CS$Criteria)
```

**Output:**
```text
--- CS Model Criteria ---
$`Name of Gene`
[1] "BCAT1"

$PPA
$PPA[[1]]
[1] 0

$PPA[[2]]
[1] 0

$log10BF
[1] -Inf

$lBFDR
[1] 1

$theta
[1] 0
```

### 4. Fit Hierarchical Spike (HS) Model
Fit the HS model and print indicators and beta parameters for each study:

```R
set.seed(321)
RES_HS <- HS(Betah, Sigmah, kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, 
         s20 = 1, m = m, K = K, niter = 20000, burnin = 10000, 
         nthin = 1, nchains = 1, 
         a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, 
         c1 = 1, c2 = 1, e2 = 1, snpnames, genename 
) 

cat("--- HS Model Indicator ---\n")
print(RES_HS$Indicator$`Significant studies and Variable pleiotropic effect based on median`)

cat("--- HS Model Beta (Study 1) ---\n")
print(RES_HS$Summary$Beta[[1]])

cat("--- HS Model Beta (Study 2) ---\n")
print(RES_HS$Summary$Beta[[2]])
```

**Output:**
```text
--- HS Model Indicator ---
           Study 1 Study 2 Total Pleiotropic effect
rs7299477        0       0     0                 No
rs10743520       1       0     1                 No
rs11047648       1       0     1                 No
rs10219610       0       1     1                 No
rs4255624        1       1     2                Yes
rs7979723        0       0     0                 No
... [Prints indicator outputs for remaining SNPs] ...

--- HS Model Beta (Study 1) ---
   Name of SNP         Mean         SD    val2.5pc        Median  val97.5pc
1    rs7299477  0.012717940 0.04305323 -0.06512223  0.000000e+00 0.12187105
2   rs10743520  0.049739100 0.06899787 -0.03815849  2.705197e-02 0.22536835
3   rs11047648  0.035070030 0.06178084 -0.05124957  1.087359e-02 0.19417666
4   rs10219610  0.005460369 0.04748649 -0.09058505  0.000000e+00 0.12253191
5    rs4255624 -0.068163490 0.08480463 -0.27236231 -4.203058e-02 0.03465762
6    rs7979723  0.013071041 0.04995372 -0.07920700  0.000000e+00 0.13874668
... [Prints post-mean & SD for remaining Study 1 SNPs] ...

--- HS Model Beta (Study 2) ---
   Name of SNP          Mean         SD     val2.5pc        Median  val97.5pc
1    rs7299477 -0.0064997921 0.04464330 -0.114959651  0.0000000000 0.08439784
2   rs10743520  0.0072286316 0.04972994 -0.093104478  0.0000000000 0.13313267
3   rs11047648  0.0100586678 0.05094592 -0.086526073  0.0000000000 0.13907182
4   rs10219610  0.0419084407 0.06419654 -0.043715377  0.0187043419 0.20352962
5    rs4255624 -0.0782725228 0.08980144 -0.290599242 -0.0541653564 0.03073214
6    rs7979723  0.0029034554 0.05014317 -0.106155730  0.0000000000 0.12036795
... [Prints post-mean & SD for remaining Study 2 SNPs] ...
```

### 5. Save Results
Save fitted DS, CS, and HS models to `results/gcpbayes_results.rds`:

```R
saveRDS(
  list(
    ds_model = RES_DS,
    cs_model = RES_CS,
    hs_model = RES_HS
  ),
  here::here("results", "gcpbayes_results.rds")
)

cat("GCPBayes analysis complete. Results saved to results/gcpbayes_results.rds\n")
```

**Output:**
```text
GCPBayes analysis complete. Results saved to results/gcpbayes_results.rds
```

## Description & Key Tasks

1. **Spike-and-Slab Formulation**: Analyses the multi-variant gene-level summary statistics of the *BCAT1* gene across two traits.
2. **Three Modelling Approaches**:
   - **Dirac Spike (DS)**: Standard spike-and-slab prior where the spike is a point mass at zero. Runs for 20,000 MCMC iterations.
   - **Continuous Spike (CS)**: Approximates the point mass using a continuous normal prior with a very small variance.
   - **Hierarchical Spike (HS)**: Models study-specific and variant-specific parameters in a hierarchical structure to explore heterogeneous configurations.
3. **Save Results**: Outputs the final fitted objects (criteria, indicators, posterior summaries) to `results/gcpbayes_results.rds`.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/06_gcpbayes.R")
```
