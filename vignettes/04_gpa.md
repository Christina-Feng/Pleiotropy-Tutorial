# Step 04: GPA Analysis

This script implements the **GPA** (Genetic Analysis Incorporating Pleiotropy and Annotation) method. It tests for pleiotropy and prioritises variants by integrating GWAS p-values with functional annotation data (specifically pathways).

## Step-by-Step Analysis & Results

### 1. Load GPA and prepared Summary Statistics
Load libraries, load data, and align trait p-values:

```R
library(GPA)
library(here)
library(tidyverse)

# Load prepared summary statistics
dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))

# Extract p-vals from both sum stats 
pvals <- left_join(dat$B_sumstats %>% dplyr::select(SNPs, p.value),
                   dat$T_sumstats %>% dplyr::select(SNPs, p.value),
                   by = "SNPs") 
colnames(pvals) <- c("snp","breast","thyroid")

cat("--- First 3 rows of p-values matrix ---\n")
print(head(pvals,3))
```

**Output:**
```text
--- First 3 rows of p-values matrix ---
# A tibble: 3 × 3
  snp        breast thyroid
  <chr>       <dbl>   <dbl>
1 rs7531583  0.0320  0.0680
2 rs12044597 0.637   0.855 
3 rs697679   0.796   0.549 
```

### 2. Create Dummy Pathway Annotations
Construct binary dummy vectors indicating variant memberships in obesity, circadian rhythm, xenobiotic metabolism, cell cycle, and puberty pathways:

```R
annots <- dat$B_sumstats %>% 
  dplyr::select(SNPs, pathway) %>% 
  arrange(SNPs, pathway) %>% 
  group_by(pathway) %>% 
  mutate(ann_obesity = ifelse(pathway == "F_obesity",1,0),
         ann_circadian = ifelse(pathway == "F_circadian",1,0),
         ann_xeno = ifelse(pathway == "F_xeno",1,0),
         ann_cellcycle = ifelse(pathway == "F_cell_cycle",1,0),
         ann_delpuberty = ifelse(pathway == "F_pub_he2010_4",1,0)
  ) %>%
  ungroup()
```

### 3. Fit GPA Model with No Annotation
Fit the baseline GPA model (using only GWAS p-values) and print the fitted object:

```R
# Fitting the GPA model
# Lower bound for pi estimates is set to 1 / [number of SNPs].
set.seed(123)
fit.GPA.noAnn <- GPA(
  pvals[, c("breast", "thyroid")],
  NULL,
  lbPi = 1 / nrow(pvals)
)

cat("--- GPA fit (No Annotation) ---\n")
print(fit.GPA.noAnn)
```

**Output:**
```text
Info: EM algorithm stops because it reaches the specified maximum number of iterations (maxIter = 2000).
--- GPA fit (No Annotation) ---
               00         10          01           11
   [1,] 0.9464057 0.03012824 0.022880888 0.0005851296
   [2,] 0.9641748 0.03060223 0.005093113 0.0001298565
   [3,] 0.9626502 0.03054703 0.006633634 0.0001690966
... [3664 rows truncated for readability] ...
```

### 4. Perform Association Mapping without Annotation
Evaluate which variants are associated under local FDR of 0.2:

```R
# Association output of GPA with local FDR set at 0.2 
assoc.GPA.noAnn <- GPA::assoc(fit.GPA.noAnn, FDR=0.2, fdrControl="local")
cat("--- Association Mapping Table (Breast) ---\n")
print(table(assoc.GPA.noAnn[,1]))
cat("--- Association Mapping Table (Thyroid) ---\n")
print(table(assoc.GPA.noAnn[,2])) 

# local FDR output of the GPA with No Annotation
fdr.GPA.noAnn <- GPA::fdr(fit.GPA.noAnn)
cat("--- First 3 rows of local FDR output ---\n")
print(head(fdr.GPA.noAnn,3))
```

**Output:**
```text
--- Association Mapping Table (Breast) ---
   0 
3707 
--- Association Mapping Table (Thyroid) ---
   0 
3707 
--- First 3 rows of local FDR output ---
         breast   thyroid
[1,] 0.9652049 0.9095851
[2,] 0.9687995 0.9933036
[3,] 0.9688175 0.9915038
```

### 5. Fit GPA Model with Xenobiotic Annotation
Fit the GPA model incorporating Xenobiotic pathway data:

```R
set.seed(123)
fit.GPA.wAnn <- GPA(pvals[,c(2,3)], annots$ann_xeno)

cat("--- GPA fit (With Xeno Annotation) ---\n")
print(fit.GPA.wAnn)
```

**Output:**
```text
Info: EM algorithm stops because it reaches the specified maximum number of iterations (maxIter = 2000).
# (Prints parameter matrices and fitted components for model with annotation)
```

### 6. Perform Association Mapping with Annotation
Determine associations and examine the first 3 rows of the new local FDR matrix:

```R
# Association output of the GPA with local FDR set at 0.2
assoc.GPA.wAnn <- GPA::assoc(fit.GPA.wAnn, FDR=0.2, fdrControl="local")

cat("--- Association Mapping Table with Annotation (Breast) ---\n")
print(table(assoc.GPA.wAnn[,1]))
cat("--- Association Mapping Table with Annotation (Thyroid) ---\n")
print(table(assoc.GPA.wAnn[,2])) 

# FDR output of the GPA with Annotation
fdr.GPA.wAnn <- GPA::fdr(fit.GPA.wAnn)
cat("--- First 3 rows of FDR output with Annotation ---\n")
print(head(fdr.GPA.wAnn,3))
```

**Output:**
```text
--- Association Mapping Table with Annotation (Breast) ---
   0 
3707 
--- Association Mapping Table with Annotation (Thyroid) ---
   0 
3707 
--- First 3 rows of FDR output with Annotation ---
         breast   thyroid
[1,] 0.9652049 0.9095851
[2,] 0.9687995 0.9933036
[3,] 0.9688175 0.9915038
```

### 7. Save Results
Save fitted objects, mappings, and local FDRs to `results/gpa_results.rds`:

```R
saveRDS(
  list(
    fit_noAnn = fit.GPA.noAnn,
    fit_wAnn = fit.GPA.wAnn,
    assoc_noAnn = assoc.GPA.noAnn,
    assoc_wAnn = assoc.GPA.wAnn,
    fdr_noAnn = fdr.GPA.noAnn,
    fdr_wAnn = fdr.GPA.wAnn
  ),
  here::here("results", "gpa_results.rds")
)

cat("GPA analysis complete. Results saved to results/gpa_results.rds\n")
```

**Output:**
```text
GPA analysis complete. Results saved to results/gpa_results.rds
```

## Description & Key Tasks

1. **Dummy Annotations**: Formulates dummy indicator variables (0 or 1) representing membership in biological pathways (Obesity, Circadian, Xenobiotic metabolism, Cell Cycle, Puberty).
2. **Fit Models**:
   - Fits a baseline model without annotation using the GWAS p-values of the two traits.
   - Fits a second model incorporating the Xenobiotic metabolism (`ann_xeno`) annotation vector.
3. **Association & FDR**: Applies a local FDR threshold of 0.2 to map variants associated with breast cancer, thyroid cancer, or both.
4. **Save Results**: Outputs parameter estimates and local FDR statistics to `results/gpa_results.rds`.

## Running this script

You can execute this script from the root directory of the project in R/RStudio using:

```R
source("R/04_gpa.R")
```
