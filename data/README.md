# Data Folder

This folder contains the example datasets used in the tutorial, stored as R Data (`.RData`) files.

## Files

### 1. `GWAS_sumstats.RData`
Contains aligned variant-level Genome-Wide Association Study (GWAS) summary statistics for breast and thyroid cancer.

Loading this file (`load("data/GWAS_sumstats.RData")`) loads the following objects:
* **`B_sumstats`**: A dataframe with breast cancer summary statistics for 3,707 SNPs. Key columns:
  * `SNPs`: The rsID of the variant.
  * `estimate`: Estimated Effect size ($\hat{\beta}$).
  * `std.error`: Standard error ($\hat{\sigma}$).
  * `statistic`: Z-statistic.
  * `p.value`: p-value.
  * `chr`, `pos`, `minor`, `major`, `gene`, `pathway`: Genomic coordinates, alleles, and pathway annotations.
* **`T_sumstats`**: A dataframe with thyroid cancer summary statistics (same structure as `B_sumstats`).
* **`annot`**: A dataframe with genes and biological pathway annotations for the 3,707 SNPs.
* **`N11`**: Case sample sizes (Breast: 1,125; Thyroid: 1,129).
* **`N00`**: Control sample sizes (Breast: 1,172; Thyroid: 1,174).

---

### 2. `BCAT1_sumstats.RData`
Contains gene-level joint summary statistics and covariance structure for the *BCAT1* gene, used for the GCPBayes analysis (Step 06).

Loading this file (`load("data/BCAT1_sumstats.RData")`) loads the following objects:
* **`Betah`**: A list of 2 vectors containing joint estimated effect sizes ($\hat{\beta}$) for 50 SNPs in the *BCAT1* gene region for Study 1 (Breast) and Study 2 (Thyroid).
* **`Sigmah`**: A list of 2 covariance matrices $\widehat{\text{Cov}}(\hat{\beta})$ (50 by 50) representing the covariance structure among the 50 SNPs in both studies.
* **`snpnames`**: rsIDs of the 50 SNPs.
* **`genename`**: The gene symbol (`"BCAT1"`).
