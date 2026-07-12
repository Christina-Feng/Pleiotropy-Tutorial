############################################################
# 01_prepare_data.R
# Load and prepare example GWAS summary statistics
############################################################

library(tidyverse)

# Load raw summary statistics
load(here::here("data", "GWAS_sumstats.RData"))

traits_lab <- c("Breast", "Thyroid")

snps <- annot$snp %>% as.character()

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

# Save the structured dataset
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

# Print a preview of the structured matrices for the audience
cat("--- First 3 rows of Aligned Estimates (beta_hat) ---\n")
print(head(beta_hat, 3))
cat("--- First 3 rows of Aligned Z-scores (z_matrix) ---\n")
print(head(z_matrix, 3))

cat("Data preparation complete. Saved to results/prepared_summary_statistics.rds\n")
