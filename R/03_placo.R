###############
#### PLACO ####
###############
library(devtools)
library(here)
library(tidyverse)
library(MASS)

# Source PLACO utility functions directly from GitHub
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")

# Load prepared summary statistics
dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))

# Step 01 Initial correlation investigation 
Z.matrix <- dat$z_matrix
P.matrix <- dat$p_matrix
CorZ <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
cat("--- Z-score Correlation (CorZ) ---\n")
print(CorZ)

# Step 02 Use 'placo' function if CorZ is less than than 0.5, 
#         indicating uncorrelated traits.
# For correlated traits (CorZ greater than equal to 0.5), 
#         use 'placo.plus' in Step 2 instead.

# Step 1 Estimate the variance parameters
VarZ <- var.placo(Z.matrix, P.matrix, p.threshold=1e-4)

# Step 2 Pleiotropic association test and its asymptotic p-value
PLACO_res <- sapply(1:nrow(Z.matrix), 
                    function(i) placo(Z=Z.matrix[i,], VarZ=VarZ)) #placo.plus for correlated traits
res <- t(PLACO_res) 

# Step 3 Output of specific variant, for example, SNP rs7531583 with index 1
cat("--- Output of specific variant (rs7531583, index 1) ---\n")
cat("T.placo:", res[1,]$T.placo, "\n")
cat("p.placo:", res[1,]$p.placo, "\n")

PLACO_res_df <- res %>% as.data.frame() %>% 
  mutate(SNPs = dat$annot$snp) %>% 
  left_join(dat$B_sumstats, by = join_by(SNPs == SNPs)) 

# Display selected SNPs reported in the manuscript
target_snps <- c("rs7531583", "rs12049503", "rs6676419")

cat("--- PLACO results for selected SNPs ---\n")
print(
  PLACO_res_df %>%
    filter(SNPs %in% target_snps) %>%
    arrange(match(SNPs, target_snps)) %>%
    dplyr::select(SNPs, T.placo, p.placo, chr, pos, minor, major, gene, pathway)
)

# Save results
saveRDS(
  list(
    CorZ = CorZ,
    VarZ = VarZ,
    results_df = PLACO_res_df
  ),
  here::here("results", "placo_results.rds")
)

cat("PLACO analysis complete. Results saved to results/placo_results.rds\n")


