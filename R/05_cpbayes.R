#################
#### CPBayes ####
#################

# Load the R package and prepared summary statistics
library(CPBayes)
library(here)
library(tidyverse)

dat <- readRDS(here::here("results", "prepared_summary_statistics.rds"))

# Find the pleiotropic association of a genetic variant
# For example: SNP rs12914272 - showed pleiotropic signals

rs12914272 <- rbind(dat$B_sumstats %>% filter(SNPs == "rs12914272"), 
                    dat$T_sumstats %>% filter(SNPs == "rs12914272"))

# Step 1: Specifies Betas, SEs, traits and SNPs of interest
Beta_Hat <- rs12914272$estimate
SE <- rs12914272$std.error
Traits <- c("Breast","Thyroid")
SNP <- "rs12914272"

# Step 2: Computes locFDR and log10_BF
result <- analytic_locFDR_BF_uncor(Beta_Hat, SE)
cat("--- Analytic locFDR and log10_BF ---\n")
print(result)

# Step 3: Implement CPBayes for unrelated traits 
# Burin 500 samples, and draw 5500 posterior samples 
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

# For multiple SNPs: Run CPBayes SNP-by-SNP for all SNPs in the dataset
# To be computationally efficient, we filter variants using working criteria: 
# locFDR < 0.05 and log10(BF) > 1

cat("\n--- Running CPBayes for all SNPs (Threshold: locFDR < 0.05 & log10(BF) > 1) ---\n")
all_snps <- dat$snps
beta_mat <- dat$beta_hat
se_mat <- dat$sigma_hat

# Compute fast analytical locFDR and log10_BF for all 3,707 SNPs
bfs <- sapply(1:length(all_snps), function(i) {
  res <- tryCatch(analytic_locFDR_BF_uncor(beta_mat[i,], se_mat[i,]), error = function(e) NULL)
  if(!is.null(res)) c(res$locFDR, res$log10_BF) else c(NA, NA)
})
bfs <- t(bfs) # Matrix of 2 columns: locFDR, log10_BF

# Filter SNPs passing threshold and matching the illustrative SNPs
target_snps <- c("rs1713411", "rs10812800", "rs12914272", "rs1482057", "rs878156")
sig_indices <- which(bfs[,1] < 0.05 & bfs[,2] > 1 & all_snps %in% target_snps)

for (i in sig_indices) {
  snp_name <- all_snps[i]
  post <- cpbayes_uncor(beta_mat[i,], se_mat[i,], 
                        Phenotypes = Traits, Variant = snp_name, 
                        MCMCiter = 5500, Burnin = 500)
  post_summ <- post_summaries(post)
  
  cat("\nSNP:", snp_name, " (locFDR =", round(bfs[i,1],4), ", log10_BF =", round(bfs[i,2],4), ")\n")
  print(post_summ$important_traits)
}

# Save results
saveRDS(
  list(
    analytic = result,
    cpbayes_model = result_CP,
    post_summary = PleioSumm
  ),
  here::here("results", "cpbayes_results.rds")
)

cat("CPBayes analysis complete. Results saved to results/cpbayes_results.rds\n")

