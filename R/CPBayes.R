#################
#### CPBayes ####
#################

# Load the R package and import the data
library(CPBayes)
library(tidyverse)
B_sum_stats<- read_csv(here::here("data","B_sumstats.csv"))
T_sum_stats<- read_csv(here::here("data","T_sumstats.csv"))

# Find the pleiotropic association of a genetic variant
# For example: SNP rs12914272 - showed pleiotropic signals

rs12914272<- rbind(B_sum_stats %>% filter(SNPs == "rs12914272"), 
                   T_sum_stats %>% filter(SNPs == "rs12914272"))

# Step 1: Specifies Betas, SEs, traits and SNPs of interest
Beta_Hat <- rs12914272$estimate
SE <- rs12914272$std.error
Traits <- c("Breast","Thyroid")
SNP <- "rs12914272"

# Step 2: Computes locFDR and log10_BF
result <- analytic_locFDR_BF_uncor(Beta_Hat, SE)
result

# Step 3: Implement CPBayes for unrelated traits 
# Burin 500 samples, and draw 5500 posterior samples 
result_CP <- cpbayes_uncor(Beta_Hat, SE, 
                           Phenotypes = Traits, Variant = SNP, 
                           MCMCiter = 5500, Burnin = 500)
result_CP$important_traits

# Summary of posterior probability of associations (PPAj)
PleioSumm <- post_summaries(result_CP, level = 0.05) 

PleioSumm$important_traits
PleioSumm$poste_summary_beta
PleioSumm$poste_summary_OR


# for multiple SNPs, loop iterations (do not execute)
for(i in 1:K){
  post <- cpbayes_uncor(Beta_Hat[i,], SE[i,], Phenotypes = Traits, 
                        Variant = SNPs[i,], MCMCiter = 5500, Burnin = 500)
  PleioSumm <- post_summaries(post) 
  print(PleioSumm$important_traits)
} 