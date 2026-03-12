###############
#### PLACO ####
###############
library(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")

# Load data
library(tidyverse)
library(MASS)
B_sum_stats<- read_csv(here::here("data","B_sumstats.csv"))
T_sum_stats<- read_csv(here::here("data","T_sumstats.csv"))

# Step 01 Initial correlation investigation 
Z.matrix<- as.matrix(cbind(B_sum_stats$statistic,
                           T_sum_stats$statistic))
P.matrix <- as.matrix(cbind(B_sum_stats$p.value,
                            T_sum_stats$p.value))
CorZ <- cor.pearson(Z.matrix, P.matrix, p.threshold=1e-4)
CorZ

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
res[1,]$T.placo
res[1,]$p.placo

PLACO_res <- res %>% as.data.frame() %>% 
  mutate(SNPs = B_sum_stats$SNPs) %>% 
  left_join(B_sum_stats, by = join_by(SNPs == SNPs)) 

head(PLACO_res,3)
