#############
#### GPA ####
#############

# library(devtools)
# if(!requireNamespace("BiocManager", quietly = TRUE)) + install.packages("BiocManager") 
# BiocManager::install("GPA")
# install_github("dongjunchung/GPA", force = TRUE)
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

# Create dummy variables for pathways 
# x 5 detected pleiotropic effect in other methods
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

#############################################
######## Fit GPA with No Annotation #########
#############################################

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

# Association output of the GPA with a local FDR set at 0.2 
assoc.GPA.noAnn <- GPA::assoc(fit.GPA.noAnn, FDR=0.2, fdrControl="local" )
cat("--- Association Mapping Table (Breast) ---\n")
print(table(assoc.GPA.noAnn[,1]))
cat("--- Association Mapping Table (Thyroid) ---\n")
print(table(assoc.GPA.noAnn[,2])) 

# local FDR output of the GPA with No Annotation
fdr.GPA.noAnn <- GPA::fdr(fit.GPA.noAnn)
cat("--- First 3 rows of local FDR output ---\n")
print(head(fdr.GPA.noAnn,3))

###############################################################
######## GPA: With one vector of annotation (xeno) data #######
###############################################################

# Fitting the GPA model
set.seed(123)
fit.GPA.wAnn <- GPA(pvals[,c(2,3)], annots$ann_xeno)
cat("--- GPA fit (With Xeno Annotation) ---\n")
print(fit.GPA.wAnn)

# estimates(fit.GPA.wAnn)
# se(fit.GPA.wAnn)

##########################################
######## GPA: Association Mapping  #######
##########################################

#‘assoc’ method returns a binary matrix indicating association of each SNP, 
# where 1 indicates that a SNP is associated with the phenotype and 0 otherwise.

# Association output of the GPA with a local FDR set at 0.2
assoc.GPA.wAnn <- GPA::assoc(fit.GPA.wAnn, FDR=0.2, fdrControl="local")

cat("--- Association Mapping Table with Annotation (Breast) ---\n")
print(table(assoc.GPA.wAnn[,1]))
cat("--- Association Mapping Table with Annotation (Thyroid) ---\n")
print(table(assoc.GPA.wAnn[,2])) 

# FDR output of the GPA with Annotation
fdr.GPA.wAnn <- GPA::fdr(fit.GPA.wAnn)
cat("--- First 3 rows of FDR output with Annotation ---\n")
print(head(fdr.GPA.wAnn,3))

# Save results
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
