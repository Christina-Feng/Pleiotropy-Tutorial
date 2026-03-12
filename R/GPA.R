#############
#### GPA ####
#############

# library(devtools)
# if(!requireNamespace("BiocManager", quietly = TRUE)) + install.packages("BiocManager") 
# BiocManager::install("GPA")
# install_github("dongjunchung/GPA", force = TRUE)
library(GPA)

# Load Data
library(tidyverse)
B_sum_stats<- read_csv(here::here("data","B_sumstats.csv"))
T_sum_stats<- read_csv(here::here("data","T_sumstats.csv"))

# Extract p-vals from both sum stats 
pvals <- left_join(B_sum_stats %>% dplyr::select(1,5),
                   T_sum_stats %>% dplyr::select(1,5),
                   by = "SNPs") 
colnames(pvals) <- c("snp","breast","thyroid")
head(pvals,3)

# Create dummy variables for pathways 
# x 5 detected pleiotropic effect in other methods
annots <- B_sum_stats %>% 
  dplyr::select(SNPs, pathway) %>% 
  arrange(SNPs, pathway) %>% 
  group_by(pathway) %>% 
  mutate(ann_obesity = ifelse(pathway == "F_obesity",1,0),
         ann_circadian = ifelse(pathway == "F_circadian",1,0),
         ann_xeno = ifelse(pathway == "F_xeno",1,0),
         ann_cellcycle = ifelse(pathway == "F_cell_cycle",1,0),
         ann_delpuberty = ifelse(pathway == "F_pub_he2010_4",1,0)
  )

#############################################
######## Fit GPA with No Annotation #########
#############################################

# Fitting the GPA model
# Lower bound for pi estimates is set to 1 / [number of SNPs].

fit.GPA.noAnn <- GPA(pvals[,c(2,3)], NULL)
fit.GPA.noAnn

# Association output of the GPA with a local FDR set at 0.5 (exploratory) 
assoc.GPA.noAnn <- GPA::assoc(fit.GPA.noAnn, FDR=0.2, fdrControl="local" )
table(assoc.GPA.noAnn[,1])
table(assoc.GPA.noAnn[,2]) #FDR = 0.5, got >=3 SNPs associated with thyroid

# local FDR output of the GPA with No Annotation
fdr.GPA.noAnn <- GPA::fdr(fit.GPA.noAnn)
head(fdr.GPA.noAnn,3)

###############################################################
######## GPA: With one vector of annotation (xeno) data #######
###############################################################

# Fitting the GPA model
fit.GPA.wAnn <- GPA(pvals[,c(2,3)], annots$ann_xeno)
fit.GPA.wAnn

# estimates(fit.GPA.wAnn)
# se(fit.GPA.wAnn)

##########################################
######## GPA: Association Mapping  #######
##########################################

#‘assoc’ method returns a binary matrix indicating association of each SNP, 
# where 1 indicates that a SNP is associated with the phenotype and 0 otherwise.

# Association output of the GPA with a global FDR set at 0.05
# assoc.GPA.wAnn <- GPA::assoc(fit.GPA.wAnn, FDR=0.05, fdrControl="global")

# Association output of the GPA with a local FDR set at 0.5
assoc.GPA.wAnn <- GPA::assoc(fit.GPA.wAnn, FDR=0.5, fdrControl="local")

table(assoc.GPA.wAnn[,1])
table(assoc.GPA.wAnn[,2]) #fdr=0.5 indicates 3 SNPs associated with Thyroid

# FDR output of the GPA with Annotation
fdr.GPA.wAnn<- GPA::fdr(fit.GPA.wAnn)
head(fdr.GPA.wAnn,3)

