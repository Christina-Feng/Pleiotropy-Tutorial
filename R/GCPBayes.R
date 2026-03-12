##################
#### GCPBayes ####
##################

# Load the R package and import the data
# install.packages("GCPBayes")
library(GCPBayes)
library(tidyverse)
load("BCAT1_sumstats.RData")

# Define Beta vectors and Diagonal Error matrices (m by m dimension)
set.seed(321)
K <- 2 # Two Studies
m <- length(snpnames) %>% as.numeric() # 50 SNPs

########################
### DS - Dirac Spike ###
########################
RES_DS <- DS(Betah, Sigmah, kappa0 = 0.5, sigma20 = 1, 
             m = m, K = K, niter = 20000, burnin = 10000, nthin = 2, nchains = 1, 
             a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, snpnames, genename 
             )

# Criteria log10BF, locFDR and theta
RES_DS$Criteria

#############################
### CS - Continuous Spike ###
#############################

# Initialise zinit for the CS prior
set.seed(321)
p_vals <- matrix(0, K, m) 
for (k in 1:K) {
  p_vals[k, ] <- 2*pnorm(-abs(Betah[[k]] / 
                                sqrt(diag(Sigmah[[k]]))))
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

RES_CS <- CS(Betah, Sigmah,
             kappa0 = 0.5, tau20 = 1, zeta0 = zinit,
             m = m, K = K, niter = 20000, burnin = 10000, nthin = 2, nchains = 1, 
             a1 = 0.1, a2 = 0.1, c1 = 0.1, c2 = 0.1, 
             sigma2 = 10^-3, snpnames = snpnames, genename = genename
             )

# Criteria log10BF, locFDR & theta
RES_CS$Criteria

#############################
### HS Hierarchical Spike ###
#############################
RES_HS <- HS(Betah, Sigmah, kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, 
         s20 = 1, m = m, K = K, niter = 20000, burnin = 10000, 
         nthin = 1, nchains = 1, 
         a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, 
         c1 = 1, c2 = 1, e2 = 1, snpnames, genename 
) 

RES_HS$Indicator

