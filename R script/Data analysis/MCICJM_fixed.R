library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
library(future)
source("R script/Functions/MCICJM fit_fixed.R")
load("Data/biopsy_sub.RData")
load("Data/pass_id.RData")
load("Data/pass.RData")
load("Output/Data analysis/inits/mixed model_inits.RData")

plan(multisession, workers = 3)

# Sens = 100% -------------------------------------------------------------
mcicjm.model <- mcicjm.fit(data = pass, biopsy = biopsy, n.adapt = 10000, 
                           n.burnin = 10000, n.iter = 10000, seed = 2023,
                           sens = 1, GLMM_fit = mm_inits)
save(mcicjm.model,
     file = "Output/Data analysis/MCICJM_1.RData")

# Sens = 80% -------------------------------------------------------------
mcicjm.model <- mcicjm.fit(data = pass, biopsy = biopsy, n.adapt = 10000, 
                           n.burnin = 10000, n.iter = 10000, seed = 2023,
                           sens = 0.8, GLMM_fit = mm_inits)
save(mcicjm.model,
     file = "Output/Data analysis/MCICJM_0.8.RData")

# Sens = 60% -------------------------------------------------------------
mcicjm.model <- mcicjm.fit(data = pass, biopsy = biopsy, n.adapt = 10000, 
                           n.burnin = 10000, n.iter = 10000, seed = 2023,
                           sens = 0.6, GLMM_fit = mm_inits)
save(mcicjm.model,
     file = "Output/Data analysis/MCICJM_0.6.RData")


# Sens = 75% (for simulation) ---------------------------------------------
mcicjm.model <- mcicjm.fit(data = pass, biopsy = biopsy, n.adapt = 10000, 
                           n.burnin = 10000, n.iter = 10000, seed = 2023,
                           sens = 0.75, GLMM_fit = mm_inits)
save(mcicjm.model,
     file = "Output/Data analysis/MCICJM_0.75.RData")