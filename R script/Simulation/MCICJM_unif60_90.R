# Load function and package -----------------------------------------------
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
library(future)
source("R script/Functions/MCICJM fit_prior.R")

# Preparation -------------------------------------------------------------
load("Output/Simulation/Seeds/DGM_seed.RData")
plan(
  list(
    tweak(multisession, workers = 5),
    tweak(multisession, workers = 3)
  )
)


out <- lapply(1:200, function(datanum) {
  future({
    
    #load data file
    load(paste0("Output/Simulation/Training sets/trainbiopsy_",
                datanum, ".RData"))
    load(paste0("Output/Simulation/Training sets/traindata_",
                datanum, ".RData"))
    load(paste0("Output/Simulation/Fitted mixed model/MM_inits_",
                datanum, ".RData"))
    
    # start modelling
    mcicjm <- mcicjm.fit(data = train.data, biopsy = train.biopsy, n.adapt = 10000,
                         n.burnin = 10000, n.iter = 5000, seed = seed[datanum],
                         prior = "unif", sens_a = 0.6, sens_b = 0.9,
                         qp = 7, df_pspline = 5,
                         GLMM_fit = mm_inits)
    save(mcicjm,
         file = paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_",
                       datanum, ".RData"))
    #
    plan(
      list(
        tweak(multisession, workers = 5),
        tweak(multisession, workers = 3)
      )
    )
    
  })
})
res <- lapply(out, future::value)
