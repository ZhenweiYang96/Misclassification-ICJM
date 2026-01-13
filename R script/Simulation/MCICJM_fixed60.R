# Load function and package -----------------------------------------------
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
library(future)
source("R script/Functions/MCICJM fit_fixed.R")

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
    mcicjm <- mcicjm.fit(data = train.data, biopsy = train.biopsy, n.adapt = 5000,
                         n.burnin = 5000, n.iter = 5000, seed = seed[datanum],
                         sens = 0.6,
                         qp = 7, df_pspline = 12,
                         GLMM_fit = mm_inits)
    save(mcicjm,
         file = paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_",
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

