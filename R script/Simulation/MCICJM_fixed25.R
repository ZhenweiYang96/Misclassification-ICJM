# Load function and package -----------------------------------------------
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
library(future)
source("R script/Functions/MCICJM fit_fixed.R")
ncore <- detectCores() - 1
ncore_div3 <- floor(ncore/3)


# Preparation -------------------------------------------------------------
load("Output/Simulation/Seeds/DGM_seed.RData")
plan(
  list(
    tweak(multisession, workers = ncore_div3),
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
                         sens = 0.25,
                         qp = 7, df_pspline = 12,
                         GLMM_fit = mm_inits)
    save(mcicjm,
         file = paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_",
                       datanum, ".RData"))
    #
    plan(
      list(
        tweak(multisession, workers = ncore_div3),
        tweak(multisession, workers = 3)
      )
    )
    
  })
})
res <- lapply(out, future::value)

# Examine the convergence -------------------------------------------------
for (i in 1:200) {
  load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_",
              i, ".RData"))
  dir.create(paste0("Output/Simulation/MCICJM/Fixed_25/convergence/", i))
  mcmcplots::mcmcplot(mcicjm$mcmc, 
                      dir = paste0("Output/Simulation/MCICJM/Fixed_25/convergence/", i,"/"),
                      filename = paste0("mcmcplot_", i), 
                      parms = c("alpha", "betaL", "tau", "gamma"))
  print(i)
}
