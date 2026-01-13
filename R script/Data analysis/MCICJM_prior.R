# Load packages -----------------------------------------------------------
library(splines)
library(GLMMadaptive)
library(rjags)
library(mcmcse)
library(future)
source("R script/Functions/MCICJM fit_prior.R")
load("Data/biopsy_sub.RData")
load("Data/pass_id.RData")
load("Data/pass.RData")
load("Output/Data analysis/inits/mixed model_t.RData")

plan(multisession, workers = 3)

# Modelling ---------------------------------------------------------------
# initial values
mcicjm.model <- mcicjm.fit(data = pass, biopsy = biopsy, n.adapt = 10000,
                           n.burnin = 10000, n.iter = 10000, seed = 2023,
                           prior = "unif", sens_a = 0.5, sens_b = 0.9,
                           qp = 7, df_pspline = 5,
                           GLMM_fit = mm_inits)
save(mcicjm.model, file = "Output/Data analysis/MCICJM_unif5090.RData")

## Summary
load("Output/Data analysis/MCICJM_unif5090.RData")
mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))

# HRs
# gammas
round(exp(colMeans(mcmc))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("gamma", colnames(mcmc))], 2)

# alphas, column is risk
round(exp(colMeans(mcmc))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("alpha", colnames(mcmc))], 2)
