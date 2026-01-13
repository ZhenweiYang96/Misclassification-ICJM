# Load packages -----------------------------------------------------------
library(rjags)
library(tidyverse)
library(splines)


# Table 1: HRs in data analysis -------------------------------------------
### Fixed sensitivity of 0.6
load("Output/Data analysis/MCICJM_0.6.RData")
mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
# gammas
round(exp(colMeans(mcmc))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("gamma", colnames(mcmc))], 2)

# alphas, column is risk
round(exp(colMeans(mcmc))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("alpha", colnames(mcmc))], 2)

### Fixed sensitivity of 0.8
load("Output/Data analysis/MCICJM_0.8.RData")
mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
# gammas
round(exp(colMeans(mcmc))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("gamma", colnames(mcmc))], 2)

# alphas, column is risk
round(exp(colMeans(mcmc))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("alpha", colnames(mcmc))], 2)

### Fixed sensitivity of 1
load("Output/Data analysis/MCICJM_1.RData")
mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
# gammas
round(exp(colMeans(mcmc))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("gamma", colnames(mcmc))], 2)

# alphas, column is risk
round(exp(colMeans(mcmc))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("alpha", colnames(mcmc))], 2)

### Sensitivity prior Unif(0.5, 0.9)
load("Output/Data analysis/MCICJM_unif5090.RData")
mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
# gammas
round(exp(colMeans(mcmc))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("gamma", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("gamma", colnames(mcmc))], 2)

# alphas, column is risk
round(exp(colMeans(mcmc))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.025)}))[grep("alpha", colnames(mcmc))], 2)
round(exp(apply(mcmc, 2, function(x) {quantile(x, prob=0.975)}))[grep("alpha", colnames(mcmc))], 2)


# Table 2: MSE of parameters ----------------------------------------------
load("Output/Plots/summary storage/gamma_pr.RData")
load("Output/Plots/summary storage/alpha_pr.RData")

# Gammas
pr.gamma_sim %>% 
  mutate(sqerror = (est - true)^2) %>% 
  select(-type) %>% 
  group_by(gamma, Model) %>% 
  summarise(mean = round(mean(sqerror, na.rm = T),3))

# Alphas
pr.alpha_sim %>% 
  mutate(sqerror = (est - true)^2) %>% 
  select(-type) %>% 
  group_by(alpha, Model) %>% 
  summarise(mean = round(mean(sqerror, na.rm = T),3))


# Table 3: MSE of baseline hazards ----------------------------------------
# load("Output/Plots/summary storage/bh_pr.RData")
# 
# # True hazard
# time_vec <- 1:6
# len <- length(time_vec)
# knts_true <- c(0, 0, 0, 0, 1.386667, 2.773333, 4.160000, 5.546667, 6.933333, 8.320000, 9.706667, 11.093333, 12.48, 12.48, 12.48, 12.48)
# gambh_true <- matrix(c(-3.022865564, -2.574942908, -2.165002586, -1.868394333, -1.780851418, -1.866034743, -2.043196689, -2.260330582, -2.506676980, -2.740856927, -2.954514094, -3.146719267,
#                        -5.128957296, -4.551501961, -4.263715050, -4.308050583, -4.474338093, -4.645207487, -4.802948164, -4.912553004, -5.113348423, -5.372396769, -5.660759223, -5.967330241), 12, 2)
# 
# bh_true <- data.frame(
#   time = rep(time_vec, 2),
#   true_hazard = exp(c(splineDesign(knts_true, time_vec) %*% 
#                    gambh_true)),
#   event = rep(c("Progression", "Treatment"), each = len)
# )
# 
# bh <- merge(bh_sim %>% filter(time %in% time_vec), bh_true)
# 
# summ <- bh %>% 
#   mutate(sqerror = (hazard - true_hazard)^2) %>% 
#   group_by(event, time, model) %>% 
#   summarise(mean = round(mean(sqerror, na.rm = T),4)) %>% 
#   filter(model %in% c("Sens = 100%", "Sens = 75%"))
