load("Output/Data analysis/MCICJM_0.75.RData")
library(ggplot2)
library(tidyverse)
library(splines)
library(ggpubr)
library(ggh4x)
library(rjags)


# Figure S3: baseline hazard real data ------------------------------------
res <- data.frame(
  sens = rep(c(0.6, 0.8, 1), each = 100),
  time = rep(seq(0, 10, length.out = 100), 3),
  hazard_prg = 0,
  hazard_prg.lower = 0,
  hazard_prg.upper = 0,
  hazard_trt = 0,
  hazard_trt.lower = 0,
  hazard_trt.upper = 0
)

prg <- NULL
trt <- NULL
for (i in c(0.6, 0.8, 1)) {
  load(paste0("Output/Data analysis/MCICJM_", i, ".RData"))
  knts <- mcicjm.model$model_info$knots$knot.surv
  time.bs <- splineDesign(knts, seq(0, 10, length.out = 100), ord = 4L, outer.ok = T)
  
  mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
  mcmc <- mcmc[,grep("gambh", colnames(mcmc))]
  prg <- rbind(prg,
               exp(time.bs %*% t(mcmc[,1:12])))
  trt <- rbind(trt,
               exp(time.bs %*% t(mcmc[,13:24])))
}

res$hazard_prg <- rowMeans(prg)
res$hazard_prg.lower <- apply(prg, 1, function(x) {quantile(x, 0.025)})
res$hazard_prg.upper <- apply(prg, 1, function(x) {quantile(x, 0.975)})

res$hazard_trt <- rowMeans(trt)
res$hazard_trt.lower <- apply(trt, 1, function(x) {quantile(x, 0.025)})
res$hazard_trt.upper <- apply(trt, 1, function(x) {quantile(x, 0.975)})
save(res, file= "Output/Plots/summary storage/plotdata.RData")


# Uniform
res_unif <- data.frame(
  sens = rep("Unif(0.5, 0.9)", each = 100),
  time = seq(0, 10, length.out = 100),
  hazard_prg = 0,
  hazard_prg.lower = 0,
  hazard_prg.upper = 0,
  hazard_trt = 0,
  hazard_trt.lower = 0,
  hazard_trt.upper = 0
)
load("Output/Data analysis/MCICJM_unif5090.RData")
knts <- mcicjm.model$model_info$knots$knot.surv
time.bs <- splineDesign(knts, seq(0, 10, length.out = 100), ord = 4L, outer.ok = T)

mcmc <- do.call(rbind, lapply(mcicjm.model$mcmc, as.matrix))
mcmc <- mcmc[,grep("gambh", colnames(mcmc))]
prg <- exp(time.bs %*% t(mcmc[,1:5]))
trt <- exp(time.bs %*% t(mcmc[,6:10]))

res_unif$hazard_prg <- rowMeans(prg)
res_unif$hazard_prg.lower <- apply(prg, 1, function(x) {quantile(x, 0.025)})
res_unif$hazard_prg.upper <- apply(prg, 1, function(x) {quantile(x, 0.975)})

res_unif$hazard_trt <- rowMeans(trt)
res_unif$hazard_trt.lower <- apply(trt, 1, function(x) {quantile(x, 0.025)})
res_unif$hazard_trt.upper <- apply(trt, 1, function(x) {quantile(x, 0.975)})


load("Output/Plots/summary storage/plotdata.RData")
# Progression part
plot.prg <- ggplot(res) + 
  geom_line(aes(x = time, y = hazard_prg, color = as.factor(sens)), size = 1.5) + 
  geom_ribbon(aes(x = time, ymin = hazard_prg.lower, ymax = hazard_prg.upper, fill = as.factor(sens)), alpha = 0.2) + 
  geom_line(data = res_unif, aes(x = time, y = hazard_prg), size = 1.5, alpha = 0.5) +
  geom_line(data = res_unif, aes(x = time, y = hazard_prg.lower), size = 1.5, alpha = 0.5, linetype = "dashed") +
  geom_line(data = res_unif, aes(x = time, y = hazard_prg.upper), size = 1.5, alpha = 0.5, linetype = "dashed") +
  #facet_wrap(~ as.factor(sens), ncol = 3) + 
  ylab("Progression-specific \n baseline hazard") + 
  xlab("Time (years)") + 
  theme_bw() + 
  theme(text = element_text(face = "bold", size = 28, family = "LM Roman 10"),
        axis.text = element_text(face = "bold", size = 24, family = "LM Roman 10")) +
  scale_colour_manual(breaks = c(0.6, 0.8, 1), labels = c(0.6, 0.8, 1),
                      values = c("#F8766D", "#00BA38", "#619CFF"), name = "Sensitivity") +
  scale_fill_manual(breaks = c(0.6, 0.8, 1), labels = c(0.6, 0.8, 1),
                    values = c("#F8766D", "#00BA38", "#619CFF"), name = "Sensitivity")
plot.prg
# ggsave(plot.prg, file="Output/Plots/Supplementary/FigureS3a_baseline hazard_progression.png",
#        width = 12, height = 8)

# Treatment part
plot.trt <- ggplot(res) + 
  geom_line(aes(x = time, y = hazard_trt, color = as.factor(sens)), size = 1.5) + 
  geom_ribbon(aes(x = time, ymin = hazard_trt.lower, ymax = hazard_trt.upper, fill = as.factor(sens)), alpha = 0.2) + 
  geom_line(data = res_unif, aes(x = time, y = hazard_trt), size = 1.5, alpha = 0.5) +
  geom_line(data = res_unif, aes(x = time, y = hazard_trt.lower), size = 1.5, alpha = 0.5, linetype = "dashed") +
  geom_line(data = res_unif, aes(x = time, y = hazard_trt.upper), size = 1.5, alpha = 0.5, linetype = "dashed") +
  ylab("Treatment-specific \n baseline hazard") + 
  xlab("Time (years)") + 
  theme_bw() + 
  theme(text = element_text(face = "bold", size = 28, family = "LM Roman 10"),
        axis.text = element_text(face = "bold", size = 24, family = "LM Roman 10")) +
  scale_colour_manual(breaks = c(0.6, 0.8, 1), labels = c(0.6, 0.8, 1),
                      values = c("#F8766D", "#00BA38", "#619CFF"), name = "Sensitivity") +
  scale_fill_manual(breaks = c(0.6, 0.8, 1), labels = c(0.6, 0.8, 1),
                    values = c("#F8766D", "#00BA38", "#619CFF"), name = "Sensitivity")
plot.trt
# ggsave(plot.trt, file="Output/Plots/Supplementary/FigureS3b_baseline hazard_treatment.png",
#        width = 12, height = 8)

plot.pool <- ggarrange(plot.prg, plot.trt, ncol=2, common.legend = TRUE, legend="top")
ggsave(plot.pool, file="Output/Plots/Supplementary/FigureS3_baseline hazard_pool.png",
       width = 15, height = 8)


# Figure S4 & S5: PSA Trajectory ---------------------------------------------------------
lknts_true <- mcicjm.model$model_info$knots$knot.longi
betaL_true <- mcicjm.model$summary$coef[grepl("betaL", names(mcicjm.model$summary$coef))]
time_vec <- sort(c(seq(0, 6, length.out = 48), 2, 4))
psa_true <- data.frame(
  time = time_vec,
  psa = 2^(cbind(1, 
                 ns(time_vec, knots = lknts_true[2:3], B = lknts_true[c(1,4)]),
                 0) %*% betaL_true) - 1
)

model_order <- factor(c("Sens = 100%", 
                         "Sens = 75%", "Sens = 60%",
                         "Sens = 25%",
                         "Sens ~ Unif(0.6, 0.9)",
                         "Sens ~ Unif(0.5, 0.8)"),
                       levels = c("Sens = 100%", 
                                  "Sens = 75%", "Sens = 60%", 
                                  "Sens = 25%",
                                  "Sens ~ Unif(0.6, 0.9)",
                                  "Sens ~ Unif(0.5, 0.8)"),
                       ordered = T)

psa_sim <- data.frame(
  dat_id = rep(rep(1:200, each = 50), length(model_order)),
  time = rep(time_vec, length(model_order)*200),
  psa = NA,
  psa_lower = NA,
  psa_upper = NA,
  model = rep(model_order, each = 50 * 200)
)

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  lknts <- mcicjm$model_info$knots$knot.longi
  idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
  spline_mat <- cbind(1, 
                      ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                      0)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
  
  psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
  psa_sim[psa_sim$model=="Sens = 75%" & 
            psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
  psa_sim[psa_sim$model=="Sens = 75%" & 
            psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
  psa_sim[psa_sim$model=="Sens = 75%" & 
            psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    lknts <- mcicjm$model_info$knots$knot.longi
    idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
    spline_mat <- cbind(1, 
                        ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                        0)
    coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
    
    psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
    psa_sim[psa_sim$model=="Sens = 25%" & 
              psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
    psa_sim[psa_sim$model=="Sens = 25%" & 
              psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
    psa_sim[psa_sim$model=="Sens = 25%" & 
              psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    lknts <- mcicjm$model_info$knots$knot.longi
    idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
    spline_mat <- cbind(1, 
                        ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                        0)
    coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
    
    psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
    psa_sim[psa_sim$model=="Sens = 60%" & 
              psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
    psa_sim[psa_sim$model=="Sens = 60%" & 
              psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
    psa_sim[psa_sim$model=="Sens = 60%" & 
              psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  }
  
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  lknts <- mcicjm$model_info$knots$knot.longi
  idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
  spline_mat <- cbind(1, 
                      ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                      0)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
  
  psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
  psa_sim[psa_sim$model=="Sens = 100%" & 
            psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
  psa_sim[psa_sim$model=="Sens = 100%" & 
            psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
  psa_sim[psa_sim$model=="Sens = 100%" & 
            psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  lknts <- mcicjm$model_info$knots$knot.longi
  idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
  spline_mat <- cbind(1, 
                      ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                      0)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
  
  psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
  psa_sim[psa_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
            psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
  psa_sim[psa_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
            psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
  psa_sim[psa_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
            psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  lknts <- mcicjm$model_info$knots$knot.longi
  idx_betaL <- grep("betaL", names(mcicjm$summary$coef))
  spline_mat <- cbind(1, 
                      ns(time_vec, knots = lknts[2:3], B = lknts[c(1,4)]),
                      0)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_betaL]
  
  psa_mat <- 2^(spline_mat %*% t(coef.mat)) - 1
  psa_sim[psa_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
            psa_sim$dat_id == i, "psa"] = rowMeans(psa_mat)
  psa_sim[psa_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
            psa_sim$dat_id == i, "psa_lower"] = apply(psa_mat, 1, function(x) {quantile(x, 0.025)})
  psa_sim[psa_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
            psa_sim$dat_id == i, "psa_upper"] = apply(psa_mat, 1, function(x) {quantile(x, 0.975)})
  
}

save(psa_sim, file="Output/Plots/summary storage/psa_pr.RData")
load("Output/Plots/summary storage/psa_pr.RData")

### Figure S4
plot_psa_est <- psa_sim %>% 
  ggplot() + 
  geom_line(aes(x = time, y = psa, 
                group = dat_id, color=model), alpha = 0.3) +
  geom_line(data = psa_true,
            aes(x = time, y = psa), color = "red", linewidth = 1.2) +
  facet_wrap(~ model, scales = "free_y") + 
  scale_color_viridis_d() +
  xlab("Time (years)") + 
  ylab("PSA trajectory (population average)") + 
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        legend.position = "none")
plot_psa_est
ggsave(plot_psa_est, file="Output/Plots/Supplementary/FigureS4_sim_psa_est.png", 
       width = 10)

### Figure S5 Width of PSA trajectory at time 2, 4, 6
ciw_psa_sim <- psa_sim %>% 
  filter(time %in% c(2, 4, 6)) %>% 
  mutate(ciw = psa_upper - psa_lower,
         time_label = paste("Year", time))

plot_psa_unc_pts <- ciw_psa_sim %>% 
  mutate(fill = case_when(model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(x = model, y = ciw, group = model, color=model,
             fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar") +  
  geom_boxplot(outlier.fill = NA, outlier.shape = NA) +
  ylim(0,1.25)+
  xlab("Time (years)") + 
  scale_color_viridis_d(guide = "none") +
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  facet_wrap(~time_label)+
  ylab("CI width in PSA trajectory") +
  theme_bw() + 
  theme(text = element_text(size = 25, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 20, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 20, face = "bold", family = "LM Roman 10"),
        plot.margin = unit(c(0.1,0.1,0.1,1), "cm"),
        strip.text = element_text(size = 30, face = "bold", family = "LM Roman 10"))
plot_psa_unc_pts
ggsave(plot_psa_unc_pts, file="Output/Plots/Supplementary/FigureS5_sim_psa_unc_points.png", width=16, height = 8)


# Figure S6: estimated baseline hazard ------------------------------------
time_vec <- seq(0, 6, 0.1)
len <- length(time_vec)
knts_true <- c(0, 0, 0, 0, 1.386667, 2.773333, 4.160000, 5.546667, 6.933333, 8.320000, 9.706667, 11.093333, 12.48, 12.48, 12.48, 12.48)
gambh_true <- matrix(c(-3.022865564, -2.574942908, -2.165002586, -1.868394333, -1.780851418, -1.866034743, -2.043196689, -2.260330582, -2.506676980, -2.740856927, -2.954514094, -3.146719267,
                       -5.128957296, -4.551501961, -4.263715050, -4.308050583, -4.474338093, -4.645207487, -4.802948164, -4.912553004, -5.113348423, -5.372396769, -5.660759223, -5.967330241), 12, 2)

bh_true <- data.frame(
  time = rep(time_vec, 2),
  hazard = exp(c(splineDesign(knts_true, time_vec) %*% 
                   gambh_true)),
  event = rep(c("Progression", "Treatment"), each = len)
)
load("Output/Plots/summary storage/bh_pr.RData")
# Line plot
plot_bh_logest <- bh_sim %>%
  ggplot() +
  geom_line(aes(x = time, y = log(hazard),
                group = dat_id), alpha = 0.1) +
  geom_line(data = bh_true,
            aes(x = time, y = log(hazard)), linewidth = 1.2) +
  facet_wrap(~ event * model, scales = "free_y") +
  scale_color_viridis_d() +
  xlab("Time (years)") +
  ylab("log(baseline hazard)") +
  theme_bw() +
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10"))
plot_bh_logest
ggsave(plot_bh_logest, file="Output/Plots/Supplementary/FigureS6_sim_bh_logest.png",
       width=15, height=8)

