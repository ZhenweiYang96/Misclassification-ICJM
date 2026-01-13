load("Output/Data analysis/MCICJM_0.75.RData")
library(ggplot2)
library(tidyverse)
library(splines)
library(ggpubr)
library(ggh4x)

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

# Figure 1: Parameter recovery and uncertainty ----------------------------

## Alpha -------------------------------------------------------------------
## Association parameters
# ESTIMATES
factor_order <- factor(c("Prg: PSA value", "Prg: PSA yearly change", 
                         "Trt: PSA value", "Trt: PSA yearly change"),
                       levels = c("Prg: PSA value", "Prg: PSA yearly change", 
                                  "Trt: PSA value", "Trt: PSA yearly change"),
                       ordered = T)
pr.alpha_true <- data.frame(alpha = 1:4,
                            value = c(0.162897891, 1.788904469, 0.399214083, 2.222056696))
pr.alpha_sim <- data.frame(dat_id = rep(rep(1:200, each = 4), length(model_order)),
                           alpha = rep(rep(factor_order, 200), length(model_order)),
                           est = NA,
                           true = rep(rep(c(0.162897891, 1.788904469, 0.399214083, 2.222056696), 
                                          200), length(model_order)),
                           Model = rep(model_order, each = 200*4),
                           type = "Fixed effects")


for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  idx_alpha <- grep("alpha", names(mcicjm$summary$coef))
  pr.alpha_sim[pr.alpha_sim$Model=="Sens = 75%" & 
                 pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    idx_alpha <- grep("alpha", names(mcicjm$summary$coef))
    pr.alpha_sim[pr.alpha_sim$Model=="Sens = 25%" & 
                   pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    pr.alpha_sim[pr.alpha_sim$Model=="Sens = 60%" &
                   pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  }
  
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  pr.alpha_sim[pr.alpha_sim$Model=="Sens = 100%" & 
                 pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  idx_alpha <- grep("alpha", names(mcicjm$summary$coef))
  pr.alpha_sim[pr.alpha_sim$Model=="Sens ~ Unif(0.6, 0.9)" & 
                 pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  idx_alpha <- grep("alpha", names(mcicjm$summary$coef))
  pr.alpha_sim[pr.alpha_sim$Model=="Sens ~ Unif(0.5, 0.8)" & 
                 pr.alpha_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_alpha]
  
}

save(pr.alpha_sim, file="Output/Plots/summary storage/alpha_pr.RData")
load("Output/Plots/summary storage/alpha_pr.RData")

# RELATIVE BIAS
pr.alpha_sim <- pr.alpha_sim %>% 
  mutate(rbias = (est - true)/true)

## check
# pr.alpha_sim %>% filter(alpha %in% c("Prg: PSA value", "Prg: PSA yearly change")) %>%
#   group_by(Model, alpha) %>%
#   summarise(median = median(rbias, na.rm = T))

plot_rbias_alpha <- pr.alpha_sim %>% 
  mutate(fill = case_when(Model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(y = rbias, group = Model, 
             x = Model, fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar") +  
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.3, linewidth= 1.1) +
  facet_wrap(~ alpha, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  # scale_color_viridis_d() +
  ylab("Relative bias (estimate/true - 1)") +
  theme_bw() +
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) 
plot_rbias_alpha
# ggsave(plot_rbias_alpha, file="Output/Plots/Main text/Figure1a_sim_rbias_alpha.png", height = 8, width=13)

## Gamma -------------------------------------------------------------------
## Coefficients for baseline covariate (PSA density)
# ESTIMATES
pr.gamma_true <- data.frame(gamma = 1:2,
                            value = c(0.414315298, 0.251801149))
pr.gamma_sim <- data.frame(dat_id = rep(rep(1:200, each = 2), length(model_order)),
                           gamma = rep(rep(c("Prg: PSA density", "Trt: PSA density"), 200), length(model_order)),
                           est = NA,
                           true = rep(rep(c(0.414315298, 0.251801149), 200), length(model_order)),
                           Model = rep(model_order, each = 200*2),
                           type = "Fixed effects")

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  idx_gamma <- grep("gamma", names(mcicjm$summary$coef))
  pr.gamma_sim[pr.gamma_sim$Model=="Sens = 75%" & 
                 pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    idx_gamma <- grep("gamma", names(mcicjm$summary$coef))
    pr.gamma_sim[pr.gamma_sim$Model=="Sens = 25%" & 
                   pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    pr.gamma_sim[pr.gamma_sim$Model=="Sens = 60%" & 
                   pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
  }
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  pr.gamma_sim[pr.gamma_sim$Model=="Sens = 100%" & 
                 pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  idx_gamma <- grep("gamma", names(mcicjm$summary$coef))
  pr.gamma_sim[pr.gamma_sim$Model=="Sens ~ Unif(0.6, 0.9)" & 
                 pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  idx_gamma <- grep("gamma", names(mcicjm$summary$coef))
  pr.gamma_sim[pr.gamma_sim$Model=="Sens ~ Unif(0.5, 0.8)" & 
                 pr.gamma_sim$dat_id == i, "est"] = mcicjm$summary$coef[idx_gamma]
}

save(pr.gamma_sim, file="Output/Plots/summary storage/gamma_pr.RData")
load("Output/Plots/summary storage/gamma_pr.RData")

# RELATIVE BIAS
pr.gamma_sim <- pr.gamma_sim %>% 
  mutate(rbias = (est - true)/true)

plot_rbias_gamma <- pr.gamma_sim %>% 
  mutate(fill = case_when(Model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(y = rbias, group = Model, 
             x = Model, fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.shape = NA) +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.3, size= 1.1) +
  facet_wrap(~ gamma, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  # scale_color_viridis_d() +
  ylab("Relative bias (estimate/true - 1)") +
  theme_bw() +
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) 
plot_rbias_gamma
# ggsave(plot_rbias_gamma, file="Output/Plots/Main text/Figure1a_sim_rbias_gamma.png",
#        height = 8, width = 13)

## Alpha and Gamma ---------------------------------------------------------
plot_rbias_ap <- ggarrange(plot_rbias_alpha, plot_rbias_gamma)+
  theme(plot.margin = margin(0.1,0.1,0.1,1, "cm")) 
plot_rbias_ap
ggsave(plot_rbias_ap, file="Output/Plots/Main text/Figure1a_sim_rbias_alphagamma.png", width=16, height = 8)

# Figure 1b: Width of CI -------------------------------------------------------------

## Alpha -------------------------------------------------------------------
## Association parameters
factor_order <- factor(c("Prg: PSA value", "Prg: PSA yearly change", 
                         "Trt: PSA value", "Trt: PSA yearly change"),
                       levels = c("Prg: PSA value", "Prg: PSA yearly change", 
                                  "Trt: PSA value", "Trt: PSA yearly change"),
                       ordered = T)
ciw.alpha_sim <- data.frame(dat_id = rep(rep(1:200, each = 4), length(model_order)),
                            alpha = rep(rep(factor_order, 200), length(model_order)),
                            est = NA,
                            Model = rep(model_order, each = 200*4),
                            type = "Fixed effects")

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_alpha <- grep("alpha", names(ciw))
  ciw.alpha_sim[ciw.alpha_sim$Model=="Sens = 75%" & 
                  ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
    idx_alpha <- grep("alpha", names(ciw))
    ciw.alpha_sim[ciw.alpha_sim$Model=="Sens = 25%" & 
                    ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
    ciw.alpha_sim[ciw.alpha_sim$Model=="Sens = 60%" & 
                    ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
  }
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  ciw.alpha_sim[ciw.alpha_sim$Model=="Sens = 100%" & 
                  ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_alpha <- grep("alpha", names(ciw))
  ciw.alpha_sim[ciw.alpha_sim$Model=="Sens ~ Unif(0.6, 0.9)" & 
                  ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_alpha <- grep("alpha", names(ciw))
  ciw.alpha_sim[ciw.alpha_sim$Model=="Sens ~ Unif(0.5, 0.8)" & 
                  ciw.alpha_sim$dat_id == i, "est"] = ciw[idx_alpha]
}

save(ciw.alpha_sim, file="Output/Plots/summary storage/alpha_ciw.RData")
load("Output/Plots/summary storage/alpha_ciw.RData")

## Remove outliers
ciw.alpha_sim_no_outliers <- ciw.alpha_sim %>%
  group_by(alpha, Model) %>%
  mutate(
    q1 = quantile(est, 0.25, na.rm = TRUE),
    q3 = quantile(est, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  ) %>%
  filter(est >= lower & est <= upper) %>%
  ungroup()
##

plot_ciw_alpha <- ciw.alpha_sim_no_outliers %>% 
  mutate(fill = case_when(Model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(y = est, group = Model, x = Model,
             fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ alpha, scales = "free_y") + 
  # scale_color_viridis_d() +
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  ylab("Width of CI") + 
  theme_bw() +
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) 
plot_ciw_alpha
# ggsave(plot_ciw_alpha, file="Output/Plots/Main text/Figure1b_sim_ciw_alpha.png", 
#        width = 13, height = 8)

#### NUMBERS
summary(ciw.alpha_sim)
ciw.alpha_sim %>% 
  filter(alpha=="Prg: PSA value" & Model %in% c("Sens = 100%", "Sens = 75%")) %>%
  select(-type) %>% 
  spread(Model, est) %>% 
  mutate(prop = (`Sens = 100%`/`Sens = 75%` - 1)*100,
         diff= `Sens = 100%` - `Sens = 75%`) %>% 
  select(prop, diff) %>% 
  colMeans() %>% 
  round(2)

ciw.alpha_sim %>% 
  filter(alpha=="Prg: PSA value" & Model %in% c("Sens = 25%", "Sens = 75%")) %>%
  select(-type) %>% 
  spread(Model, est) %>% 
  mutate(prop = (`Sens = 25%`/`Sens = 75%` - 1)*100,
         diff= `Sens = 25%` - `Sens = 75%`) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

ciw.alpha_sim %>% 
  filter(alpha=="Prg: PSA value" & Model %in% c("Sens = 60%", "Sens = 75%")) %>%
  select(-type) %>% 
  spread(Model, est) %>% 
  mutate(prop = (`Sens = 60%`/`Sens = 75%` - 1)*100,
         diff= `Sens = 60%` - `Sens = 75%`) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

ciw.alpha_sim %>% 
  filter(alpha=="Prg: PSA value" & Model %in% c("Sens ~ Unif(0.6, 0.9)", "Sens = 75%")) %>%
  select(-type) %>% 
  spread(Model, est) %>% 
  mutate(prop = (`Sens ~ Unif(0.6, 0.9)`/`Sens = 75%` - 1)*100,
         diff= `Sens ~ Unif(0.6, 0.9)` - `Sens = 75%`) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

ciw.alpha_sim %>% 
  filter(alpha=="Prg: PSA value" & Model %in% c("Sens ~ Unif(0.5, 0.8)", "Sens = 75%")) %>%
  select(-type) %>% 
  spread(Model, est) %>% 
  mutate(prop = (`Sens ~ Unif(0.5, 0.8)`/`Sens = 75%` - 1)*100,
         diff= `Sens ~ Unif(0.5, 0.8)` - `Sens = 75%`) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

## Gamma -------------------------------------------------------------------
## Coefficients for baseline covariate (PSA density)
ciw.gamma_sim <- data.frame(dat_id = rep(rep(1:200, each = 2), length(model_order)),
                            gamma = rep(rep(c("Prg: PSA density", "Trt: PSA density"), 200), length(model_order)),
                            est = NA,
                            Model = rep(model_order, each = 200*2),
                            type = "Fixed effects")

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_gamma <- grep("gamma", names(ciw))
  ciw.gamma_sim[ciw.gamma_sim$Model=="Sens = 75%" & 
                  ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
    idx_gamma <- grep("gamma", names(ciw))
    ciw.gamma_sim[ciw.gamma_sim$Model=="Sens = 25%" & 
                    ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
    ciw.gamma_sim[ciw.gamma_sim$Model=="Sens = 60%" & 
                    ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
  }
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  ciw.gamma_sim[ciw.gamma_sim$Model=="Sens = 100%" & 
                  ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_gamma <- grep("gamma", names(ciw))
  ciw.gamma_sim[ciw.gamma_sim$Model=="Sens ~ Unif(0.6, 0.9)" & 
                  ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  ciw <- mcicjm$summary$coef_upper - mcicjm$summary$coef_lower
  idx_gamma <- grep("gamma", names(ciw))
  ciw.gamma_sim[ciw.gamma_sim$Model=="Sens ~ Unif(0.5, 0.8)" & 
                  ciw.gamma_sim$dat_id == i, "est"] = ciw[idx_gamma]
}

save(ciw.gamma_sim, file="Output/Plots/summary storage/gamma_ciw.RData")
load("Output/Plots/summary storage/gamma_ciw.RData")

## Remove outliers
ciw.gamma_sim_no_outliers <- ciw.gamma_sim %>%
  group_by(gamma, Model) %>%
  mutate(
    q1 = quantile(est, 0.25, na.rm = TRUE),
    q3 = quantile(est, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  ) %>%
  filter(est >= lower & est <= upper) %>%
  ungroup()
##

plot_ciw_gamma <- ciw.gamma_sim_no_outliers %>% 
  mutate(fill = case_when(Model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(y = est, group = Model, x = Model,
             fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar")+
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ gamma, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  ylab("Width of CI") + 
  theme_bw() +
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) 
plot_ciw_gamma
# ggsave(plot_ciw_gamma, file="Output/Plots/Main text/Figure1b_sim_ciw_gamma.png",
#        width = 13, height = 8)


## Alpha and Gamma ---------------------------------------------------------
plot_ciw_ap <- ggarrange(plot_ciw_alpha, plot_ciw_gamma)
plot_ciw_ap
ggsave(plot_ciw_ap, file="Output/Plots/Main text/Figure1b_sim_ciw_alphagamma.png", 
       width = 16, height = 8)


# Figure 2: baseline hazard -----------------------------------------------
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

bh_sim <- data.frame(
  dat_id = rep(rep(rep(1:200, each = len), length(model_order)), 2),
  time = rep(time_vec, length(model_order)*200*2),
  hazard = NA,
  hazard_lower = NA,
  hazard_upper = NA,
  model = rep(rep(model_order, each = len * 200), 2),
  event = rep(c("Progression", "Treatment"), each = 200*len*length(model_order))
)

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  knts <- mcicjm$model_info$knots$knot.surv
  idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
  # coef <- matrix(mcicjm$summary$coef[idx_gambh], ncol = 2)
  spline_mat <- splineDesign(knts, time_vec)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens = 75%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    knts <- mcicjm$model_info$knots$knot.surv
    idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
    # coef <- matrix(mcicjm$summary$coef[idx_gambh], ncol = 2)
    spline_mat <- splineDesign(knts, time_vec)
    coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
    
    bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
    
    bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
    bh_sim[bh_sim$model=="Sens = 25%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  }
  
  # Sens = 60%
  if (!(i %in% c(137,186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    knts <- mcicjm$model_info$knots$knot.surv
    idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
    spline_mat <- splineDesign(knts, time_vec)
    coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
    
    bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
    
    bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
    bh_sim[bh_sim$model=="Sens = 60%" & 
             bh_sim$dat_id == i & 
             bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  }
  
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  knts <- mcicjm$model_info$knots$knot.surv
  idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
  spline_mat <- splineDesign(knts, time_vec)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens = 100%" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  knts <- mcicjm$model_info$knots$knot.surv
  idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
  spline_mat <- splineDesign(knts, time_vec)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens ~ Unif(0.6, 0.9)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  knts <- mcicjm$model_info$knots$knot.surv
  idx_gambh <- grep("gambh", names(mcicjm$summary$coef))
  spline_mat <- splineDesign(knts, time_vec)
  coef.mat <- do.call("rbind", lapply(mcicjm$mcmc, "[[", 1))[,idx_gambh]
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,1:(ncol(coef.mat)/2)]))
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Progression", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
  
  bh_mat <- exp(splineDesign(knts, time_vec, outer.ok = T) %*% t(coef.mat[,-(1:(ncol(coef.mat)/2))]))
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard"] = rowMeans(bh_mat)
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_lower"] = apply(bh_mat, 1, function(x) {quantile(x, 0.025)})
  bh_sim[bh_sim$model=="Sens ~ Unif(0.5, 0.8)" & 
           bh_sim$dat_id == i & 
           bh_sim$event == "Treatment", "hazard_upper"] = apply(bh_mat, 1, function(x) {quantile(x, 0.975)})
}

save(bh_sim, file="Output/Plots/summary storage/bh_pr.RData")
load("Output/Plots/summary storage/bh_pr.RData")
## Line plot
# plot_bh_logest <- bh_sim %>% 
#   ggplot() + 
#   geom_line(aes(x = time, y = log(hazard), 
#                 group = dat_id), alpha = 0.1) +
#   geom_line(data = bh_true,
#             aes(x = time, y = log(hazard)), size = 1.2) +
#   facet_wrap(~ event * model, scales = "free_y") + 
#   scale_color_viridis_d() +
#   xlab("Time (years)") + 
#   ylab("log(baseline hazard)") + 
#   theme_bw() + 
#   theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
#         axis.text.x = element_text(size = 16, face = "bold", family = "LM Roman 10"),
#         axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10"))
# plot_bh_logest
# ggsave(plot_bh_logest, file="Output/Plots/Main text/Figure2a_sim_bh_logest.png", 
#        width=15, height=8)

### Version 1: Connected boxplot of baseline hazard relative bias at 6 time points
rbias_bh_sim <- bh_sim %>% 
  mutate(true_hazard = c(rep(bh_true$hazard[bh_true$event == "Progression"], length(model_order)*200),
                         rep(bh_true$hazard[bh_true$event == "Treatment"], length(model_order)*200))) %>% 
  filter(time %in% 1:6) %>% 
  mutate(rbias = (hazard - true_hazard)/true_hazard,
         time_label = factor(time))

rbias_bh_sim$model <- factor(
  rbias_bh_sim$model,
  levels = c("Sens = 100%", "Sens = 75%", "Sens = 60%", 
             "Sens ~ Unif(0.6, 0.9)", "Sens ~ Unif(0.5, 0.8)", 
             "Sens = 25%")
)

# Get median
median_rbias_bh_sim <- rbias_bh_sim %>% 
  group_by(model, time_label, event) %>% 
  summarize(median = median(rbias, na.rm = T))

## Facet level 
limits_list <- rbias_bh_sim %>%
  filter(time_label == 6) %>% 
  # Ensure columns exist and are not NA
  filter(!is.na(event), !is.na(time_label), !is.na(rbias)) %>%
  group_by(event, model) %>%
  summarise(
    lower = -1,
    upper = quantile(rbias, 0.75, na.rm = TRUE) + 2 * IQR(rbias, na.rm = TRUE),
    .groups = "drop"
  ) 

# Connected boxplots
plot_bh_rbias_pts <- rbias_bh_sim %>% 
  mutate(fill = case_when(model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(x = time_label, y = rbias)) +
  geom_line(data = median_rbias_bh_sim, aes(x = time_label, y = median, group = interaction(model, event)), linewidth = 1.5) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot(aes(group = interaction(time_label, model, event), fill = as.factor(fill)), outlier.shape = NA) +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.3, linewidth= 1.1) +
  facet_wrap(~ event * model, nrow = 2, ncol = 6, scales = "free_y") +
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  xlab("Time points (year)") + 
  ylab("Relative bias in baseline hazards (estimate/true - 1)") + 
  ylim(-1,3) +
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) +
  facetted_pos_scales(
    y = lapply(split(limits_list, seq(nrow(limits_list))), function(lim) {
      scale_y_continuous(limits = c(lim$lower, lim$upper))
    })
  )
plot_bh_rbias_pts
ggsave(plot_bh_rbias_pts, file="Output/Plots/Main text/Figure2a_sim_bh_rbias_points_connected.png", 
       width=18, height=10)

### Version 2: per event and per time points
rbias_bh_sim <- bh_sim %>% 
  mutate(true_hazard = c(rep(bh_true$hazard[bh_true$event == "Progression"], 5*200),
                         rep(bh_true$hazard[bh_true$event == "Treatment"], 5*200))) %>% 
  filter(time %in% c(2, 4, 6)) %>% 
  mutate(rbias = (hazard - true_hazard)/true_hazard,
         time_label = paste("Year", time))

# Create a list of limits for each facet
limits_list <- rbias_bh_sim %>%
  filter(time_label == 6) %>% 
  # Ensure columns exist and are not NA
  filter(!is.na(event), !is.na(time_label), !is.na(rbias)) %>%
  group_by(event, time_label) %>%
  summarise(
    lower = quantile(rbias, 0.75, na.rm = TRUE) - 2 * IQR(rbias, na.rm = TRUE),
    upper = quantile(rbias, 0.75, na.rm = TRUE) + 3 * IQR(rbias, na.rm = TRUE),
    .groups = "drop"
  )

plot_bh_rbias_pts <- rbias_bh_sim %>% 
  mutate(fill = case_when(model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(x = model, y = rbias, group = model, 
             fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  geom_hline(aes(yintercept = 0), 
             alpha = 0.3, linewidth= 1.1) +
  facet_wrap(~ event * time_label, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  xlab("Models") + 
  ylab("Relative bias (estimate/true - 1)") +
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) +
  facetted_pos_scales(
    y = lapply(split(limits_list, seq(nrow(limits_list))), function(lim) {
      scale_y_continuous(limits = c(lim$lower, lim$upper))
    })
  )
plot_bh_rbias_pts
ggsave(plot_bh_rbias_pts, file="Output/Plots/Main text/Figure2a_sim_bh_rbias_points.png", 
       width=15, height=8)


##### NUMBERS
bh_sim %>% 
  filter(time == 2 & event == "Progression" & model %in% c("Sens = 100%", "Sens = 75%")) %>%
  mutate(hazard = log(hazard)) %>% 
  select(-hazard_lower, -hazard_upper) %>% 
  spread(key = model, value = hazard) %>%
  mutate(prop = -(`Sens = 100%`/`Sens = 75%` - 1)*100,
         diff= `Sens = 100%` - `Sens = 75%`) %>%
  # mutate(prop = (log(`Sens = 100%`)/log(`Sens = 75%`) - 1)*100,
  #        diff= log(`Sens = 100%`) - log(`Sens = 75%`)) %>%
  select(prop, diff) %>% 
  colMeans() %>% 
  round(2)

bh_sim %>% 
  filter(time == 2 & event == "Progression" & model %in% c("Sens = 60%", "Sens = 75%")) %>%
  select(-hazard_lower, -hazard_upper) %>% 
  spread(model, hazard) %>%
  mutate(prop = -(log(`Sens = 60%`)/log(`Sens = 75%`) - 1)*100,
         diff= log(`Sens = 60%`) - log(`Sens = 75%`)) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

bh_sim %>% 
  filter(time == 2 & event == "Progression" & model %in% c("Sens = 25%", "Sens = 75%")) %>%
  select(-hazard_lower, -hazard_upper) %>% 
  spread(model, hazard) %>%
  mutate(prop = -(log(`Sens = 25%`)/log(`Sens = 75%`) - 1)*100,
         diff= log(`Sens = 25%`) - log(`Sens = 75%`)) %>% 
  select(prop, diff) %>% 
  colMeans(na.rm = T) %>% 
  round(2)

bh_sim %>% 
  filter(time == 2 & event == "Progression" & model %in% c("Sens ~ Unif(0.6, 0.9)", "Sens = 75%")) %>%
  select(-hazard_lower, -hazard_upper) %>% 
  spread(model, hazard) %>%
  mutate(prop = -(log(`Sens ~ Unif(0.6, 0.9)`)/log(`Sens = 75%`) - 1)*100,
         diff= log(`Sens ~ Unif(0.6, 0.9)`) - log(`Sens = 75%`)) %>% 
  select(prop, diff) %>% 
  colMeans() %>% 
  round(2)

bh_sim %>% 
  filter(time == 2 & event == "Progression" & model %in% c("Sens ~ Unif(0.5, 0.8)", "Sens = 75%")) %>%
  select(-hazard_lower, -hazard_upper) %>% 
  spread(model, hazard) %>%
  mutate(prop = -(log(`Sens ~ Unif(0.5, 0.8)`)/log(`Sens = 75%`) - 1)*100,
         diff= log(`Sens ~ Unif(0.5, 0.8)`) - log(`Sens = 75%`)) %>% 
  select(prop, diff) %>% 
  colMeans() %>% 
  round(2)

#################################################################################

## Version 1: Width of baseline hazard at 6 time points
ciw_bh_sim <- bh_sim %>% 
  filter(time %in% 1:6) %>% 
  mutate(ciw = hazard_upper - hazard_lower,
         time_label = factor(time))

ciw_bh_sim$model <- factor(
  ciw_bh_sim$model,
  levels = c("Sens = 100%", "Sens = 75%", "Sens = 60%", 
             "Sens ~ Unif(0.6, 0.9)", "Sens ~ Unif(0.5, 0.8)", 
             "Sens = 25%")
)


# Get median
median_ciw_bh_sim <- ciw_bh_sim %>% 
  group_by(model, time_label, event) %>% 
  summarize(median = median(ciw, na.rm = T))

# Create a list of limits for each facet
limits_list <- ciw_bh_sim %>%
  filter(time_label == 6) %>% 
  # Ensure columns exist and are not NA
  filter(!is.na(event), !is.na(time_label), !is.na(ciw)) %>%
  group_by(event, model) %>%
  summarise(
    lower = 0,
    upper = quantile(ciw, 0.75, na.rm = TRUE) + 2.4* IQR(ciw, na.rm = TRUE),
    .groups = "drop"
  )

plot_bh_ciw_pts_connected <- ciw_bh_sim %>% 
  mutate(fill = case_when(model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(x = time_label, y = ciw)) +
  geom_line(data = median_ciw_bh_sim , aes(x = time_label, y = median, group = model), linewidth = 1.5) +
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot(aes(group = interaction(time_label, model), fill = as.factor(fill)), outlier.shape = NA) +
  facet_wrap(~ event * model, nrow = 2, ncol = 6, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  xlab("Time points (year)") + 
  ylab("CI width in baseline hazard") +
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10")) +
  facetted_pos_scales(
    y = lapply(split(limits_list, seq(nrow(limits_list))), function(lim) {
      scale_y_continuous(limits = c(lim$lower, lim$upper))
    })
  )
plot_bh_ciw_pts_connected
ggsave(plot_bh_ciw_pts_connected, file="Output/Plots/Main text/Figure2b_sim_bh_unc_points_connected.png", 
       width=18, height=10)

## Version 2: Width of baseline hazard at time 2, 4, 6
ciw_bh_sim <- bh_sim %>% 
  filter(time %in% c(2, 4, 6)) %>% 
  mutate(ciw = hazard_upper - hazard_lower,
         time_label = paste("Year", time))

## Remove outliers
ciw_bh_sim_no_outliers <- ciw_bh_sim %>%
  group_by(time_label, model) %>%
  mutate(
    q1 = quantile(ciw, 0.25, na.rm = TRUE),
    q3 = quantile(ciw, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    lower = q1 - 1.5 * iqr,
    upper = q3 + 1.5 * iqr
  ) %>%
  filter(ciw >= lower & ciw <= upper) %>%
  ungroup()
##

plot_bh_unc_pts <- ciw_bh_sim_no_outliers %>% 
  mutate(fill = case_when(model=="Sens = 75%" ~ 1, T~0)) %>% 
  ggplot(aes(x = model, y = ciw, group = model, 
             fill = as.factor(fill))) + 
  stat_boxplot(geom = "errorbar") + 
  geom_boxplot(outlier.shape = NA, coef = 1.5) +
  facet_wrap(~ event * time_label, scales = "free_y") + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  xlab("Models") + 
  ylab("CI width in baseline hazard") +
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 25, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10"))

plot_bh_unc_pts
ggsave(plot_bh_unc_pts, file="Output/Plots/Main text/Figure2b_sim_bh_unc_points.png", 
       width=15, height=8)



# Figure 3: coverage of coefficients --------------------------------------
true_parameters <- mcicjm.model$summary$coef
true_parameters <- true_parameters[-c(grep("gambh", names(true_parameters)), 
                                      grep("betaL", names(true_parameters)), 
                                      grep("D", names(true_parameters)))]

cov_sim_long <- data.frame(
  dat_id = rep(1:200, each = length(model_order)*length(true_parameters)),
  Model = rep(rep(model_order, length(true_parameters)), 200),
  parameter = rep(rep(names(true_parameters), each = length(model_order)), 200),
  cover = 0
)

for (i in 1:200) {
  # Sens = 75%
  load(paste0("Output/Simulation/MCICJM/Fixed_75/MCICJM_", i, ".RData"))
  param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                    grep("D", names(mcicjm$summary$coef)),
                    grep("betaL", names(mcicjm$summary$coef)))
  upper <- mcicjm$summary$coef_upper[-param_no_idx]
  lower <- mcicjm$summary$coef_lower[-param_no_idx]
  cov_sim_long[cov_sim_long$Model=="Sens = 75%" & 
                 cov_sim_long$dat_id == i, "cover"] <- 
    upper >= true_parameters & 
    lower <= true_parameters
  
  # Sens = 25%
  if (!(i %in% c(25,40,42,47,67,81,109,111,128,134,136,159,187,189))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_25/MCICJM_", i, ".RData"))
    param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                      grep("D", names(mcicjm$summary$coef)),
                      grep("betaL", names(mcicjm$summary$coef)))
    upper <- mcicjm$summary$coef_upper[-param_no_idx]
    lower <- mcicjm$summary$coef_lower[-param_no_idx]
    cov_sim_long[cov_sim_long$Model=="Sens = 25%" & 
                   cov_sim_long$dat_id == i, "cover"] <- 
      upper >= true_parameters & 
      lower <= true_parameters
  }
  
  # Sens = 60%
  if (!(i %in% c(137, 186))) {
    load(paste0("Output/Simulation/MCICJM/Fixed_60/MCICJM_", i, ".RData"))
    param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                      grep("D", names(mcicjm$summary$coef)),
                      grep("betaL", names(mcicjm$summary$coef)))
    upper <- mcicjm$summary$coef_upper[-param_no_idx]
    lower <- mcicjm$summary$coef_lower[-param_no_idx]
    cov_sim_long[cov_sim_long$Model=="Sens = 60%" & 
                   cov_sim_long$dat_id == i, "cover"] <- 
      upper >= true_parameters & 
      lower <= true_parameters
  }
  
  # Sens = 100%
  load(paste0("Output/Simulation/MCICJM/Fixed_100/MCICJM_", i, ".RData"))
  param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                    grep("D", names(mcicjm$summary$coef)),
                    grep("betaL", names(mcicjm$summary$coef)))
  upper <- mcicjm$summary$coef_upper[-param_no_idx]
  lower <- mcicjm$summary$coef_lower[-param_no_idx]
  cov_sim_long[cov_sim_long$Model=="Sens = 100%" & 
                 cov_sim_long$dat_id == i, "cover"] <- 
    upper >= true_parameters & 
    lower <= true_parameters
  
  # Sens ~ unif(0.6, 0.9)
  load(paste0("Output/Simulation/MCICJM/Unifprior_60_90/MCICJM_", i, ".RData"))
  param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                    grep("D", names(mcicjm$summary$coef)),
                    grep("betaL", names(mcicjm$summary$coef)),
                    grep("sens", names(mcicjm$summary$coef)))
  upper <- mcicjm$summary$coef_upper[-param_no_idx]
  lower <- mcicjm$summary$coef_lower[-param_no_idx]
  cov_sim_long[cov_sim_long$Model=="Sens ~ Unif(0.6, 0.9)" & 
                 cov_sim_long$dat_id == i, "cover"] <- 
    upper >= true_parameters & 
    lower <= true_parameters
  # Sens ~ unif(0.5, 0.8)
  load(paste0("Output/Simulation/MCICJM/Unifprior_50_80/MCICJM_", i, ".RData"))
  param_no_idx <- c(grep("gambh", names(mcicjm$summary$coef)),
                    grep("D", names(mcicjm$summary$coef)),
                    grep("betaL", names(mcicjm$summary$coef)),
                    grep("sens", names(mcicjm$summary$coef)))
  upper <- mcicjm$summary$coef_upper[-param_no_idx]
  lower <- mcicjm$summary$coef_lower[-param_no_idx]
  cov_sim_long[cov_sim_long$Model=="Sens ~ Unif(0.5, 0.8)" & 
                 cov_sim_long$dat_id == i, "cover"] <- 
    upper >= true_parameters & 
    lower <= true_parameters
}
save(cov_sim_long, file = "Output/Plots/summary storage/coverage.RData")
load("Output/Plots/summary storage/coverage.RData")

cov_sim <- cov_sim_long %>% 
  group_by(Model, parameter) %>% 
  summarise(coverage = ifelse(first(Model %in% c("Sens = 60%", "Sens = 25%")), 
                              ifelse(first(Model) == "Sens = 60%", sum(cover)/198, sum(cover)/187), 
                              sum(cover)/200)) #ifelse(Model == "Sens = 60%", sum(cover)/198, sum(cover)/200)

plot_cov <- cov_sim %>% 
  mutate(
    new_parameter = case_when(parameter == "alpha[1,1]" ~ "Prg: PSA value",
                              parameter == "alpha[2,1]" ~ "Prg: PSA yearly change",
                              parameter == "alpha[1,2]" ~ "Trt: PSA value",
                              parameter == "alpha[2,2]" ~ "Trt: PSA yearly change",
                              parameter == "gamma[1]" ~ "Prg: PSA density",
                              parameter == "gamma[2]" ~ "Trt: PSA density",
                              parameter == "tau" ~ "Residual \n standard deviation",
                              T ~ parameter),
    fill = case_when(Model=="Sens = 75%" ~ 1, T~0)) %>% 
  mutate(new_parameter = fct_relevel(new_parameter, c("Prg: PSA value", "Prg: PSA yearly change",
                                                      "Trt: PSA value", "Trt: PSA yearly change",
                                                      "Prg: PSA density", "Trt: PSA density",
                                                      "Residual \n standard deviation"))) %>% 
  ggplot() + 
  geom_point(aes(x = new_parameter, y = coverage,
                 shape = Model, group = Model, fill=as.factor(fill)), size = 5, position=position_dodge(width=0.3)) +
  geom_hline(yintercept = 0.95, color = "gray") + 
  # facet_wrap(~ Group, scales = "free_y") + 
  scale_shape(solid = FALSE) +
  xlab("Important Parameters") + 
  ylab("Coverage (%)") + 
  scale_shape_manual(values = 20:25) + 
  scale_fill_manual(values = c("white", "lightgray"),
                    guide = "none") + 
  theme_bw() + 
  theme(text = element_text(size = 22, face = "bold", family = "LM Roman 10"),
        axis.text.x = element_text(angle = 30, hjust = 1, 
                                   size = 16, face = "bold", family = "LM Roman 10"),
        axis.text.y = element_text(size = 16, face = "bold", family = "LM Roman 10"),
        legend.position = "top")
plot_cov
ggsave(plot_cov, file="Output/Plots/Main text/Figure3_sim_cov.png", width = 15, height = 8)


