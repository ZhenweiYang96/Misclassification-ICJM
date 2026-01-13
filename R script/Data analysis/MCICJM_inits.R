
# Packages ----------------------------------------------------------------
library(GLMMadaptive)
library(splines)
load("data/pass.RData")

# Inits of mixed model ----------------------------------------------------
mm <- mixed_model(fixed = PSAValue ~ ns(TimeSince_Dx, 3, B = c(0.00, 12.48)) + DxAge,
                  data = pass, random = ~ ns(TimeSince_Dx, 3, B = c(0.00, 12.48)) | CISNET_ID,
                  family = students.t(df = 3), n_phis = 1,
                  initial_values = list("betas" = gaussian()))
save(mm, file = "Output/Data analysis/inits/mixed model_t.RData")

mm_inits <- list(betaL = fixef(mm),
                 b = ranef(mm),
                 inv_D = solve(mm$D),
                 tau = 1/exp(mm$phis)^2)
save(mm_inits, file = "Output/Data analysis/inits/mixed model_inits.RData")
