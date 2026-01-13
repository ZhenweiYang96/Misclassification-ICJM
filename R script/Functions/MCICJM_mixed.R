mcicjm_mixed <- function(df = 3, data,
                       mixed_path) {
  
  # Load packages -----------------------------------------------------------
  source("R script/Functions/function.R")
  
  # Start with the mixed model ----------------------------------------------
  bound <- range(data$TimeSince_Dx, data$time.cmp1, data$time.cmp2)
  data[,(ncol(data)+1):(ncol(data)+3)] <- ns(data$TimeSince_Dx, 3, B = bound)
  colnames(data)[(ncol(data)-2):ncol(data)] <- c("TimeSince_Dx.1",
                                                 "TimeSince_Dx.2",
                                                 "TimeSince_Dx.3")
  mm <- mixed_model(fixed = PSAValue ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 + DxAge,
                    data = data, 
                    random = ~ TimeSince_Dx.1 + TimeSince_Dx.2 + 
                      TimeSince_Dx.3 | CISNET_ID,
                    family = students.t(df = df), n_phis = 1,
                    initial_values = list("betas" = gaussian()))
  save(mm, file=mixed_path)
}