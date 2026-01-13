mcicjm.fit <- function(df = 3, data, biopsy,
                       summary = T, n.adapt = 3000,
                       n.burnin = 3000, n.iter = 10000, 
                       prior = "unif",
                       sens_a = 0.6, sens_b = 0.9,
                       qp = 7, kntbh_spec = "equal", 
                       df_pspline = 5,
                       reduced_mcmc = FALSE, itersaved_perchain = 100,
                       GLMM_fit = NULL,
                       seed) {
  
  # Check parameters --------------------------------------------------------
  if (!(qp %in% c(7, 15))) {
    stop("The number of quadrature points can either be 7 or 15!")
  }
  
  if (!(kntbh_spec %in% c("equal", "quantile"))) {
    stop("The knots in the baseline hazard function can either be 'equal' or 'quantile'!")
  }
  
  if (!(prior %in% c("unif", "beta"))) {
    stop("The sensitivity prior can either be 'unif' or 'beta'!")
  }
  
  # Load packages -----------------------------------------------------------
  source("R script/Functions/function.R")
  
  # Start with the mixed model ----------------------------------------------
  if (is.null(GLMM_fit)) {
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
    betaL_inits <- fixef(mm)
    b_inits <- ranef(mm)
    inv_D_inits <- solve(mm$D)
    tau_inits <-  1/exp(mm$phis)^2
  } else {
    betaL_inits <- GLMM_fit$betaL
    b_inits <- GLMM_fit$b
    inv_D_inits <- GLMM_fit$inv_D
    tau_inits <-  GLMM_fit$tau
  }
  
  # Data preparation --------------------------------------------------------
  
  ## quadrature points
  f.org <- JMbayes:::gaussKronrod(qp)$sk
  w <- JMbayes:::gaussKronrod(qp)$wk
  w <- w[order(f.org)]
  f.org <- f.org[order(f.org)]
  
  ## sample size and times
  data.id <- data[!duplicated(data$CISNET_ID),]
  N_pt <- nrow(data.id)
  N_bp <- sapply(1:N_pt, function(i) {
    sum(biopsy$CISNET_ID == i)
  })
  lgt <- data$TimeSince_Dx
  evt <- matrix(NA, nrow(data.id), max(table(biopsy$CISNET_ID)))
  for (i in 1:N_pt) {
    if (N_bp[i] != 0) {
      evt[i,1:N_bp[i]] <- biopsy$TimeSince_Dx[biopsy$CISNET_ID == i]
    } else {
      evt[i,1] <- 0
    }
    
  }
  
  ### Longitudinal
  age <- data$DxAge
  num <- as.vector(c(1, 1+ cumsum(tapply(data$CISNET_ID, data$CISNET_ID, length))))
  Y <- data$PSAValue
  XL <- array(1, dim = c(nrow(data), 5))
  knot.longi <- c(0,
                  as.vector(quantile(lgt, probs = seq(0,1,length.out=4)))[2:3],
                  max(data.id$time.cmp2))
  
  XL[,2:4] <- ns(lgt, knots = knot.longi[2:3], B = knot.longi[c(1,4)])
  XL[,5] <- data$DxAge
  ZL <- XL[,1:4]
  
  ### Time-to-event
  delta <- cbind(as.integer(data.id$status.cmp == 1),
                 as.integer(data.id$status.cmp == 2),
                 as.integer(data.id$status.cmp == 0))
  X <- data.id$density
  age.id <- data.id$DxAge
  knot.surv <- get_knots(df_pspline - 4, c(data.id$time.cmp2, data.id$time.cmp1), 
                         f.org, bh = kntbh_spec)
  
  ## Treatment part for all three scenarios
  evt.trt <- data.id$time.cmp2
  evt.trt.qt <- sapply(evt.trt, function(x) {Qt(x, f.org)})
  T.trt.s <- array(NA, c(N_pt, qp, df_pspline)) # BH in the survival
  for (i in 1:N_pt) {
    T.trt.s[i,,] <- bsp_dm(evt.trt.qt[,i], knot.surv)
  }
  
  T.trt.h <- bsp_dm(evt.trt, knot.surv) # BH in the hazard function: 833*df_pspline
  
  XL.trt.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
  # N_patient, quadrature, covartiates, slope
  for (i in 1:N_pt) {
    XL.trt.s[i,,2:4,1] <- ns_dm(evt.trt.qt[,i], knot.longi)
    XL.trt.s[i,,2:4,2] <- ns_dm(evt.trt.qt[,i] - 1, knot.longi)
  }
  XL.trt.s[,,5,] <- age.id
  ZL.trt.s <- XL.trt.s[,,1:4,]  # ZL in the survival
  
  XL.trt.h <- array(1, c(N_pt, 5, 2)) # XL in the hazard
  # N_patient, covartiates, slope
  for (i in 1:N_pt) {
    XL.trt.h[i,2:4,1] <- ns_dm(evt.trt[i], knot.longi)
    XL.trt.h[i,2:4,2] <- ns_dm(evt.trt[i] - 1, knot.longi)
  }
  XL.trt.h[,5,] <- age.id
  ZL.trt.h <- XL.trt.h[,1:4,]
  
  ## Progression for the exact survival
  evt.prg <- data.id$time.cmp1
  evt.prg.qt <- sapply(evt.prg, function(x) {Qt(x, f.org)})
  T.prg.s <- array(NA, c(N_pt, qp, df_pspline)) # BH in the survival
  for (i in 1:N_pt) {
    T.prg.s[i,,] <- bsp_dm(evt.prg.qt[,i], knot.surv)
  }
  
  XL.prg.s <- array(1, c(N_pt, qp, 5, 2)) # XL in the survival
  # N_patient, quadrature, covariates, slope
  for (i in 1:N_pt) {
    XL.prg.s[i,,2:4,1] <- ns_dm(evt.prg.qt[,i], knot.longi)
    XL.prg.s[i,,2:4,2] <- ns_dm(evt.prg.qt[,i] - 1, knot.longi)
  }
  XL.prg.s[,,5,] <- age.id
  ZL.prg.s <- XL.prg.s[,,1:4,]  # ZL in the survival
  
  ## Progression for the interval censoring part
  
  # t.diff
  evt.prginc.diff <- evt - cbind(0,evt)[,1:ncol(evt)]
  
  # hazard part
  evt.prginc.qt <- array(NA, c(N_pt, max(N_bp), qp)) # quadrature timepoints in the outer integral, bs_{l-1} to bs_l
  T.prginc.h <- array(0, c(N_pt, max(N_bp), qp, df_pspline))
  XL.prginc.h <- array(0, c(N_pt, max(N_bp), qp, 5, 2))
  N_bp <- sapply(N_bp, function(x) {ifelse(x ==0, 
                                           x + 1,
                                           x)}) # take care of patients who do not do any biopsy
  for (i in 1:N_pt) {
    for (j in 1:N_bp[i]) {
      evt.prginc.qt[i,j,] <- Qt(rt = evt[i,j], lt = ifelse(j==1, 0, evt[i,j-1]), qkq=f.org)
      T.prginc.h[i,j,,] <- bsp_dm(evt.prginc.qt[i,j,], knot.surv)
      XL.prginc.h[i,j,,2:4,1] <- ns_dm(evt.prginc.qt[i,j,], knot.longi)
      XL.prginc.h[i,j,,2:4,2] <- ns_dm(evt.prginc.qt[i,j,] - 1, knot.longi)
      XL.prginc.h[i,j,,1,] <- 1
      XL.prginc.h[i,j,,5,] <- age.id[i]
    }
  }
  ZL.prginc.h <- XL.prginc.h[,,,1:4,]
  
  # survival part
  evt.prginc.qt.qt <- array(NA, c(N_pt, max(N_bp), qp, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
  T.prginc.s <- array(0, c(N_pt, max(N_bp), qp, df_pspline, qp)) # first qp is the qds from the outer integral; second is from [0,qds]
  XL.prginc.s <- array(0, c(N_pt, max(N_bp), qp, 5, qp, 2)) # first qp is the qds from the outer integral; second is from [0,qds]
  for (i in 1:N_pt) {
    for (j in 1:N_bp[i]) {
      for (m in 1:qp) {
        evt.prginc.qt.qt[i,j,m,] <- Qt(rt = evt.prginc.qt[i,j,m], qkq=f.org)
        T.prginc.s[i,j,m,,] <- t(bsp_dm(evt.prginc.qt.qt[i,j,m,], knot.surv))
        XL.prginc.s[i,j,m,2:4,,1] <- t(ns_dm(evt.prginc.qt.qt[i,j,m,], knot.longi))
        XL.prginc.s[i,j,m,2:4,,2] <- t(ns_dm(evt.prginc.qt.qt[i,j,m,] - 1, knot.longi))
        XL.prginc.s[i,j,,1,,] <- 1
        XL.prginc.s[i,j,,5,,] <- age.id[i]
      }
    }
  }
  ZL.prginc.s <- XL.prginc.s[,,,1:4,,]
  
  
  # Model specification -----------------------------------------------------
  
  ### Model script
  if (prior == "unif") {
    prior_script <- "sens ~ dunif(prior_a, prior_b)\n"
  } else if (prior == "beta") {
    prior_script <- "sens ~ dbeta(prior_a, prior_b)\n"
  }
  Model_1 <- "model {
    #### longitudinal part
    for (i in 1:N_pt) {
      for (j in num[i]:(num[i+1]-1)) {
        Y[j] ~ dt(mu[j], tau, 3)
        mu[j] <- inprod(XL[j,], betaL[]) + inprod(ZL[j,], b[i,])
      }
    
      #### Time-to-event part
      
      ### Treatment part
      ## Survival part
      log.s.trt[i] <- - evt.trt[i]/2 * 
                          (t(
                            exp(1) ^ (T.trt.s[i,1:qp,] %*% gambh[,2] + X[i] * gamma[2] +
                                  alpha[1,2] * (
                                        XL.trt.s[i,1:qp,,1] %*% betaL[] + 
                                        ZL.trt.s[i,1:qp,,1] %*% b[i,]) + 
                                  alpha[2,2] * (
                                        XL.trt.s[i,1:qp,,1] %*% betaL[] + 
                                        ZL.trt.s[i,1:qp,,1] %*% b[i,] - 
                                        XL.trt.s[i,1:qp,,2] %*% betaL[] - 
                                        ZL.trt.s[i,1:qp,,2] %*% b[i,]))
                          ) %*% w )
      
      ## Hazard part
      log.h.trt[i] <- inprod(T.trt.h[i, ], gambh[,2]) + X[i] * gamma[2] + 
                    alpha[1,2] * (
                      inprod(XL.trt.h[i,,1], betaL[]) + inprod(ZL.trt.h[i,,1], b[i,])) +
                    alpha[2,2] * (
                      inprod(XL.trt.h[i,,1], betaL[]) + inprod(ZL.trt.h[i,,1], b[i,]) - 
                      inprod(XL.trt.h[i,,2], betaL[]) - inprod(ZL.trt.h[i,,2], b[i,]))
                      
      ### Exact progression part
      ## Survival part
      s.prg[i] <- exp(- evt.prg[i]/2 * 
                          (t(
                            exp(1) ^ (T.prg.s[i,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +
                                  alpha[1,1] * (
                                        XL.prg.s[i,1:qp,,1] %*% betaL[] + 
                                        ZL.prg.s[i,1:qp,,1] %*% b[i,]) + 
                                  alpha[2,1] * (
                                        XL.prg.s[i,1:qp,,1] %*% betaL[] + 
                                        ZL.prg.s[i,1:qp,,1] %*% b[i,] - 
                                        XL.prg.s[i,1:qp,,2] %*% betaL[] - 
                                        ZL.prg.s[i,1:qp,,2] %*% b[i,]))
                          ) %*% w))
      
      ### Interval-censored progression part
      ## Survival part
      for (j in 1:N_bp[i]) {
      
        ### Biopsy sensitivity
        prob.det[i,j]  <- ifelse(delta[i,1] == 1,
          sens * (1-sens)^(N_bp[i] - j),
          (1-sens)^(N_bp[i] - j + 1))
      
        for (m in 1:qp) {
          s.prginc.eint.gk[i,j,m] <- exp(1) ^ 
                                      (- evt.prginc.qt[i,j,m]/2 * 
                                          inprod(w, 
                                            exp(1) ^ (t(T.prginc.s[i,j,m,,1:qp]) %*% gambh[,1] + X[i] * gamma[1] + 
                                                  alpha[1,1] * (t(XL.prginc.s[i,j,m,,1:qp,1]) %*% betaL[] + 
                                                t(ZL.prginc.s[i,j,m,,1:qp,1]) %*% b[i,]) + 
                                                  alpha[2,1] * (
                                                t(XL.prginc.s[i,j,m,,1:qp,1]) %*% betaL[] + 
                                                t(ZL.prginc.s[i,j,m,,1:qp,1]) %*% b[i,] - 
                                                t(XL.prginc.s[i,j,m,,1:qp,2]) %*% betaL[] - 
                                                t(ZL.prginc.s[i,j,m,,1:qp,2]) %*% b[i,]))
                                          ))
        }
                                       
        prginc.eint[i,j] <- evt.prginc.diff[i,j]/2 * 
                              inprod(w, 
                                     s.prginc.eint.gk[i,j,] * 
                                        exp(1) ^ (T.prginc.h[i,j,1:qp,] %*% gambh[,1] + X[i] * gamma[1] +  
                                     alpha[1,1] * (
                                       XL.prginc.h[i,j,1:qp,,1] %*% betaL[] +
                                       ZL.prginc.h[i,j,1:qp,,1] %*% b[i,]) + 
                                     alpha[2,1] * (
                                       XL.prginc.h[i,j,1:qp,,1] %*% betaL[] +
                                       ZL.prginc.h[i,j,1:qp,,1] %*% b[i,] - 
                                       XL.prginc.h[i,j,1:qp,,2] %*% betaL[] -
                                       ZL.prginc.h[i,j,1:qp,,2] %*% b[i,]))
                                     ) * prob.det[i,j]  
      }
      prginc[i] <- ifelse(sum(prginc.eint[i,1:N_bp[i]]) < 1e-16,
                          1e-16,
                          sum(prginc.eint[i,1:N_bp[i]]))
      
      #### LIKELIHOOD
      loglike[i] <- delta[i,1] * (log(prginc[i]) + log.s.trt[i]) + 
                    delta[i,2] * (log.h.trt[i] + log.s.trt[i] + log(s.prg[i] + prginc[i])) +
                    delta[i,3] * (log.s.trt[i] + log(s.prg[i] + prginc[i]))
      phi[i] <- 5000 - loglike[i]
      zeros[i] ~ dpois(phi[i])
      
      b[i,1:Nb] ~ dmnorm(mub[], inv_D[,])
    }
    
    #### PRIORS
    tau ~ dgamma(0.01,0.01)
    inv_D[1:Nb,1:Nb] ~ dwish(4 * diagtb[,], Nb+1)
    for (l in 1:Nb) {
      for (k in 1:Nb) {diagtb[l,k] <- ifelse(l == k, tb[l], 0)}
      tb[l] ~ dgamma(0.5, 0.01)
    }
  
    D[1:Nb,1:Nb] <- inverse(inv_D[,]) \n"
  
  Model_2 <- 
    "for (l in 1:Nbeta) {betaL[l] ~ dnorm(0, 0.01)}
    for (k in 1:Nrisk) {
      gambh[1:Ngambh,k] ~ dmnorm(mu.gambh[], tau.smooth[k] * tau.gambh[,])
      tau.smooth[k] ~ dgamma(5, 0.5)
      alpha[1,k] ~ dnorm(0, 0.01)
      alpha[2,k] ~ dnorm(0, 0.01)
      gamma[k] ~ dnorm(0, 0.01)
    }
  }"
  
  Model <- paste0(Model_1, prior_script, Model_2)
  
  # Data input --------------------------------------------------------------
  DD <- diag(df_pspline)
  
  data <- list(N_pt = N_pt, N_bp = N_bp, num = num, qp = qp,
               Y = Y, XL = XL, ZL = ZL,
               T.trt.s = T.trt.s, X = X, XL.trt.s = XL.trt.s, ZL.trt.s = ZL.trt.s,
               evt.trt = evt.trt, w = w, 
               T.trt.h = T.trt.h, XL.trt.h = XL.trt.h, ZL.trt.h = ZL.trt.h,
               T.prg.s = T.prg.s, XL.prg.s = XL.prg.s, ZL.prg.s = ZL.prg.s, 
               evt.prg= evt.prg,
               T.prginc.s = T.prginc.s, XL.prginc.s = XL.prginc.s, ZL.prginc.s = ZL.prginc.s,
               evt.prginc.qt = evt.prginc.qt,
               T.prginc.h = T.prginc.h, XL.prginc.h = XL.prginc.h, ZL.prginc.h = ZL.prginc.h,
               evt.prginc.diff = evt.prginc.diff, #prob.det = prob.det,
               prior_a = sens_a, prior_b = sens_b,
               delta = delta, zeros = rep(0, N_pt),
               mub = rep(0,dim(ZL)[2]), Nb = dim(ZL)[2], Nbeta = dim(XL)[2], 
               Nrisk = ncol(delta) - 1, Ngambh = df_pspline,
               mu.gambh = rep(0,df_pspline), tau.gambh = crossprod(diff(DD, diff = 2)) + 1e-6* DD)
  #b = ranef(mm), betaL = fixef(mm))
  
  
  # Inits -------------------------------------------------------------------
  set.seed(seed)
  rng.name <- sample(c("base::Wichmann-Hill",
                       "base::Marsaglia-Multicarry",
                       "base::Super-Duper",
                       "base::Mersenne-Twister"), 3, replace = T)
  rng.num <- sample(1:100000, 3, replace = F)
  
  initials <- list(
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= betaL_inits,
         gambh = matrix(rnorm(df_pspline*2),df_pspline,2),
         gamma = rnorm(2),
         b = b_inits,
         inv_D = inv_D_inits,
         tau = tau_inits,
         .RNG.name = rng.name[1],
         .RNG.seed = rng.num[1]),
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= betaL_inits,
         gambh = matrix(rnorm(df_pspline*2),df_pspline,2),
         gamma = rnorm(2),
         b = b_inits,
         inv_D = inv_D_inits,
         tau = tau_inits,
         .RNG.name = rng.name[2],
         .RNG.seed = rng.num[2]),
    list(alpha = matrix(rnorm(4), 2, 2),
         betaL= betaL_inits,
         gambh = matrix(rnorm(df_pspline*2),df_pspline,2),
         gamma = rnorm(2),
         b = b_inits,
         inv_D = inv_D_inits,
         tau = tau_inits,
         .RNG.name = rng.name[3],
         .RNG.seed = rng.num[3])
  )
  
  # run the model -----------------------------------------------------------
  out <- lapply(1:3, function(i) {
    future::future({
      mcicjm <- jags.model(file = textConnection(Model), data = data,
                           inits = initials[[i]], n.chains = 1, n.adapt = n.adapt,
                           quiet = T)  #
      update(mcicjm, n.iter = n.burnin)
      mcmc <- coda.samples(mcicjm, variable.names = c("alpha", "gambh", "betaL",
                                                      "D","gamma", "tau", "sens"), # ,  
                           n.iter = n.iter, thin = 10)
      return(list(model = mcicjm,
                  mcmc = mcmc))
    })
  })
  
  samples.mcicjm <- lapply(out, future::value)
  
  #print("JAGS model done!")
  iter.perchain <- nrow(as.matrix(samples.mcicjm[[1]]$mcmc))
  if (reduced_mcmc) {
    mcmcoutput <- lapply(samples.mcicjm, function(x) {
      as.mcmc(x$mcmc[(iter.perchain-itersaved_perchain+1):iter.perchain,])
    })
  } else {
    mcmcoutput <- lapply(samples.mcicjm, function(x) {x$mcmc})
  }
  
  if (summary == F) {
    return(list(mcmc = mcmcoutput,
                model_info = list(var_names = list(id = "CISNET_ID",
                                                   longitime = "TimeSince_Dx",
                                                   survtime = c("time.cmp1", "time.cmp2"),
                                                   longi.bsVar = "DxAge",
                                                   surv.bsVar = "density"),
                                  knots = list(knot.longi = knot.longi,
                                               knot.surv = knot.surv)),
                inits = initials))
  } else if (summary == T) {
    mcmc.mat <- do.call(rbind, lapply(samples.mcicjm, 
                                      function(x) {as.matrix(x$mcmc)}))
    mcmc_all <- as.mcmc.list(lapply(1:3, function(x) 
      as.mcmc(samples.mcicjm[[x]]$mcmc)))
    return(list(mcmc = mcmcoutput,
                model_info = list(var_names = list(id = "CISNET_ID",
                                                   longitime = "TimeSince_Dx",
                                                   survtime = c("time.cmp1", "time.cmp2"),
                                                   longi.bsVar = "DxAge",
                                                   surv.bsVar = "density"),
                                  knots = list(knot.longi = knot.longi,
                                               knot.surv = knot.surv)),
                inits = initials,
                summary = list(coef = apply(mcmc.mat, 2, mean),
                               coef_lower = apply(mcmc.mat, 2, function(x) quantile(x, 0.025)),
                               coef_upper = apply(mcmc.mat, 2, function(x) quantile(x, 0.975)),
                               coef.sd = apply(mcmc.mat, 2, sd),
                               gelman_rubin = my.gelman.diag(mcmc_all),
                               se = mcse.mat(mcmc.mat)[,2])))
  } 
  
}