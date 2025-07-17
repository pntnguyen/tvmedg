data = dat
montecarlo = 10000
g_boot <- function(data, length, montecarlo= 100, seed = 123) {
  
  df <- as.data.table(df_prep(data))
  
  set.seed(seed)
  
  # Resampling based on id and store in `boot` dataset
  clusters <- unique(df$id)
  samples  <- sample(clusters, length(clusters), replace = TRUE)
  bb       <- table(samples)
  
  #— bootstrap
  if (seed == 0) {
    # no bootstrap
    boot <- copy(df)
  } else {
    maxbb   <- max(bb)
    out_list <- vector("list", maxbb)
    
    for (zzz in seq_len(maxbb)) {
      # IDs drawn at least zzz times
      ids_zzz <- names(bb)[bb >= zzz]
      cc      <- df[id %in% ids_zzz]
      cc[, bid := paste0(id, zzz)]
      out_list[[zzz]] <- cc
    }
    
    # one single bind of all “layers”
    boot <- rbindlist(out_list, use.names = TRUE)
  }
  
  
  # Check boot data
  cat("Resampled data identical to original for seed =", seed, "?:", identical(boot, df),'\n')
  
  
  # # Centering and scaling the time-scale variable
  # boot$jj <- scale(boot$j)
  # mean_j <- attributes(boot$jj)$`scaled:center`
  # sd_j <- attributes(boot$jj)$`scaled:scale`
  # boot$jj <- as.numeric(boot$jj)
  # 
  #----- fit parametric models for
  #--- Mediator models
  
  fitM1 <- glm(M1 ~ A + Al1 + Al2 + M1l1 + M1l2 + T1l1 + T1l2 + T2l1 + T2l2 + T3l1 + T3l2 + 
                 v1 + v2 + v3 + v4 + bs(j, df = 3), 
               family = binomial, data = boot)
  
  # PseudoR2(fitM1)
  #--- Covariate models
  
  fitT1 <- glm(T1 ~ A + Al1 + Al2 + M1 + M1l1 + M1l2 + T1l1 + T1l2 + T2l1 + T2l2 + T3l1 + T3l2 + 
                 v1 + v2 + v3 + v4 + bs(j, df = 3), 
               family = binomial, data = boot)
  
  # PseudoR2(fitT1)
  
  fitT2 <- lm(T2 ~ A + Al1 + Al2 + M1 + M1l1 + M1l2 + T1 + T1l1 + T1l2 + T2l1 + T2l2 + T3l1 + T3l2 + 
                v1 + v2 + v3 + v4 + bs(j, df = 3), data = boot)
  
  fitT3 <- lm(T3 ~ A + Al1 + Al2 + M1 + M1l1 + M1l2 + T1 + T1l1 + T1l2 + T2 + T2l1 + T2l2 + T3l1 + T3l2 + 
                v1 + v2 + v3 + v4 + bs(j, df = 3), data = boot)
  
  # Outcome model: E(Y|a, m, l, v)
  
  fitY <- glm(Y ~ A + Al1 + Al2 + M1 + M1l1 + M1l2 + T1 + T1l1 + T1l2 + T2 + T2l1 + T2l2 + 
                T3 + T3l1 + T3l2 +  v1 + v2 + v3 + v4 + bs(j, df = 3), 
              family = binomial, data = boot)
  
  fitR <- list(fitM1, fitY, fitT1, fitT2, fitT3)
  
  
  # Avoid growing object
  
  df0 <- boot[j == 1]
  df0[, idn := .I]
  samples <- sample(df0$idn, size = montecarlo, replace = TRUE)
  bb <- table(samples)
  
  MC_list <- lapply(as.integer(names(bb)), function(idn_val) {
    reps <- bb[as.character(idn_val)]
    dt  <- df0[idn == idn_val]
    dt_rep <- dt[rep(1, reps), ]
    dt_rep[, rep := seq_len(reps)]
    dt_rep
  })
  
  MC <- rbindlist(MC_list, idcol = "idsim")
  MC[, idsim := seq_len(.N)]
  
  
  
  # split once (Key improvement)
  MC_list2 <- split(MC, by = "idsim", keep.by = TRUE)
  
  # Monte Carlo in parallel
  result <- foreach(
    subdat   = MC_list2,
    .combine = rbind,
    .packages = c("splines", "data.table"),
    .export   = c("gform_single")
  ) %dopar% {
    rbind(
      gform_single(subdat, length, fitR, ay = 1, am = 1),
      gform_single(subdat, length, fitR, ay = 1, am = 0),
      gform_single(subdat, length, fitR, ay = 0, am = 0)
    )
  }
  
  
}


ay = 1
am = 1
length = 12


gform_single <- function(subdat = MC_list2, length, fitR, ay, am) {
  
  rFunc <- function(mod, ndat) {
    pred_prob <- predict(mod, newdata = ndat, type = "response")
    return(rbinom(1, size = 1, prob = pred_prob))
  }
  
  id     <- subdat$idsim[1]
  id_ori <- subdat$id[1]
  
  
  mm <- seq_len(length)
  
  Yp   <- numeric(length)
  
  M1p  <- numeric(length)
  
  T1p  <- numeric(length)
  T2p  <- numeric(length)
  T3p  <- numeric(length)
  
  T1mp <- numeric(length)
  T2mp <- numeric(length)
  T3mp <- numeric(length)
  
  # Baseline covariates
  Vp      <- subdat[1, c("v1","v2","v3","v4")]
  
  M1p[1] <- 0
  Yp[1]  <- 0
  
  # Time-varying covariates contribute to outcome model
  T1p[1] <- subdat$T1
  T2p[1] <- subdat$T2
  T3p[1] <- subdat$T3
  
  # Time-varying covariates contribute to mediator model
  T1mp[1] <- subdat$T1
  T2mp[1] <- subdat$T2
  T3mp[1] <- subdat$T3
  
  

  for (l in 2:length) {
    
    if (Yp[l-1] == 1) {
      # event occurred at time (l-1) → stop here
      actual_length <- l - 1
      break
    }
    
    # compute 1- and 2-step lags
    i1  <- max(1, l-1)
    i2  <- max(1, l-2)
    
    T1l1  <- T1p[i1]
    T2l1  <- T2p[i1]
    T3l1  <- T3p[i1]
    
    T1l2  <- T1p[i2] 
    T2l2  <- T2p[i2]
    T3l2  <- T3p[i2]
    
    T1ml1 <- T1mp[i1]
    T2ml1 <- T2mp[i1]
    T3ml1 <- T3mp[i1]
    
    T1ml2 <- T1mp[i2]
    T2ml2 <- T2mp[i2]
    T3ml2 <- T3mp[i2]
    
    M1l1  <- M1p[i1] 
    M1l2  <- M1p[i2]
    
    # Predict mediator
    dM1p <- data.frame(Vp, A = am, Al1 = am, Al2 = am, M1l1, M1l2, T1l1 = T1ml1, T1l2 = T1ml2, 
                       T2l1 = T2ml1, T2l2 = T2ml2, T3l1 = T3ml1, T3l2 = T3ml2,
                       j = l)
    
    M1p[l] <- if (M1p[l-1] == 0) rFunc(fitR[[1]], dM1p) else 1
    
    
    # (2) covariates for mediator model
    dT1mp <- transform(dM1p, M1 = M1p[l])
    T1mp[l] <- rFunc(fitR[[3]], dT1mp)
    
    dT2mp <- transform(dT1mp, T1 = T1mp[l])
    T2mp[l] <- predict(fitR[[4]], newdata = dT2mp)
    
    dT3mp <- transform(dT2mp, T2 = T2mp[l])
    T3mp[l] <- predict(fitR[[5]], newdata = dT3mp)
    
    if (ay != am) {
      dT1p   <- transform(dM1p, A = ay, Al1 = ay, Al2 = ay, M1 = M1p[l],
                          T1l1  = T1p[i1],T2l1  = T2p[i1],T3l1  = T3p[i1],
                          T1l2  = T1p[i2],T2l2  = T2p[i2],T3l2  = T3p[i2]) 
        
      T1p[l] <- rFunc(fitR[[3]], dT1p)
      
      dT2p   <- transform(dT1p, T1 = T1p[l])
      T2p[l] <- predict(fitR[[4]], newdata = dT2p)
      
      dT3p   <- transform(dT2p, T2 = T2p[l])
      T3p[l] <- predict(fitR[[5]], newdata = dT3p)
      
    } else {
      T1p[l] <- T1mp[l]
      T2p[l] <- T2mp[l]
      T3p[l] <- T3mp[l]
    }
    
    # Y
    dYp <- data.frame(Vp, A = ay, Al1 = ay, Al2 = ay, M1 = M1p[l], M1l1, M1l2, 
                      T1 = T1p[l], T1l1, T1l2, 
                      T2 = T2p[l], T2l1, T2l2,
                      T3 = T3p[l], T3l1, T3l2, 
                      j = l)
    Yp[l] <- rFunc(fitR[[2]], dYp)
  }
  
  if (!exists("actual_length")) {
    actual_length <- length
  }
  
  
  mm   <- mm[1:actual_length]
  M1p  <- M1p[1:actual_length]
  Yp   <- Yp[1:actual_length]
  T1p  <- T1p[1:actual_length]
  T2p  <- T2p[1:actual_length]
  T3p  <- T3p[1:actual_length]
  T1mp <- T1mp[1:actual_length]
  T2mp <- T2mp[1:actual_length]
  T3mp <- T3mp[1:actual_length]
  
  gdat <- data.frame(id, id_ori, mm, Ay = ay, Am = am, M1p, Yp,
                     T1mp, T1p, T2mp, T2p, T3mp, T3p, Vp)
  gdat$lastid <- as.numeric(!duplicated(gdat$id, fromLast = T))
  
  return(gdat)
}
