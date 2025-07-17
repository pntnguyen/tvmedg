
gform_single <- function(subdat, length, fitR, ay, am) {
  
  rFunc <- function(mod, ndat) {
    pred_prob <- predict(mod, newdata = ndat, type = "response")
    return(rbinom(1, size = 1, prob = pred_prob))
  }
  
  id     <- subdat$idsim[1]
  id_ori <- subdat$id[1]
  
  
  mm <- seq_len(length)
  
  M1p  <- numeric(length)
  Yp   <- numeric(length)
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
      dT1p   <- transform(dM1p, A = ay, Al1 = ay, Al2 = ay, M1 = M1p[l])
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

gform_single <- cmpfun(gform_single)

