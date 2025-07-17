
g_boot <- function(data, length, montecarlo, seed = 0) {
  
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