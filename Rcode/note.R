fitg <- function(data, seed = 0,
                 mreg = "binomial",
                 lreg = c("binomial","gaussian","gaussian"),
                 yreg = "binomial",dof = 3){
  
  cat("g-formula for single mediator",'\n')
  
  res_df <- resamp(data = data$df,seed = seed)
  
  fitR <- list()
  
  fitR$df <- res_df

  #----- fit parametric models for
  #--- Mediator models
  
  if(length(mreg) != length(data$fm)){
    stop("the defined regression of M is not equal")
  }
  
  for (i in 1:length(data$fm)){
    fitM <- paste0(data$fm[[i]],"+","bs(jj,df=",dof,")")
    fitR$M[[i]] <- glm(fitM ,family = mreg[i], data = res_df) 
  }
  
  
  #--- Covariate models
  if(length(lreg) != length(data$fl)){
    stop("the defined regression of L is not equal")
  }
  
  for (i in 1:length(data$fl)){
    fitL <- paste0(data$fl[[i]],"+","bs(jj,df=",dof,")")
    fitR$L[[i]] <- glm(fitL ,family = lreg[i], data = res_df) 
  }
  
  #--- Outcome model: 
  fitY <- paste0(data$fy,"+","bs(jj,df=",dof,")")
  
  
  fitR$Y <-  glm(fitY ,family = yreg, data = res_df) 
  
  fitR
}

fitR2 <- process_data(
  fix = c("age","sex","ow","risk"),
  expo = c("Ap"),
  med = c("Mp"),
  tvar = c("L1","L2","L3"),
  outc = c("Yp"),
  lag = 2,
  time = c("mm"),
  data = dat0
) %>% fitg(mreg = "binomial",
           lreg = c("binomial","gaussian","gaussian"),
           yreg = "binomial",dof = 3,seed = 0)

data <- process_data(
  fix = c("age","sex","ow","risk"),
  expo = c("Ap"),
  med = c("Mp"),
  tvar = c("L1","L2","L3"),
  outc = c("Yp"),
  lag = 2,
  time = c("mm"),
  data = dat0
)
fitR2$df

## pre_g

boot_g <- function(data,n_boot = 10000,repdf = 1){
  
  boot <- data$df
  # Select baseline visit
  df0 <- boot[boot$j==1, ]
  
  df0$idn <- 1:nrow(df0)
  
  MC <- NULL
  
  if (n_boot == 0) {
    MC <- do.call(rbind, replicate(repdf, df0, simplify = FALSE))
    MC$idsim <- 1:nrow(MC)
    
  } else {
    samples <- sample(df0$idn, size = n_boot, replace = T)
    bb <- table(samples)
    
    for(zzz in 1:max(bb)) {
      cc <- df0[df0$idn %in% names(bb[bb %in% c(zzz:max(bb))]), ]
      cc$bid <- paste0(cc$idn, zzz)
      MC <- rbind(MC, cc)
    }
    
    MC$idsim <- 1:n_boot
  }
  
  data$res_df <- MC %>% as_tibble()
  data
}

fitR2 <- process_data(
  fix = c("age","sex","ow","risk"),
  expo = c("Ap"),
  med = c("Mp"),
  tvar = c("L1","L2","L3"),
  outc = c("Yp"),
  lag = 2,
  time = c("mm"),
  data = dat0
) %>% fitg(mreg = "binomial",
           lreg = c("binomial","gaussian","gaussian"),
           yreg = "binomial",dof = 3,seed = 0) %>% 
  boot_g(n_boot = 0,repdf = 1)

## gform

pFunc <- function(mod, ndat) {
  as.numeric(predict(mod, newdata = ndat, type = "response") > runif(1))
}
  
rFunc <- function(mod, ndat) {
  pred_prob <- predict(mod, newdata = ndat, type = "response")
  return(rbinom(1, size = 1, prob = pred_prob))
}
  
g_form <- function(data = fitR2 , ii = 2, length = 12, am = 1, ay = 0){
  
  dddd <- data$res_df %>% data.frame()
  
  lagg <- dddd %>% dplyr::select(contains("L1l")) %>% ncol() 
  
  d2 <- dddd[dddd$idsim==ii, ]
  
  id <- d2$idsim
  id_ori <- d2$id
  
  length <- length
  
  # Baseline covariates
  Vp <- d2 %>% select(starts_with("v"))
  
  Yp <- mm <- numeric()
  
  mm[1:lagg-1] <- j <- 1
  
  Yp[1:lagg-1] <- 0
  
  timee <- cen_and_scale(data$df$j)
  
  # mediator
  Mp <- matrix(ncol = length(data$M)) %>% data.frame()
  names(Mp) <- paste0("M",1:length(data$M))
  Mp[1:lagg-1,] <- d2 %>% select(names(Mp)) 
  
  
  # time-varying covariates (contribute to mediator models)
  Lmp <- matrix(ncol = length(data$L)) %>% data.frame()
  names(Lmp) <- paste0("L",1:length(data$L))
  Lmp[1:lagg-1,] <- d2 %>% select(names(Lmp)) 
  # time-varying covariates (contribute to outcome models)
  Lp <- Lmp
  
  
  for (l in lagg:length) {
    l = 3
    if (Yp[l-1]==1) {
      break
    } else{
      
      # Predict mediator
      var_fm <- attr(data$M[[1]]$terms, "term.labels")
      var_fm <- var_fm[-length(var_fm)] 
      
      dfMp <- d2 %>% select(matches(var_fm)) %>% 
        mutate(jj = (l-timee$mean_j)/timee$sd_j)
      
      dfMp[startsWith(colnames(dfMp), "A")] <- am
      
      
      if (l > lagg){
        
        for (zz in 1:(lagg)){
          term <- paste0("l",zz)
          # L lag
          dfMp[startsWith(colnames(dfMp), "L") & endsWith(colnames(dfMp), term)] <- Lmp[l-zz,]
          # M lag
          dfMp[startsWith(colnames(dfMp), "M") & endsWith(colnames(dfMp), term)] <- Mp[l-zz,]
        }
        
      }
      
      
      for (x in 1:length(data$M)){
        M_reg <- data$M[[x]]$family$family
        
        if (Mp[l-1,x] == 0) {
          Mp[l,x] <- case_when(
            M_reg  == "binomial" ~ pFunc(data$M[[x]], dfMp),
            M_reg  == "gaussian" ~ predict(data$M[[x]], dfMp)
          )
        } else {
          Mp[l,x] <- 1
        }
        
      }
      
      # Predict time-varying covariates (contribute to mediator models)
      # L
      var_fl <- attr(data$L[[1]]$terms, "term.labels")
      var_fl <- var_fl[-length(var_fl)] 
      dfLmp <- d2 %>% select(matches(var_fl)) %>% 
        mutate(jj = (l-timee$mean_j)/timee$sd_j)
      dfLmp[startsWith(colnames(dfLmp), "A")] <- am
      dfLmp[colnames(dfLmp) == colnames(Mp)] <- Mp[l,]
      
      
      if (l > lagg){
        
        for (zz in 1:lagg){
          term <- paste0("l",zz)
          # L lag
          dfLmp[startsWith(colnames(dfLmp), "L") & endsWith(colnames(dfLmp), term)] <- Lmp[l-zz,]
          # M lag
          dfLmp[startsWith(colnames(dfLmp), "M") & endsWith(colnames(dfLmp), term)] <- Mp[l-zz,]
        }
        
      }
      
      for (x in 1:length(data$L)){
        
        L_reg <- data$L[[x]]$family$family
        
        Lmp[l,x] <- case_when(
          L_reg  == "binomial" ~ pFunc(data$L[[x]], dfLmp),
          L_reg  == "gaussian" ~ predict(data$L[[x]], dfLmp)
        )
      }
      
      
      # Predict time-varying covariates (contribute to outcome models, if ay != am)
      if (ay != am) {
        dfLp <- d2 %>% select(matches(var_fl)) %>% 
          mutate(jj = (l-timee$mean_j)/timee$sd_j)
        dfLp[startsWith(colnames(dfLp), "A")] <- ay
        dfLp[colnames(dfLp) == colnames(Mp)] <- Mp[l,]
        
        if (l > lagg){
          
          for (zz in 1:lagg){
            term <- paste0("l",zz)
            # L lag
            dfLp[startsWith(colnames(dfLp), "L") & endsWith(colnames(dfLp), term)] <- Lp[l-zz,]
            # M lag
            dfLp[startsWith(colnames(dfLp), "M") & endsWith(colnames(dfLp), term)] <- Mp[l-zz,]
          }
          
        }
        
        
        for (x in 1:length(data$L)){
          
          L_reg <- data$L[[x]]$family$family
          
          Lp[l,x] <- case_when(
            L_reg  == "binomial" ~ pFunc(data$L[[x]], dfLp),
            L_reg  == "gaussian" ~ predict(data$L[[x]], dfLp)
          )
        }
        
      } else{
        Lp <- Lmp
      }
      
      # Y
      var_y <- attr(data$Y$terms, "term.labels")
      var_y <- var_y[-length(var_y)] 
      dfYp <- d2 %>% select(matches(var_y)) %>% 
        mutate(jj = (l-timee$mean_j)/timee$sd_j)
      dfYp[startsWith(colnames(dfYp), "A")] <- ay
      dfYp[colnames(dfYp) %in% colnames(Mp)] <- Mp[l,]
      dfYp[colnames(dfYp) %in% colnames(Lp)] <- Lp[l,] 
      
      if (l > lagg){
        
        for (zz in 1:lagg){
          term <- paste0("l",zz)
          # L lag
          dfYp[startsWith(colnames(dfYp), "L") & endsWith(colnames(dfYp), term)] <- Lp[l-zz,]
          # M lag
          dfYp[startsWith(colnames(dfYp), "M") & endsWith(colnames(dfYp), term)] <- Mp[l-zz,]
        }
        
      }
      
      Yp[l] <- pFunc(data$Y, dfYp)
      
    }
    
    colnames(Lmp) <- paste0("Lmp",1:length(data$L))
    colnames(Lp) <- paste0("Lp",1:length(data$L))
    
    mm[l] <- l
  }
  boot_num <- seed
  gdat <- data.frame(boot_num,id, id_ori, mm, Ay = ay, Am = am, Mp, Yp,
                     Lmp, Lp, Vp)
  gdat$lastid <- as.numeric(!duplicated(gdat$id, fromLast = T))
  return(gdat)
}

process_data(
  fix = c("age","sex","ow","risk"),
  expo = c("Ap"),
  med = c("Mp"),
  tvar = c("L1","L2","L3"),
  outc = c("Yp"),
  lag = 2,
  time = c("mm"),
  data = dat0
) %>% fitg(mreg = "binomial",
           lreg = c("binomial","gaussian","gaussian"),
           yreg = "binomial",dof = 3,seed = 0) %>% 
  boot_g(n_boot = 0,repdf = 1)  %>% 
  g_form(ii = 1, length = 12, ay = 1, am = 1)

gform(ii = 1, pgdat = MC, length = 12, ay = 1, am = 1) %>% View()


# install.packages("furrr")
library(furrr)

plan(multisession, workers = 3)

maxiter <- ifelse(montecarlo == 0, nrow(MC), montecarlo)

cat("Number of replicates for baseline", maxiter, "\n")

resultDat1M <- foreach (iii = 1:maxiter, .combine = rbind) %dopar% {
  outdat11 <- gform(ii = iii, pgdat = MC, length = length, ay = 1, am = 1)
  outdat10 <- gform(ii = iii, pgdat = MC, length = length, ay = 1, am = 0)
  outdat00 <- gform(ii = iii, pgdat = MC, length = length, ay = 0, am = 0)
  rbind(outdat11, outdat10, outdat00)
}

gform_wrapper <- function(iii, pgdat, length) {
  outdat11 <- gform(ii = iii, pgdat = pgdat, length = length, ay = 1, am = 1)
  outdat10 <- gform(ii = iii, pgdat = pgdat, length = length, ay = 1, am = 0)
  outdat00 <- gform(ii = iii, pgdat = pgdat, length = length, ay = 0, am = 0)
  bind_rows(outdat11, outdat10, outdat00)
}

length = 100
plan(multisession, workers = 3)

resultDat1M <- future_map_dfr(
  1:nrow(MC),
  ~ gform_wrapper(.x, pgdat = MC, length = length)
)

## 

## benchmark

# install.packages("rbenchmark")
library(rbenchmark)

gform_wrapper2 <- function(iii, data, length) {
  outdat11 <- g_form(ii = iii,data=data, length = length, ay = 1, am = 1)
  outdat10 <- g_form(ii = iii,data=data, length = length, ay = 1, am = 0)
  outdat00 <- g_form(ii = iii,data=data, length = length, ay = 0, am = 0)
  bind_rows(outdat11, outdat10, outdat00)
}

resultDat2M <- future_map_dfr(
  1:100,
  ~ gform_wrapper2(ii = .x, data = fitR2, length = 12)
)
library(doParallel)
resultDat1M <- foreach (iii = 1:100, .combine = rbind) %dopar% {
  outdat11a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 1)
  outdat10a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 0)
  outdat00a <- g_form(ii = iii, data=fitR2, length = 12, ay = 0, am = 0)
  rbind(outdat11a, outdat10a, outdat00a)
}

benchmark("furrr" = {
  gform_wrapper2 <- function(iii, data, length) {
    outdat11 <- g_form(ii = iii,data=data, length = length, ay = 1, am = 1)
    outdat10 <- g_form(ii = iii,data=data, length = length, ay = 1, am = 0)
    outdat00 <- g_form(ii = iii,data=data, length = length, ay = 0, am = 0)
    bind_rows(outdat11, outdat10, outdat00)
  }
  
  resultDatM <- future_map_dfr(
    1:100,
    ~ gform_wrapper2(ii = .x, data = fitR2, length = 12)
  )
},
"dopar" = {
  resultDatM <- foreach (iii = 1:100, .combine = rbind) %dopar% {
    outdat11a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 1)
    outdat10a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 0)
    outdat00a <- g_form(ii = iii, data=fitR2, length = 12, ay = 0, am = 0)
    rbind(outdat11a, outdat10a, outdat00a)
  }
},
replications = 10,
columns = c("test", "replications", "elapsed",
            "relative", "user.self", "sys.self"))

library(bench)
mark(
  future_map_dfr(1:100,  ~ gform_wrapper2(ii = .x, data = fitR2, length = 12)),
  foreach (iii = 1:100, .combine = rbind) %dopar% {
    outdat11a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 1)
    outdat10a <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 0)
    outdat00a <- g_form(ii = iii, data=fitR2, length = 12, ay = 0, am = 0)
    rbind(outdat11a, outdat10a, outdat00a)
  })

library(foreach)
library(doParallel)
library(DescTools)
library(tidyverse)

foreach (iii = 1:10, .combine = rbind) %dopar% {
  outdat11 <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 1)
  outdat10 <- g_form(ii = iii, data=fitR2, length = 12, ay = 1, am = 0)
  outdat00 <- g_form(ii = iii, data=fitR2, length = 12, ay = 0, am = 0)
  rbind(outdat11, outdat10, outdat00)
}

g_form(ii = 1, data=fitR2, length = 12, ay = 1, am = 1)
