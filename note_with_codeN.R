set.seed(123)
fitR2 <- process_data(
  basec = c("age","sex","ow","risk"),
  expo = c("Ap"),
  med = c("Mp"),
  tvar = c("L1","L2","L3"),
  outc = c("Yp"),
  lag = 2,
  time = c("mm"),
  norev = c("Mp"),
  tvar_to_med = F,
  cont_exp = F,
  cont_exp_std = F,
  sp_list = c("mm"),
  sp_type = c("bs"),
  sp_df= c(3),
  data = dat0
) %>% fitg(boot=T,
            mreg = "binomial",
            lreg = c("binomial","gaussian","gaussian"),
            yreg = "binomial") %>% 
  baseline_mc(montecarlo = 100)


rbindlist(MC_list2)
rbindlist(fitR2$res_df)

## run parallel

g_form <- function(data = fitR2$res_df, model = fitR2, followup = 12, am = 1, ay = 0){
  
  norev_var <- model$norev_var
  
  dddd <- data %>% as.data.frame()
  id     <- dddd$idsim[1]
  id_ori <- dddd$id[1]
  
  lagg <- dddd %>% dplyr::select(contains("L1l")) %>% ncol() 
  
  # followup <- followup
  
  # Baseline covariates
  Vp <- dddd %>% select(starts_with("v"))
  
  Yp2 <- mm <- numeric()
  
  mm[1:lagg-1] <- j <- 1
  
  Yp2[1:lagg-1] <- 0
  
  
  # mediator
  Mp <- matrix(ncol = length(model$M)) %>% data.frame()
  names(Mp) <- paste0("M",1:length(model$M))
  Mp[1:lagg-1,] <- dddd %>% select(names(Mp)) 
  
  
  # time-varying covariates (contribute to mediator models)
  Lmp <- matrix(ncol = length(model$L)) %>% data.frame()
  names(Lmp) <- paste0("L",1:length(model$L))
  Lmp[1:lagg-1,] <- dddd %>% select(names(Lmp)) 
  # time-varying covariates (contribute to outcome models)
  Lp <- Lmp
  
  for (l in lagg:followup) {
    
    if (Yp2[l-1]==1) {
      break
    } else{
      
      # Predict mediator
      var_fm <- attr(model$M[[1]]$terms, "term.labels")
      
      var_fm <- gsub(".*\\(([^,]+),.*", "\\1", var_fm)
      
      var_fm <- var_fm[-length(var_fm)] 
      
      dfMp <- dddd %>% select(matches(var_fm)) %>% 
        mutate(j = l)
      
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
      
      
      for (x in 1:length(model$M)){
        
        M_reg <- model$M[[x]]$family$family
        
        if (names(Mp[x]) %in% norev_var){
          
          if (M_reg == "binomial" & Mp[l-1,x] == 1) {
            Mp[l,x] <- 1
          } else {
            Mp[l,x] <- case_when(
              M_reg  == "binomial" ~ rFunc(model$M[[x]], dfMp),
              M_reg  == "gaussian" ~ predict(model$M[[x]], dfMp)
            )
          }
          
        } else {
          
          Mp[l,x] <- case_when(
            M_reg  == "binomial" ~ rFunc(model$M[[x]], dfMp),
            M_reg  == "gaussian" ~ predict(model$M[[x]], dfMp)
          )
          
        }
        
      }
      
      # Predict time-varying covariates (contribute to mediator models)
      # L
      var_fl <- attr(model$L[[1]]$terms, "term.labels")
      var_fl <- gsub(".*\\(([^,]+),.*", "\\1", var_fl)
      var_fl <- var_fl[-length(var_fl)] 
      dfLmp <- dddd %>% select(matches(var_fl)) %>% 
        mutate(j = l)
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
      
      for (x in 1:length(model$L)){
        
        L_reg <- model$L[[x]]$family$family
        
        if (names(Lmp[x]) %in% norev_var){
          if (L_reg == "binomial" & Lmp[l-1,x] == 1){
            
            Lmp[l,x] <- 1
            
          } else {
            
            Lmp[l,x] <- case_when(
              L_reg  == "binomial" ~ rFunc(model$L[[x]], dfLmp),
              L_reg  == "gaussian" ~ predict(model$L[[x]], dfLmp)
            )
            
          }
          
        } else{
          
          Lmp[l,x] <- case_when(
            L_reg  == "binomial" ~ rFunc(model$L[[x]], dfLmp),
            L_reg  == "gaussian" ~ predict(model$L[[x]], dfLmp)
          )
          
        }
        
      }
      
      # Predict time-varying covariates (contribute to outcome models, if ay != am)
      if (ay != am) {
        dfLp <- dddd %>% select(matches(var_fl)) %>% 
          mutate(j = l)
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
        
        for (x in 1:length(model$L)){
          
          L_reg <- model$L[[x]]$family$family
          
          if (names(Lp[x]) %in% norev_var){
            if (L_reg == "binomial" & Lp[l-1,x] == 1){
              
              Lp[l,x] <- 1
              
            } else {
              
              Lp[l,x] <- case_when(
                L_reg  == "binomial" ~ rFunc(model$L[[x]], dfLp),
                L_reg  == "gaussian" ~ predict(model$L[[x]], dfLp)
              )
              
            }
            
          } else {
            
            Lp[l,x] <- case_when(
              L_reg  == "binomial" ~ rFunc(model$L[[x]], dfLp),
              L_reg  == "gaussian" ~ predict(model$L[[x]], dfLp)
            )
            
          }
          
        }
        
      } else{
        Lp <- Lmp
      }
      
      # Y
      var_y <- attr(model$Y$terms, "term.labels")
      var_y <- gsub(".*\\(([^,]+),.*", "\\1", var_y)
      var_y <- var_y[-length(var_y)] 
      dfYp <- dddd %>% select(matches(var_y)) %>% 
        mutate(j = l)
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
      
      Yp2[l] <- rFunc(model$Y, dfYp)
      
    }
    
    mm[l] <- l
  }
  
  colnames(Lmp) <- paste0("Lmp",1:length(model$L))
  colnames(Lp) <- paste0("Lp",1:length(model$L))
  
  # boot_num <- seed
  gdat2 <- data.frame(id, id_ori, mm, Ay = ay, Am = am, Mp, Yp2,
                      Lmp, Lp, Vp)
  gdat2$lastid <- as.numeric(!duplicated(gdat2$id, fromLast = T))
  
  return(gdat2)
}


result %>% filter(Ay == 1 & Am == 0) %>% View()

rbindlist(fitR2$res_df)
gdat2
# library(furrr)
future_map_dfr(fitR2$res_df,  
               ~ g_form(data = .,model = fitR2,followup = 12,ay = 1,am = 0)) 


tvmedg <- function(data,basec,expo,med,tvar,outc,time,lag = 2,
                   norev = NULL, cont_exp = NULL,cont_exp_std = F,
                   tvar_to_med = F,
                   mreg = "binomial",
                   lreg = c("binomial","gaussian","gaussian"),
                   yreg = "binomial",
                   sp_list = NULL,sp_type = NULL,sp_df= NULL,
                   followup = 12,
                   seed = 0,montecarlo = 10,boot = FALSE,nboot = 1,ci = .95,
                   parallel=TRUE){
  
  set.seed(seed)
  
  start_time <- Sys.time()
  
  qqq <- matrix(ncol = 3) %>% data.frame()
  colnames(qqq) <- c("mQ11","mQ10","mQ00")
  
  qqq_ci <- matrix(ncol = 3) %>% data.frame()
  colnames(qqq_ci) <- c("mQ11","mQ10","mQ00")
  
  
  ## point estimate
  
  fitR2 <- process_data(
    basec = basec,
    expo = expo,
    med = med,
    tvar = tvar,
    outc = outc,
    lag = lag,
    time = time,
    norev = norev,
    tvar_to_med = tvar_to_med,
    cont_exp = cont_exp,
    cont_exp_std = cont_exp_std,
    sp_list = sp_list,
    sp_type = sp_type,
    sp_df = sp_df,
    data = data
  )  %>% 
    fitg(boot=boot,
         mreg = mreg,
         lreg = lreg,
         yreg = yreg) %>% 
    baseline_mc(montecarlo = montecarlo)
  
  am <- fitR2$am
  
  
  if (parallel == TRUE){
    
    gform_wrapper2 <- function(data, model) {
      outdat11 <- g_form(data=data, model = model,followup = followup, ay = am+1, am = am+1)
      outdat10 <- g_form(data=data, model = model,followup = followup, ay = am+1, am = am)
      outdat00 <- g_form(data=data, model = model,followup = followup, ay = am, am = am)
      bind_rows(outdat11, outdat10, outdat00)
    }
    
    resultDatM <- future_map_dfr(fitR2$res_df,  
                                 ~ gform_wrapper2(data = .,model = fitR2))
  } else {
    resultDatM <- data.frame()
    for (iii in 1:montecarlo){
      outdat11 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am+1, am = am+1)
      outdat10 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am+1, am = am)
      outdat00 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am, am = am)
      resultDatM2 <- rbind(outdat11, outdat10, outdat00)
      resultDatM <- rbind(resultDatM,resultDatM2)
    }
  }
  
  qqq <- ExtResult2(resultDatM,am = am) %>% mutate(
    rIE_b = mQ11 - mQ10,
    rDE_b = mQ10 - mQ00,
    rTE_b = mQ11 - mQ00,
    rPE_b = rIE_b/ rTE_b
  )
  
  if (boot == TRUE){
    
    for (it in 1:nboot){
      
      ## boostrap
      fitR2a <- process_data(
        basec = basec,
        expo = expo,
        med = med,
        tvar = tvar,
        outc = outc,
        lag = lag,
        time = time,
        norev = norev,
        tvar_to_med = tvar_to_med,
        cont_exp = cont_exp,
        cont_exp_std = cont_exp_std,
        sp_list = sp_list,
        sp_type = sp_type,
        sp_df = sp_df,
        data = data)  %>% 
        fitg(boot=boot,
             mreg = mreg,
             lreg = lreg,
             yreg = yreg) %>% 
        baseline_mc(montecarlo = montecarlo)
      
      am_ci <- fitR2a$am
      
      ## extract mean of q11,q10,q00 of the ith iter
      if (parallel == TRUE){
        
        gform_wrapper2 <- function(data, model) {
          outdat11 <- g_form(data=data, model = model,followup = followup, ay = am+1, am = am+1)
          outdat10 <- g_form(data=data, model = model,followup = followup, ay = am+1, am = am)
          outdat00 <- g_form(data=data, model = model,followup = followup, ay = am, am = am)
          bind_rows(outdat11, outdat10, outdat00)
        }
        
        resultDatM_ci <- future_map_dfr(fitR2$res_df,  
                                        ~ gform_wrapper2(data = .,model = fitR2))
      } else {
        resultDatM <- data.frame()
        for (iii in 1:montecarlo){
          outdat11 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am+1, am = am+1)
          outdat10 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am+1, am = am)
          outdat00 <- g_form(data=fitR2$res_df, model = fitR2,followup = followup, ay = am, am = am)
          resultDatM2_ci <- rbind(outdat11, outdat10, outdat00)
          resultDatM_ci <- rbind(resultDatM_ci,resultDatM2_ci)
        }
      }
      
      qqq_ci[it,] <- ExtResult2(resultDatM_ci,am = am_ci)
    }
    
    qqq_ci <- qqq_ci %>% mutate(
      rIE_b = mQ11 - mQ10,
      rDE_b = mQ10 - mQ00,
      rTE_b = mQ11 - mQ00,
      rPE_b = rIE_b/ rTE_b
    )
    
  }
  
  end_time <- Sys.time()
  elapsed_time <- end_time - start_time
  
  obj <- list()
  
  obj$ori_df <- fitR2$df
  obj$dat_MC <- resultDatM
  class(obj) <- "tvmedg"
  
  ## print result
  cat("Q(a,a):", round(qqq$mQ11, 3),cal_ci(qqq_ci$mQ11,ci,boot = boot),'\n')
  cat("Q(a,a*):", round(qqq$mQ10, 3),cal_ci(qqq_ci$mQ10,ci,boot = boot),'\n')
  cat("Q(a*,a*):", round(qqq$mQ00, 3),cal_ci(qqq_ci$mQ00,ci,boot = boot),'\n')
  
  cat("Indirect:", round(qqq$rIE_b, 3),cal_ci(qqq_ci$rIE_b,ci,boot = boot),'\n')
  cat("Direct:", round(qqq$rDE_b, 3),cal_ci(qqq_ci$rDE_b,ci,boot = boot),'\n')
  cat("Total:", round(qqq$rTE_b, 3),cal_ci(qqq_ci$rTE_b,ci,boot = boot),'\n')
  cat("Proportional explain:", 
      round(qqq$rPE_b, 3),cal_ci(qqq_ci$rPE_b,ci,boot = boot),'\n')
  
  cat("Total time elapsed:",elapsed_time,attr(elapsed_time,"units"),'\n')
  
  invisible(obj)
  
}


library(furrr)
library(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

tvmedg(data = dat0,
       basec = c("age","sex","ow","risk"),
       expo = c("Ap"),
       med = c("Mp"),
       tvar = c("L1","L2","L3"),
       outc = c("Yp"),
       time = c("mm"),
       lag = 2,
       norev = c("Mp"),
       cont_exp = F,
       cont_exp_std = F,
       tvar_to_med = F,
       mreg = "binomial",
       lreg = c("binomial","gaussian","gaussian"),
       yreg = "binomial",
       sp_list = NULL,
       sp_type = NULL,
       sp_df= NULL,
       followup = 12,
       seed = 123,
       montecarlo = 10,
       boot = T, 
       nboot = 2,
       ci = .95,
       parallel=TRUE)

followup = 12

re_df <- fitR2$res_df

### gform 2

fitR2$res_df

g_form2 <- function(data = fitR2$res_df,model = fitR2,followup = 12,am = 0,ay = 1){
  
  norev_var <- model$norev_var
  
  dddd <- data$`1` %>% as.data.frame()
  lagg <- dddd %>% dplyr::select(contains("L1l")) %>% ncol()
  
  
  id     <- dddd$idsim[1]
  id_ori <- dddd$id[1]
  
  # Baseline covariates
  
  Vp <- dddd %>% select(starts_with("v"))
  
  ## create vector for time, outcome
  mm <- seq_len(followup)
  Yp2 <- numeric(followup)
  
  ## create vector for mediator
  Mp <- numeric(followup)
  
  # create vector for time-varying covariates (contribute to mediator models and outcome)
  for (u in 1:length(model$L)) {
    ## create vector with length = followup
    assign(paste0("Lmp", u), numeric(followup))
    assign(paste0("Lp", u), numeric(followup))
    ## assign for followup = 1
    vec <- get(paste0("Lmp", u))
    vec[1] <- dddd |> pull(paste0("L", u))
    assign(paste0("Lmp", u), vec)
    assign(paste0("Lp", u), vec)
  }
  
  
  for (l in lagg:followup) {
  l = 2  
    if (Yp2[l-1]==1) {
      # event occurred at time (l-1) → stop here
      actual_length <- l - 1
      break
    } else{
      
      for (u in 1:length(model$L)){
        for (zz in 1:lagg){
          ## create lag variable
          assign(paste0("i",zz),max(1,l - zz))
          
          
          assign(paste0("Lmp",u,"l", zz), get(paste0("Lp", u))[get(paste0("i",zz))])
          assign(paste0("Lp",u,"l", zz), get(paste0("Lp", u))[get(paste0("i",zz))])
          assign(paste0("Mp","l", zz), get(paste0("Mp"))[get(paste0("i",zz))])

          # A and lag of A
          assign(paste0("A",u,"l",zz),am)
          
          dfMp <- data.frame(Vp, 
                     get(paste0("A",u,"l",zz)), 
                     get(paste0("Mp","l", zz)), 
                     get(paste0("Lmp",u,"l", zz)),
                     j = l)
          
          # var_fm <- attr(data$M[[1]]$terms, "term.labels")
          # var_fm <- gsub(".*\\(([^,]+),.*", "\\1", var_fm)
          # dfMp <- dddd[2,] %>% select(matches(var_fm))%>%
          #   mutate(j = l) %>% as.data.frame()
          # 
          # dfMp[startsWith(colnames(dfMp), "A")] <- am
          # 
          # data.frame(Vp, A1 = am, Al1 = am, Al2 = am, Mpl1, Mpl2, L1l1 = Lmp1l1, L1l2 = Lmp1l1,
          #            L2l1 = Lmp2l1, T2l2 = Lmp2l2, T3l1 = Lmp3l1, T3l2 = Lmp3l2,
          #            j = l)
        }
      }
      
    }
  
  }
  
}

library(furrr)
library(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

foreach(
  data   = fitR2$res_df,
  .combine = rbind,
  .packages = c("splines", "data.table","tidyverse")
) %dopar% {
  rbind(
    g_form(data, model = fitR2, followup = 12, am = 1, ay = 1),
    g_form(data, model = fitR2, followup = 12, am = 0, ay = 1),
    g_form(data, model = fitR2, followup = 12, am = 0, ay = 0)
  )
}

gform_wrapper2 <- function(data, model) {
  outdat11 <- g_form(data=data, model = model,followup = followup, ay = 1, am = 1)
  outdat10 <- g_form(data=data, model = model,followup = followup, ay = 1, am = 0)
  outdat00 <- g_form(data=data, model = model,followup = followup, ay = 0, am = 0)
  bind_rows(outdat11, outdat10, outdat00)
}

future_map_dfr(fitR2$res_df,
               ~ gform_wrapper2(data = .,model = fitR2))


###

library(doParallel)
cl <- makeCluster(4)  # or detectCores()
registerDoParallel(cl)

bench_foreach <- function() {
  foreach(
    data = fitR2$res_df,
    .combine = rbind,
    .packages = c("splines", "data.table", "tidyverse")
  ) %dopar% {
    rbind(
      g_form(data, model = fitR2, followup = 12, am = 1, ay = 1),
      g_form(data, model = fitR2, followup = 12, am = 0, ay = 1),
      g_form(data, model = fitR2, followup = 12, am = 0, ay = 0)
    )
  }
}


library(furrr)
plan(multisession, workers = 4)

gform_wrapper2 <- function(data, model) {
  bind_rows(
    g_form(data=data, model = model, followup = 12, ay = 1, am = 1),
    g_form(data=data, model = model, followup = 12, ay = 1, am = 0),
    g_form(data=data, model = model, followup = 12, ay = 0, am = 0)
  )
}

bench_futuremap <- function() {
  future_map_dfr(fitR2$res_df, ~ gform_wrapper2(.x, model = fitR2))
}


library(bench)

results <- bench::mark(
  foreach = foreach(
    data = fitR2$res_df,
    .combine = rbind,
    .packages = c("splines", "data.table", "tidyverse")
  ) %dopar% {
    rbind(
      g_form(data, model = fitR2, followup = 12, am = 1, ay = 1),
      g_form(data, model = fitR2, followup = 12, am = 0, ay = 1),
      g_form(data, model = fitR2, followup = 12, am = 0, ay = 0)
    )
  },
  futuremap = bench_futuremap(),
  iterations = 3,
  check = FALSE
)

print(results)
