###

dat0 <- readRDS("./Data/tvmed_dat100_12mo.RDS")

process_data <- function(basec,expo,med,tvar,lag,outc,time,
                         norev = NULL,tvar_to_med = FALSE,
                         cont_exp = FALSE,cont_exp_std = F,
                         sp_list = NULL,
                         sp_type = NULL,
                         sp_df= NULL, data){
  
  ## detect which variables are non-reversible
  if(length(which(expo %in% norev)) != 0){
    norev_expo <- paste0("A",which(expo %in% norev))
  } else{
    norev_expo <- NULL
  }
  
  if(length(which(med %in% norev)) != 0){
    norev_med <- paste0("M",which(med %in% norev))
  } else {
    norev_med <- NULL
  }
  
  if(length(which(tvar %in% norev)) != 0){
    norev_tvar <- paste0("L",which(tvar %in% norev))
  } else {
    norev_tvar <- NULL
  }
  
  norev_var <- c(norev_expo,norev_med,norev_tvar)
  
  
  sp_var <- NULL
  
  ## detect which variables in splines list
  ### time-fixed var
  if(length(which(sp_list %in% basec)) != 0){
    sp_fix <- paste0("v",which(basec %in% sp_list))
  } else{
    sp_fix <- NULL
  }
  
  sp_var[which(sp_list %in% basec)] <- sp_fix
  
  ## exposure var
  if(length(which(sp_list %in% expo)) != 0){
    sp_expo <- paste0("A",which(expo %in% sp_list))
  } else{
    sp_expo <- NULL
  }
  
  sp_var[which(sp_list %in% expo)] <- sp_expo
  
  ## mediator var
  if(length(which(sp_list %in% med)) != 0){
    sp_med <- paste0("M",which(med %in% sp_list))
  } else {
    sp_med <- NULL
  }
  
  sp_var[which(sp_list %in% med)] <- sp_med
  
  ## time-varying var
  if(length(which(sp_list %in% tvar)) != 0){
    sp_tvar <- paste0("L",which(tvar %in% sp_list))
  } else {
    sp_tvar <- NULL
  }
  
  sp_var[which(sp_list %in% tvar)] <- sp_tvar
  
  ## time var
  if(length(which(sp_list %in% time)) != 0){
    sp_time <- paste0("j")
  } else{
    sp_time <- NULL
  }
  
  sp_var[which(sp_list %in% time)] <- sp_time
  
  ## continous exposure
  if(cont_exp == F & cont_exp_std == T){
    stop("standardize only applicable when continuous exposure was TRUE, default is FALSE")
  } else if (cont_exp == F & cont_exp_std == F){
    ## binary exposure
    expo_mean <- 0
  } else if (cont_exp == T & cont_exp_std == F){
    ## continuous exposure without standardize
    expo_mean <- mean(data[,expo])
  } else if (cont_exp == T & cont_exp_std == T){
    ## continuous exposure with standardize
    expo_mean <- 0
    data[,expo] <- as.numeric(scale(data[,expo]))
  }
  
  
  ## column names of time-fixed variables
  name_v <- paste0("v",1:length(basec))
  
  out <- data.frame(id = data$id) %>% 
    mutate(data[,basec]) %>% 
    magrittr::set_colnames(c("id",name_v))
  
  ## column names of exposure variables  
  name_e <- paste0("A",1:length(expo))
  
  ## column names of mediator variables  
  name_me <- paste0("M",1:length(med))
  
  ## column names of time-varying variables  
  name_tvar <- paste0("L",1:length(tvar))
  
  
  for (i in 1:lag){
    
    ## column name of lag effect on exposure variable
    co_ex <- paste0(name_e,"l",i)
    co_me <- paste0(name_me,"l",i)
    co_tvar <- paste0(name_tvar,"l",i)
    
    for (k in 1:length(name_e)){
      
      out[,name_e[k]] <- data[,expo[k]]
      name_ep <- {{co_ex}}[k]
      out <- out %>% 
        group_by(id) %>% 
        mutate(
          {{name_ep}} := lag(!!sym(name_e[k]),n=i,default = data[1,expo[k]])
        )
      
    } 
    
    for (k in 1:length(name_me)){
      
      out[,name_me[k]] <- data[,med[k]]
      name_med <- {{co_me}}[k]
      out <- out %>% 
        group_by(id) %>% 
        mutate(
          {{name_med}} := lag(!!sym(name_me[k]),n=i,default = data[1,med[k]])
        ) 
      
    }  
    
    for (k in 1:length(name_tvar)){
      
      out[,name_tvar[k]] <- data[,tvar[k]]
      name_tv <- {{co_tvar}}[k]
      out <- out %>% 
        group_by(id) %>% 
        mutate(
          {{name_tv}} := lag(!!sym(name_tvar[k]),n=i,default = data[1,tvar[k]])
        )
      
    }
    
    
  }
  
  ## outcome variable
  out$Y <- data[,outc]
  
  ## time variable
  out$j <- data[,time]
  
  
  kq <- list()
  kq$df <- out %>% data.frame() 
  kq$norev_var <- norev_var
  kq$am <-  expo_mean
  
  ## column name
  
  eps <- kq$df %>% select(starts_with("A")) %>% colnames()
  tf <- kq$df %>% select(starts_with("v")) %>% colnames()
  tva <- kq$df %>% select(starts_with("L")) %>% colnames()
  mediator <- kq$df %>% select(starts_with("M")) %>% colnames()
  outcome <- kq$df %>% select(starts_with("Y")) %>% colnames()
  timee <- kq$df %>% select(starts_with("j")) %>% colnames()
  
  
  ### exposure variable
  if (any(eps %in% sp_var)){
    ep_sp <- sp_type[which(sp_var %in% eps)]
    ep_df <- sp_df[which(sp_var %in% eps)]
    eps <- paste0("splines::",ep_sp,"(",eps,",df=",ep_df,")")
  } 
  
  ### time-fixed variables
  if (any(tf %in% sp_var)){
    tf_sp <- sp_type[which(sp_var %in% tf)]
    tf_df <- sp_df[which(sp_var %in% tf)]
    tf_v <- keep(tf, ~ any(str_starts(.x, sp_var)))
    tf <- paste0("splines::",tf_sp,"(",tf_v,",df=",tf_df,")")
  } 
  
  ### mediator variables
  if (any(mediator %in% sp_var)){
    mediator_sp <- sp_type[which(sp_var %in% mediator)]
    mediator_df <- sp_df[which(sp_var %in% mediator)]
    mediator_v <- keep(tva, ~ any(str_starts(.x, sp_var)))
    mediator <- paste0("splines::",mediator_sp,"(",mediator_v,",df=",mediator_df,")")
  } 
  
  ### time-varying variables
  if (any(tva %in% sp_var)){
    tva_sp <- sp_type[which(sp_var %in% tva)]
    tva_df <- sp_df[which(sp_var %in% tva)]
    tva_v <- keep(tva, ~ any(str_starts(.x, sp_var)))
    tva <- paste0("splines::",tva_sp,"(",tva_v,",df=",tva_df,")")
  } 
  
  ### time variables
  if (any(timee %in% sp_var)){
    timee_sp <- sp_type[which(sp_var %in% timee)]
    timee_df <- sp_df[which(sp_var %in% timee)]
    timee <- paste0("splines::",timee_sp,"(j,df=",timee_df,")")
  } else {
    timee <- paste0("j")
  }
  
  
  ## formula for M(t)
  l_tm1 <- tva[!tva %in% name_tvar]
  m_tm1 <- mediator[!mediator %in% name_me]
  
  if (tvar_to_med == FALSE){
    formula_mt <- list()
    
    for (i in 1:length(med)){
      formula_mt[i] <- paste(name_me[i],"~",paste(c(eps,m_tm1,l_tm1,tf,timee),collapse = " + "))
    }
    
    kq$fm <- formula_mt
    
    ## formula for L(t)
    formular_lt <- list()
    
    for (i in 1:length(tvar)){
      formular_lt[i] <- paste(name_tvar[i],"~",paste(c(eps,mediator,l_tm1,tf,timee),collapse = " + "))
    }
    
    kq$fl <- formular_lt
    
    ## formula for Y(t)
    formular_y <- paste(outcome,"~",paste(c(eps,mediator,tva,tf,timee),collapse = " + "))
    kq$fy <- formular_y
    
  } else {
    formula_mt <- list()
    
    for (i in 1:length(med)){
      formula_mt[i] <- paste(name_me[i],"~",paste(c(eps,m_tm1,name_tvar,l_tm1,tf,timee),collapse = " + "))
    }
    
    kq$fm <- formula_mt
    
    ## formula for L(t)
    formular_lt <- list()
    
    for (i in 1:length(tvar)){
      formular_lt[i] <- paste(name_tvar[i],"~",paste(c(eps,m_tm1,l_tm1,tf,timee),collapse = " + "))
    }
    
    kq$fl <- formular_lt
    
    ## formula for Y(t)
    formular_y <- paste(outcome,"~",paste(c(eps,tva,mediator,tf,timee),collapse = " + "))
    kq$fy <- formular_y
  }
  
  return(kq)
}


##

resamp <- function(data,boot = FALSE){
  
  df <- as.data.table(data)
  
  # set.seed(seed)
  # cat("Running SEED", seed, "\n")
  # cat("\n")
  # cat("Resampling Data", "\n")
  
  clusters <- unique(df$id)
  samples  <- sample(clusters, length(clusters), replace = TRUE)
  bb       <- table(samples)
  
  #— bootstrap
  if (boot == F) {
    # no bootstrap
    boot_df <- copy(df)
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
    boot_df <- rbindlist(out_list, use.names = TRUE)
  }
  
  boot_df
}

fitg <- function(data,boot = FALSE,
                 mreg = "binomial",
                 lreg = c("binomial","gaussian","gaussian"),
                 yreg = "binomial"){
  
  res_df <- resamp(data = data$df,boot = boot)
  
  fitR <- list()
  
  fitR$df <- res_df %>% as_tibble() 
  
  #----- fit parametric models for
  #--- Mediator models
  
  if(length(mreg) != length(data$fm)){
    stop("the defined regression of M is not equal")
  }
  
  for (i in 1:length(data$fm)){
    fitM <- paste0(data$fm[[i]])
    fitR$M[[i]] <- glm(fitM ,family = mreg[i], data = fitR$df) 
  }
  
  
  #--- Covariate models
  if(length(lreg) != length(data$fl)){
    stop("the defined regression of L is not equal")
  }
  
  for (i in 1:length(data$fl)){
    fitL <- paste0(data$fl[[i]])
    fitR$L[[i]] <- glm(fitL ,family = lreg[i], data = fitR$df) 
  }
  
  #--- Outcome model: 
  fitY <- paste0(data$fy)
  fitR$Y <-  glm(fitY ,family = yreg, data = fitR$df) 
  
  fitR$norev_var <- data$norev_var
  
  fitR$am <- data$am
  
  fitR
}

##
baseline_mc <- function(data = fitR2,montecarlo = 10000){
  
  boot <- as.data.table(data$df)
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
  
  data$res_df <- split(MC, by = "idsim", keep.by = TRUE)
  data
}

##

rFunc <- function(mod, ndat) {
  pred_prob <- predict(mod, newdata = ndat, type = "response")
  return(rbinom(1, size = 1, prob = pred_prob))
}

#### function g_form2


g_form2 <- function(data = fitR2,followup = 12,am,ay){
  
  norev_var <- data$norev_var
  
  dddd <- data$res_df$`1` %>% as.data.frame()
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
  for (u in 1:length(data$L)) {
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
    
    if (Yp2[l-1]==1) {
      # event occurred at time (l-1) → stop here
      actual_length <- l - 1
      break
    } else{
      
      # ## create lag variable
      # # compute 1- and 2-step lags
      # i1  <- max(1, l-1)
      # i2  <- max(1, l-2)
      # 
      # 
      # for (u in 1:length(data$L)){
      #   for (zz in 1:lagg){
      #     assign(paste0("i",zz),max(1,l - zz))
      #     assign(paste0("Lmp",u,"l", zz), get(paste0("Lp", u))[get(paste0("i",zz))])
      #     assign(paste0("Lp",u,"l", zz), get(paste0("Lp", u))[get(paste0("i",zz))])
      #     assign(paste0("Mp","l", zz), get(paste0("Mp"))[get(paste0("i",zz))])
      #     
      #     # A and lag of A
      #     assign(paste0("A",zz),max(1,l - zz))
      #     var_fm <- attr(data$M[[1]]$terms, "term.labels")
      #     var_fm <- gsub(".*\\(([^,]+),.*", "\\1", var_fm)
      #     dfMp <- dddd[2,] %>% select(matches(var_fm))%>% 
      #       mutate(j = l) %>% as.data.frame()
      #     
      #     dfMp[startsWith(colnames(dfMp), "A")] <- am
      #     
      #     data.frame(Vp, A1 = am, Al1 = am, Al2 = am, Mpl1, Mpl2, L1l1 = Lmp1l1, L1l2 = Lmp1l1,
      #                L2l1 = Lmp2l1, T2l2 = Lmp2l2, T3l1 = Lmp3l1, T3l2 = Lmp3l2,
      #                j = l)
      #   }
      # }
    
      
      # Predict mediator
      var_fm <- attr(data$M[[1]]$terms, "term.labels")
      
      var_fm <- gsub(".*\\(([^,]+),.*", "\\1", var_fm)
      
      var_fm <- var_fm[-length(var_fm)] 
      
      dfMp <- dddd %>% select(matches(var_fm))%>% 
        mutate(j = l)
      
      dfMp[startsWith(colnames(dfMp), "A")] <- am
      
      dfMp
    }
  
  
  }
  
}



### test base_line MC
along_redf <- rbindlist(MC_list2) %>% as.data.frame()
test_redf <- rbindlist(fitR2$res_df) %>% as.data.frame()

test_that("baseline_mc function", {
  
  # id, v1 - v4
  expect_equal(along_redf[,1],test_redf[,1])
  expect_equal(along_redf[,2],test_redf[,2])
  expect_equal(along_redf[,3],test_redf[,3])
  expect_equal(along_redf[,4],test_redf[,4])
  expect_equal(along_redf[,5],test_redf[,5])
  
  # exposure variables
  expect_equal(along_redf[,c("A")],test_redf[,c("A1")])
  expect_equal(along_redf[,c("Al1")],test_redf[,c("A1l1")])
  expect_equal(along_redf[,c("Al2")],test_redf[,c("A1l2")])
  
  ## time-dependent covariate
  expect_equal(along_redf[,c("T1")],test_redf[,c("L1")])
  expect_equal(along_redf[,c("T2")],test_redf[,c("L2")])
  expect_equal(along_redf[,c("T3")],test_redf[,c("L3")])
  
  ## time-dependent covariate + lag
  expect_equal(along_redf[,c("T1l1")],test_redf[,c("L1l1")])
  expect_equal(along_redf[,c("T2l1")],test_redf[,c("L2l1")])
  expect_equal(along_redf[,c("T3l1")],test_redf[,c("L3l1")])
  
  expect_equal(along_redf[,c("T1l2")],test_redf[,c("L1l2")])
  expect_equal(along_redf[,c("T2l2")],test_redf[,c("L2l2")])
  expect_equal(along_redf[,c("T3l2")],test_redf[,c("L3l2")])
  
  ## Mediator
  expect_equal(along_redf[,c("M1")],test_redf[,c("M1")])
  expect_equal(along_redf[,c("M1l1")],test_redf[,c("M1l1")])
  expect_equal(along_redf[,c("M1l2")],test_redf[,c("M1l2")])
  
  ## outcome and time
  expect_equal(along_redf[,c("Y")],test_redf[,c("Y")])
  expect_equal(along_redf[,c("j")],test_redf[,c("j")])
})
