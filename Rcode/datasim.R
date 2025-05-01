
# 1. Preparation
#===============================================================================
# Logistic helper
logit_rev <- function(x) {
  1 / (1 + exp(-x))
}

# Truncated normal sampler (upper bound only)
rTruncNorm <- function(n, mean, sd, upper) {
  out <- numeric(n)
  i   <- 1
  while (i <= n) {
    cand <- rnorm(1, mean = mean, sd = sd)
    if (cand <= upper & cand >= 0) {
      out[i] <- cand
      i <- i + 1
    }
  }
  out
}



# 3. Longitudinal Simulation Loop
#===============================================================================

datgeni <- function(id, length) {
  
  # Extract time-fixed covariates
  age  <- rTruncNorm(1, mean = 8, sd = 5, upper = 20) 
  sex  <- rbinom(1, size = 1, prob = 0.47) 
  ow   <- rbinom(1, size = 1, prob = 0.25) 
  risk <- rbinom(1, size = 1, prob = 0.33) 
  
  # Initialize time-varying variables
  
  Ap <- Mp <- L1 <- L2 <- L3 <- Yp <- mm <- numeric()
  
  Ap[1] <- 0
  Mp[1] <- 0
  L1[1] <- 0
  L2[1] <- 100
  L3[1] <- 80
  Yp[1] <- 0
  mm[1] <- 1
  
  for (t in 2:length) {

    if (Yp[t-1]==0) {
      if (t == 2) {
        Apl2 <- Ap[1]
        Mpl2 <- Mp[1]
        L1l2 <- L1[1]
        L2l2 <- L2[1]
        L3l2 <- L3[1]
      } else {
        Apl2 <- Ap[t-2]
        Mpl2 <- Mp[t-2]
        L1l2 <- L1[t-2]
        L2l2 <- L2[t-2]
        L3l2 <- L3[t-2]
      }
      Apl1 <- Ap[t-1]
      Mpl1 <- Mp[t-1]
      L1l1 <- L1[t-1]
      L2l1 <- L2[t-1]
      L3l1 <- L3[t-1]
      
      
      # --- Cardiotoxicity ---
      if (Apl1 == 1) {
        Ap[t] <- 1
      } else {
        lin_A <- ba0 + ba1*Apl1 + ba2*Apl2 + ba3*Mpl1 + ba4*Mpl2 + ba5_1*L1l1 + 
          ba5_2*L2l1 + ba5_3*L3l1 + ba6_1*L1l2 + ba6_2*L2l2 + ba6_3*L3l2 +
          ba7_1*age + ba7_2*sex + ba7_3*ow + ba7_4*risk + s*log(t)
        p_A <- logit_rev(lin_A)
        Ap[t] <- rbinom(1, size = 1, prob = p_A)
      }
      
      # --- Relapse ---
      if (Mpl1 == 1) {
        Mp[t] <- 1
      } else {
        lin_M <- bm0 + bm1*Ap[t] + bm2*Apl1 + bm3*Apl2 + bm4*Mpl1 + bm5*Mpl2 + bm6_1*L1l1 + bm6_2*L2l1 + 
        bm6_3*L3l1 + bm7_1*L1l2 + bm7_2*L2l2 + bm7_3*L3l2 + bm8_1*age + bm8_2*sex + bm8_3*ow + bm8_4*risk + s*log(t)
        p_M <- logit_rev(lin_M)
        Mp[t] <- rbinom(1, size = 1, prob = p_M)
      }
      
      # --- Bloodstream Infection ---
      lin_L1 <- bl1_0 + bl1_1*Ap[t] + bl1_2*Apl1 + bl1_3*Apl2 + bl1_4*Mp[t] + bl1_5*Mpl1 + bl1_6*Mpl2 + 
      bl1_7_1*L1l1 + bl1_7_2*L2l1 + bl1_7_3*L3l1 + bl1_8_1*L1l2 + bl1_8_2*L2l2 + bl1_8_3*L3l2 + 
      bl1_9_1*age + bl1_9_2*sex + bl1_9_3*ow + bl1_9_4*risk + s*log(t)
      p_L1 <- logit_rev(lin_L1)
      L1[t] <- rbinom(1, size = 1, prob = p_L1)
      
      # --- Systolic BP ---
      mu_L2 <- bl2_0 + bl2_1*Ap[t] + bl2_2*Apl1 + bl2_3*Apl2 + bl2_4*Mp[t] + bl2_5*Mpl1 + bl2_6*Mpl2 +  
      bl2_7_1*L1l1 + bl2_7_2*L2l1 + bl2_7_3*L3l1 + bl2_8_1*L1l2 + bl2_8_2*L2l2 + bl2_8_3*L3l2 + 
      bl2_9_1*age + bl2_9_2*sex + bl2_9_3*ow + bl2_9_4*risk + s*log(t)
        
      L2[t] <- rnorm(1, mean = mu_L2, sd = sd_L2)
      
      # --- Diastolic BP ---
      mu_L3 <- bl3_0 + bl3_1*Ap[t] + bl3_2*Apl1 + bl3_3*Apl2 + bl3_4*Mp[t] + bl3_5*Mpl1 + bl3_6*Mpl2 +  
        bl3_7_1*L1l1 + bl3_7_2*L2l1 + bl3_7_3*L3l1 + bl3_8_1*L1l2 + bl3_8_2*L2l2 + bl3_8_3*L3l2 + 
        bl3_9_1*age + bl3_9_2*sex + bl3_9_3*ow + bl3_9_4*risk + s*log(t)
      
      L3[t] <- rnorm(1, mean = mu_L3, sd = sd_L3)
      
      # --- Death ---
      lin_Yp <- by0 + by1*Ap[t] + by2*Apl1 + by3*Apl2 + by4*Mp[t] + by5*Mpl1 + by6*Mpl2 + 
      by7_1*L1[t] + by7_2*L2[t] + by7_3*L3[t] + by8_1*L1l1 + by8_2*L2l1 + by8_3*L3l1 + 
      by9_1*L1l2 + by9_2*L2l2 + by9_3*L3l2 + by10_1*age + by10_2*sex + by10_3*ow + by10_4*risk + s*log(t)
      p_Yp <- logit_rev(lin_Yp)
      Yp[t] <- rbinom(1, size = 1, prob = p_Yp)
      
      # Store data

    } else {
      break
    }
    mm[t] <- t
  }
  gdat <- data.frame(id, mm, Ap, Mp, L1, L2, L3, Yp, age, sex, ow, risk)
  return(gdat)
}



set.seed(12345)


# 2. Define Parameters for Each Variable
#===============================================================================
{
  # Cardiotoxicity
  ba0<- -10; ba1<-5; ba2<-3; ba3<-3; ba4<-1; ba5_1<-0.1; ba5_2<-0.01; ba5_3<-0.01; 
  ba6_1<-0.1; ba6_2<-0.01; ba6_3<-0.01; ba7_1<--0.05; ba7_2<-0.1; ba7_3<-0.1; ba7_4<-10; s<-0.1
  
  # Relapse
  bm0<--20; bm1<-2; bm2<-1; bm3<-0.5; bm4<-5; bm5<-3; bm6_1<-0.1; bm6_2<-0.01; 
  bm6_3<-0.01; bm7_1<-0.1; bm7_2<-0.01; bm7_3<-0.01; bm8_1<--0.05; bm8_2<-0.1; bm8_3<-0.1; bm8_4<-10; s<-0.1
  
  # Bloodstream Infection 
  bl1_0<-log(0.05/0.95); bl1_1<-0.5; bl1_2<-0.3; bl1_3<-0.1; bl1_4<-0.5; bl1_5<-0.3; bl1_6<-0.1; 
  bl1_7_1<-1; bl1_7_2<--0.01; bl1_7_3<--0.01; bl1_8_1<-0.5; bl1_8_2<--0.01; bl1_8_3<--0.01; 
  bl1_9_1<--0.05; bl1_9_2<-0.1; bl1_9_3<-0.1; bl1_9_4<-0.1; s<-0.1
  
  # Systolic Blood Pressure
  bl2_0<-100; bl2_1<--5; bl2_2<--3; bl2_3<--1; bl2_4<--5; bl2_5<--3; bl2_6<--1; 
  bl2_7_1<-1; bl2_7_2<-0.1; bl2_7_3<-0.05; bl2_8_1<-1; bl2_8_2<-0.1; bl2_8_3<-0.05; 
  bl2_9_1<--0.5; bl2_9_2<-1; bl2_9_3<-1; bl2_9_4<-1; s<-0.1; sd_L2 <- 10
  
  # Diastolic Blood Pressure
  bl3_0<-80; bl3_1<--5; bl3_2<--3; bl3_3<--1; bl3_4<--5; bl3_5<--3; bl3_6<--1; 
  bl3_7_1<-1; bl3_7_2<-0.05; bl3_7_3<-0.1; bl3_8_1<-1; bl3_8_2<-0.05; bl3_8_3<-0.1; 
  bl3_9_1<--0.5; bl3_9_2<-1; bl3_9_3<-1; bl3_9_4<-1; s<-0.1; sd_L3 <- 8
  
  # Death
  by0<- -15; by1<-1; by2<-1; by3<-1; by4<-3; by5<-2; by6<-1; 
  by7_1<-0.1; by7_2<-0.01; by7_3<-0.01; by8_1<-0.1; by8_2<-0.01; by8_3<-0.01; 
  by9_1<-0.1; by9_2<-0.01; by9_3<-0.01; by10_1<--0.05; by10_2<-0.1; by10_3<-0.1; by10_4<-0.1; s<-0.1
}



df <- data.frame()

for (i in 1:100) {
  tempdat <- datgeni(id = i, length = 12)
  df <- rbind(df, tempdat)
}

df$lastid <- as.numeric(!duplicated(df$id, fromLast = T))

df_last <- df[df$lastid == 1,]

mean(df_last$Yp)
mean(df_last$Ap)
mean(df_last$Mp)

mean(df_last$Yp[df_last$Ap == 1])
mean(df_last$Yp[df_last$Ap == 0])

mean(df_last$Yp[df_last$Ap == 1 & df_last$Mp == 1])
mean(df_last$Yp[df_last$Ap == 1 & df_last$Mp == 0])
mean(df_last$Yp[df_last$Ap == 0 & df_last$Mp == 0])

summary(glm(Mp ~ Ap, family = binomial(), data = df))
summary(glm(Yp ~ Ap, family = binomial(), data = df))
summary(glm(Yp ~ Mp, family = binomial(), data = df))

# saveRDS(df, "tvmed_dat100_12mo.RDS")






