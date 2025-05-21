library(tidyverse)
library(splines)
library(foreach)
library(doParallel)
library(DescTools)
#-------------------------------------------------------------------------------
dat <- readRDS("./Data/tvmed_dat100.RDS")

dat <- dat |> 
  group_by(id) |>
  mutate(
    A_lag1 = lag(Ap, n = 1, default = 0),
    A_lag2 = lag(Ap, n = 2, default = 0),
    M_lag1 = lag(Mp, n = 1, default = 0),
    M_lag2 = lag(Mp, n = 2, default = 0),
    L1_lag1 = lag(L1, n = 1, default = 0),
    L1_lag2 = lag(L1, n = 2, default = 0),
    L2_lag1 = lag(L2, n = 1, default = 100),
    L2_lag2 = lag(L2, n = 2, default = 100),
    L3_lag1 = lag(L3, n = 1, default = 80),
    L3_lag2 = lag(L3, n = 2, default = 80)
  ) |> ungroup()


source("./Rcode/0Gformula.R")


ExtResult <- function(data) {
  Q11 <- data |> filter(lastid == 1 & Ay ==1 & Am ==1)
  Q10 <- data |> filter(lastid == 1 & Ay ==1 & Am ==0)
  Q00 <- data |> filter(lastid == 1 & Ay ==0 & Am ==0)
  
  cat("Q(1,1):", round(mean(Q11$Yp), 3),'\n')
  cat("Q(1,0):", round(mean(Q10$Yp), 3),'\n')
  cat("Q(0,0):", round(mean(Q00$Yp), 3),'\n')
  
  cat("Indirect:", round(mean(Q11$Yp) - mean(Q10$Yp), 3),'\n')
  cat("Direct:", round(mean(Q10$Yp) - mean(Q00$Yp), 3),'\n')
  cat("Total:", round(mean(Q11$Yp) - mean(Q00$Yp), 3),'\n')
  cat("Proportional explain:", 
      round((mean(Q11$Yp) - mean(Q10$Yp))/(mean(Q11$Yp) - mean(Q00$Yp)), 3),'\n')
}


#----- g-formula estimate for single mediator
#===============================================================================
# Parallel

cl <- makeCluster(8)
registerDoParallel(cl)

start_time <- Sys.time()
datMC <- g_boot(data = dat, length = 12, seed = 0, montecarlo = 10000)
end_time <- Sys.time()
(elapsed_time <- end_time - start_time)

# Stop the cluster
stopCluster(cl)

ExtResult(datMC)

# Q(1,1): 0.943 
# Q(1,0): 0.879 
# Q(0,0): 0.019 
# Indirect: 0.064 
# Direct: 0.86 
# Total: 0.923 
# Proportional explain: 0.069 


#----- Bootstrap: single mediator
#===============================================================================
start_time <- Sys.time()

cl <- makeCluster(8)
registerDoParallel(cl)

B <- 5

Q11_b <- numeric(length = B)
Q10_b <- numeric(length = B)
Q00_b <- numeric(length = B)


rIE_b <- numeric(length = B)
rDE_b <- numeric(length = B)
rTE_b <- numeric(length = B)
rPE_b <- numeric(length = B)


for(b in 1:B) {
  tempdat <- g_boot(data = dat, length = 12, seed = b, montecarlo = 5000)
  tQ11 <- tempdat |> filter(lastid == 1 & Ay == 1 & Am == 1)
  tQ10 <- tempdat |> filter(lastid == 1 & Ay == 1 & Am == 0)
  tQ00 <- tempdat |> filter(lastid == 1 & Ay == 0 & Am == 0)
  
  Q11_b[b] <- mean(tQ11$Yp)
  Q10_b[b] <- mean(tQ10$Yp)
  Q00_b[b] <- mean(tQ00$Yp)
  
  rIE_b[b] <- Q11_b[b] - Q10_b[b]
  rDE_b[b] <- Q10_b[b] - Q00_b[b]
  rTE_b[b] <- Q11_b[b] - Q00_b[b]
  rPE_b[b] <- rIE_b[b]/rTE_b[b]
}

stopCluster(cl)

end_time <- Sys.time()
(elapsed_time <- end_time - start_time)

datBoot_repl <- data.frame(id_boot = 1:B, Q11_b, Q10_b, Q00_b, rIE_b, rDE_b, rTE_b, rPE_b) 

quantile(datBoot_repl$Q11_b, probs = c(0.025, 0.975))
quantile(datBoot_repl$Q10_b, probs = c(0.025, 0.975))
quantile(datBoot_repl$Q00_b, probs = c(0.025, 0.975))

quantile(datBoot_repl$rDE_b, probs = c(0.025, 0.975))
quantile(datBoot_repl$rIE_b, probs = c(0.025, 0.975))
quantile(datBoot_repl$rTE_b, probs = c(0.025, 0.975))
quantile(datBoot_repl$rPE_b, probs = c(0.025, 0.975))










