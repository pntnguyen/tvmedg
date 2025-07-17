
df_prep <- function(data) {
  data.frame(
    id   = data$id,
    v1   = data$age,
    v2   = data$sex,
    v3   = data$ow,
    v4   = data$risk,
    A    = data$Ap,
    Al1  = data$A_lag1,
    Al2  = data$A_lag2,
    T1   = data$L1,
    T1l1 = data$L1_lag1,
    T1l2 = data$L1_lag2,
    T2   = data$L2,
    T2l1 = data$L2_lag1,
    T2l2 = data$L2_lag2,
    T3   = data$L3,
    T3l1 = data$L3_lag1,
    T3l2 = data$L3_lag2,
    M1   = data$Mp,
    M1l1 = data$M_lag1,
    M1l2 = data$M_lag2,
    Y    = data$Yp,
    j    = data$mm
  )
}

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
