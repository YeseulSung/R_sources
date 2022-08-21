library(glmnet)
library(magrittr)
library(tidyverse)

####################################################################
############################# DATADAT ##############################
####################################################################
#returns statistic (Pearson X2, Likelihood G2, Ordianl M2)
#returns row_prob, col_prob
#returns expected (MUij_hat = pi_i+ * pi_+j * N)
#returns std_resid (Standardized Residual for each cell)

datadat = function(count, I, J, byrow = F, ordinal = T, score_x = NULL, score_y = NULL){
  data = matrix(count, nrow =I, ncol = J, byrow = byrow, dimnames = list(X = score_x, Y = score_y))
  N = sum(data)
  row_prob = matrix(apply(data, 1, sum) / N, nrow = I)
  col_prob = matrix(apply(data, 2, sum) / N, ncol = J)
  
  E = row_prob %*% col_prob * N
  A = data / E
  R = (data - E) / sqrt(E * (1 - row_prob) %*% (1 - col_prob))
  X2 = sum((data - E)^2 / E)
  G2 = 2 * sum(data * log(data / E))
  df_x = df_g = (I - 1) * (J - 1)
  p_x = pchisq(X2, df = df_x, lower.tail = F)
  p_g = pchisq(G2, df = df_g, lower.tail = F)
  if (ordinal == T){
    ui_u = matrix(score_x - sum(row_prob * score_x), nrow = I)
    vi_v = matrix(score_y - sum(col_prob * score_y), ncol = J)
    r = sum(ui_u %*% vi_v * (data / N)) / sqrt(sum(ui_u ^ 2 * row_prob) * sum(vi_v ^ 2 * col_prob))
    M2 = (N - 1) * r^2
    df_m = 1
    p_m = pchisq(M2, df = df_m, lower.tail = F)
    outmat = matrix(c(X2, G2, M2, df_x, df_g, df_m, p_x, p_g, p_m), nrow = 3, byrow = T, 
                    dimnames = list(info = c("statistic", "df", "p_value"),
                                    test = c("Pearson X2", "Likelihood G2", "Ordinal M2")))
  }
  if (ordinal == F){
    outmat = matrix(c(X2, G2, df_x, df_g, p_x, p_g), nrow = 3, byrow = T, 
                    dimnames = list(info = c("statistic", "df", "p_value"),
                                    test = c("Pearson X2", "Likelihood G2")))
  }
  result = list(statistic = outmat, row_prob = row_prob, col_prob = col_prob, expected = E, std_resid = R, association = A)
  return(result)
}



####################################################################
########################### SAMPLE_ODDS ############################
####################################################################

sample_odds = function(count, I, J, byrow = F, ascending_x = T, ascending_y = T, sparse = T){
  data = matrix(count, nrow = I, ncol = J, byrow = byrow)
  if (ascending_x == F) data = data[I:1, ]
  if (ascending_y == F) data = data[, J:1]
  if (I == 2 |J == 2){
    if (I == 2 & J != 2){data = t(data) ; I = nrow(data) ; J = ncol(data)}
    odds = data[, 2] / data[, 1]
    OR = sd = numeric(I - 1)
    for (i in 2:I){
      if (sparse == T){OR[i-1] = odds[i] / odds[1]; sd[i-1] = sqrt(sum(1/data[1, ]) + sum(1/data[i, ]))}
      if (sparse == F){OR[i-1] = odds[i] / odds[i-1]; sd[i-1] = sqrt(sum(1/data[i, ]) + sum(1/data[(i-1), ]))}
    }
    odds = cbind(1:I, odds)
    #OR is in ascending x (None, Light smoker, Heavy smoker
  }
  else{
    odds = NULL
    OR = sd = matrix(0, nrow = I-1, ncol = J-1)
    for (i in 2:I){
      for (j in 2:J){
        if (sparse == T){
          OR[(i-1), (j-1)] = (data[1, 1] * data[i, j]) / (data[1, j] * data[i, 1])
          sd[(i-1), (j-1)] = sqrt(1/data[1, 1] + 1/data[1, j] + 1/data[i, 1] + 1/data[i, j])
        }
        if (sparse == F){
          OR[(i-1), (j-1)] = (data[i, j] * data[(i-1), (j-1)]) / (data[(i-1), j] * data[i, (j-1)])
          sd[(i-1), (j-1)] = sqrt(1/data[(i-1), (j-1)] + 1/data[(i-1), j] + 1/data[i, (j-1)] + 1/data[i, j])
        }
      }
    }
  }
  if (I == 2 | J == 2){CI = cbind(exp(log(OR) - 1.96 * sd), exp(log(OR) + 1.96 * sd))} else{
    CI = list(LB = exp(log(OR) - 1.96 * sd), UB = exp(log(OR) + 1.96 * sd))
  }
  return(list(odds = odds, OR = OR, CI = CI))
}


####################################################################
########################### STOCHASTIC #############################
####################################################################

stochastic = function(count, I, J, byrow = F){
  data = matrix(count, nrow = I, ncol = J, byrow = byrow)
  if (J == 2) data = t(data)
  F1 = cumsum(data[1,] / apply(data, 1, sum)[1])
  F2 = cumsum(data[2,] / apply(data, 1, sum)[2])
  result = ifelse(F1[F1 != F2] >= F2[F1 != F2], 1, 2)
  N = sum(F1 != F2)
  if (sum(result == 1) == N){
    print("2nd group is stochastically higher than 1st group") 
    print("conditional dist is not identical > X n Y NOT INDEP")
  } else if (sum(result == 2) == N){
    print("1st group is stochastically higher than 2nd group")
    print("conditional dist is not identical > X n Y NOT INDEP")
  } else {print("No stochastical relationship between 1st and 2nd groups.")}
  return(result)
}

####################################################################
############################# CONCDISC #############################
####################################################################
concdisc = function(count, I, J, ascending_x = T, ascending_y = T, byrow = F){
  data = matrix(count, nrow = I, ncol = J, byrow = byrow)
  if (ascending_x == F) data = data[I:1, ]; if (ascending_y == F) data = data[, J:1]
  idx_row = 1:I; idx_col = 1:J; N = sum(count)
  C = 0; D = 0
  if (I == 2){
    for (j in 1:(J-1)){C = C + data[1, j] * sum(data[2, idx_col > j])}
    for (j in J:2){D = D + data[1, j] * sum(data[2, idx_col < j])}
  } else if (J == 2){
    for (i in 1:(I-1)){
      C = C + data[i, 1] * sum(data[idx_row > i, 2])
      D = D + data[i, 2] * sum(data[idx_row > i, 1])}
  } else{
    for (i in 1:(I-1)){for (j in 1:(J-1)){
      C = C + data[i, j] * sum(data[idx_row > i, idx_col > j])}}
    for (i in 1:(I-1)){for (j in J:2){
      D = D + data[i, j] * sum(data[idx_row > i, idx_col < j])}}
  }
  prob_con = C / (C + D); prob_dis = D / (C + D)
  gamma = (prob_con - prob_dis) / (prob_con + prob_dis)
  if (round(gamma) == 0){print("may consider indep of X & Y")}
  if (I == 2 | J == 2){
    d = which(dim(data) == 2)
    delta = (C - D) / (apply(data, d, sum)[1] * apply(data, d, sum)[2])
    if (delta < 0){
      print("outcome of larger group tends to be smaller than outcome of smaller grp")
      print("Smaller grp stochastically higher")}
    if (delta > 0){
      print("outcome of larger group tends to be larger than outcome of smaller grp")
      print("Larger grp stochastically higher")}
    return(list(C_D = c(C, D), prob_c_d = c(prob_con, prob_dis), gamma = gamma, delta = delta))}
  return(list(C_D = c(C, D), prob_c_d = c(prob_con, prob_dis), gamma = gamma))
}


####################################################################
############################ SUBTABLES #############################
####################################################################
subtables = function(count, I, J, byrow = F){
  data = matrix(count, nrow = I, ncol = J, byrow = byrow)
  idx = cbind(rep(2:I, rep((J-1), (J-1))), rep(2:J, (I-1)))
  len = dim(idx)[1]
  result = list(statistic = list(), observed = list(), expected = list(), sresid = list())
  for (l in 1:len){
    i = idx[l, 1]; j = idx[l, 2]
    sub = matrix(c(sum(data[1:(i-1), 1:(j-1)]), sum(data[i, 1:(j-1)]), sum(data[1:(i-1), j]), data[i, j]), nrow = 2)
    N = sum(sub)
    row_prob = matrix(apply(sub, 1, sum) / N, nrow = 2)
    col_prob = matrix(apply(sub, 2, sum) / N, ncol = 2)
    
    E = row_prob %*% col_prob * N
    R = (sub - E) / sqrt(E * (1 - row_prob) %*% (1 - col_prob))
    X2 = round(sum((sub - E)^2 / E), digits = 4)
    G2 = round(2 * sum(sub * log(sub / E)), digits = 4)
    df_x = df_g = 1
    p_x = round(pchisq(X2, df = df_x, lower.tail = F), digits = 4)
    p_g = round(pchisq(G2, df = df_g, lower.tail = F), digits = 4)
    result$statistic[[l]] = matrix(c(X2, df_x, p_x, G2, df_g, p_g), nrow = 3, dimnames = list(info = c("stat", "df", "p_value"), test = c("Pearson X2", "Likelihood G2")))
    result$observed[[l]] = sub
    result$expected[[l]] = E
    result$sresid[[l]] = R
  }
  return(result)
}



CMH = function(y, n, z){
  mat = cbind(y, n-y)
  mu = var = n = M = H = numeric(length(unique(z)))
  for (i in 1:length(unique(z))){
    matt = mat[z == unique(z)[i], ]
    row = apply(matt, 1, sum)
    col = apply(matt, 2, sum)
    mu[i] = row[1] * col[1] / sum(matt)
    var[i] = prod(row) * prod(col) / (sum(matt)^2 * (sum(matt) - 1))
    n[i] = matt[1, 1]
    M[i] = matt[1, 1] * matt[2, 2] / sum(matt)
    H[i] = matt[1, 2] * matt[2, 1] / sum(matt)
  }
  CMH = sum(n - mu)^2 / sum(var)
  MH = sum(M) / sum(H)
  return(c("CMH : ", CMH, "MH OR : ", MH))
}