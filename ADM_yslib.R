library(MASS)
library(magrittr)
library(tidyverse)

##############################################
################### SUBSET ###################
##############################################
subset = function(method, x = NULL, y = NULL, data = NULL, target_idx= NULL){
  #DATA
  if (!is.null(x) & !is.null(y)){
    n = dim(x)[1]; p = dim(x)[2]
    x = as.matrix(x) ; y = as.matrix(y)
  }else if (!is.null(data)){
    n = dim(data)[1]; p = dim(data)[2] - 1
    target_idx = ifelse(is.null(target_idx), p + 1, target_idx)
    x = as.matrix(data[, -target_idx]) ; y = as.matrix(data[, target_idx])
  }else return("WARNINGS : requires pair of x & y or a data as input")
  
  result = vector("list", p + 1) #OUTPUT for p+1 models
  if (method == "best"){
    beta.hat = solve(t(as.matrix(rep(1, n))) %*% as.matrix(rep(1, n))) %*% t(as.matrix(rep(1, n))) %*% y
    y.hat = as.matrix(rep(1, n)) %*% beta.hat
    RSS = t(y - y.hat) %*% (y - y.hat)
    result[[1]] = list(idx = NULL, RSS = RSS)
    for (i in 1:p){
      idx = t(utils::combn(1:p, i)) #all pCi combination matrix
      RSS = numeric(nrow(idx))
      for (j in 1:nrow(idx)){
        mat = as.matrix(cbind(rep(1, n), x[, idx[j, ]]))
        beta.hat = solve(t(mat) %*% mat) %*% t(mat) %*% y
        y.hat = mat %*% beta.hat
        RSS[j] = t(y - y.hat) %*% (y - y.hat)
      }
      result[[i + 1]] = list(idx = idx[which.min(RSS), ], RSS = min(RSS))
    }
  }  else if (method == "forward"){
    beta.hat = solve(t(as.matrix(rep(1, n))) %*% as.matrix(rep(1, n))) %*% t(as.matrix(rep(1, n))) %*% y
    y.hat = as.matrix(rep(1, n)) %*% beta.hat
    RSS = t(y - y.hat) %*% (y - y.hat)
    result[[1]] = list(idx = NULL, RSS = RSS)
    idx = 1:p
    m_idx = NULL
    for (i in 1:p){
      RSS = numeric(p - i + 1)
      for (j in length(RSS)){
        new = c(m_idx, idx[!(idx %in% m_idx)][j])
        mat = as.matrix(cbind(rep(1, n), x[, new]))
        beta.hat = solve(t(mat) %*% mat) %*% t(mat) %*% y
        y.hat = mat %*% beta.hat
        RSS[j] = t(y - y.hat) %*% (y - y.hat)
      }
      m_idx = c(m_idx, idx[!idx %in% m_idx][which.min(RSS)])
      result[[i + 1]] = list(idx = m_idx, RSS = min(RSS))
    }
  } else if (method == "backward"){
    mat = cbind(rep(1, n), x)
    y.hat = mat %*% solve(t(mat) %*% mat) %*% t(mat) %*% y
    idx = 1:p
    result[[1]] = list(idx = idx, RSS = t(y - y.hat) %*% (y - y.hat))
    for (i in 1:p){
      RSS = numeric(p - i + 1)
      for (j in length(RSS)){
        new = idx[idx != idx[j]]
        mat = as.matrix(cbind(rep(1, n), x[,new]))
        beta.hat = solve(t(mat) %*% mat) %*% t(mat) %*% y
        y.hat = mat %*% beta.hat
        RSS[j] = t(y - y.hat) %*% (y - y.hat)
      }
      idx = idx[idx != idx[which.min(RSS)]]
      result[[i + 1]] = list(idx = idx, RSS = min(RSS))
    }
  }
  num = 1:n
  rand_idx = NULL
  comp = matrix(0, nrow = 5, ncol = p+1)
  for (i in 1:5){
    num = num[!num %in% rand_idx]
    rand_idx = sample(num, round(n / 5))
    x_train = x[-rand_idx, ] ; y_train = y[-rand_idx]
    x_test = x[rand_idx, ] ; y_test = y[rand_idx]
    yhat = matrix(rep(1, nrow(x_test))) %*% solve(t(matrix(rep(1, nrow(x_train)))) %*% matrix(rep(1, nrow(x_train)))) %*% t(matrix(rep(1, nrow(x_train)))) %*% y_train
    comp[i, 1] = mean((y_test - yhat)^2)
    for (j in 1:p){
      mat = as.matrix(cbind(rep(1, nrow(x_train)), x_train[, result[[j + 1]]$idx]))
      testmat = as.matrix(cbind(rep(1, nrow(x_test)), x_test[, result[[j + 1]]$idx]))
      comp[i, j + 1] = mean((y_test - testmat %*% solve(t(mat) %*% mat) %*% t(mat) %*% y_train)^2)
    }
  }
  best_model = result[[which.min(apply(comp, 2, mean))]]
  return(list(best_model = best_model, result = result))
}

##############################################
#################### RIDGE ###################
##############################################
ridge = function(lambda, n_train = 200, rand_idx = NULL, target_idx = NULL, x = NULL, y = NULL, data = NULL){
  #DATA
  if (!is.null(x) & !is.null(y)){
    n = dim(x)[1]; p = dim(x)[2]
    x = as.matrix(x)
  }else if (!is.null(data)){
    n = dim(data)[1]; p = dim(data)[2] - 1
    target_idx = ifelse(is.null(target_idx), p + 1, target_idx)
    x = as.matrix(data[, -target_idx]); y = data[, target_idx]
  }else return("WARNINGS : requires pair of x & y or a data as input")
  idx = 1:n
  if (is.null(rand_idx)){
    train_idx = sample(idx, n_train, replace = F)
  }else {train_idx = rand_idx}
  test_idx = idx[!(idx %in% train_idx)]
  train_x = x[train_idx, ]; test_x = x[test_idx, ]
  train_y = as.matrix(y[train_idx]); test_y = as.matrix(y[test_idx])
  mean_trx = apply(train_x, 2, mean); sd_trx = apply(train_x, 2, sd)
  train_x = (train_x - matrix(mean_trx, nrow = nrow(train_x), ncol = ncol(train_x), byrow = T)) / matrix(sd_trx, nrow = nrow(train_x), ncol = ncol(train_x), byrow = T)
  test_x = (test_x - matrix(mean_trx, nrow = nrow(test_x), ncol = ncol(test_x), byrow = T)) / matrix(sd_trx, nrow = nrow(test_x), ncol = ncol(test_x), byrow = T)
  train_y = train_y - mean(train_y)
  test_y = test_y - mean(train_y)
  beta.hat = solve(t(train_x) %*% train_x + lambda * diag(p)) %*% t(train_x) %*% train_y
  y.hat = train_x %*% beta.hat
  train_MSE = mean((train_y - y.hat)^2)
  y.hat0 = test_x %*% beta.hat
  test_MSE = mean((test_y - y.hat0)^2)
  result = list(beta.hat = beta.hat, fitted_value = y.hat0, MSE = test_MSE, training_error = train_MSE)
  return(result)
}