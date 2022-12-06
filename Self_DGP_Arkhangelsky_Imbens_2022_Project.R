rm(list = ls())

# Loading required libraries
library(mvtnorm)
library(pryr)
library(kableExtra)



est_function <- function(W,Y,Z,K){
  
  ## This function takes assignment, data, aggregate shocks and number of clusters K
  
  n <- dim(W)[1]
  T <- dim(W)[2]
  
  # Sufficient cluster-specific statistic creation
  S <- W%*%Z/T 
  # set.seed(1234)
  clust_res <- kmeans(S,K, iter.max = 20)
  clust_ind <- clust_res$cluster
  results_full <- matrix(0,ncol = 2, nrow = K)
  
  for(k in 1:K){
    
    index_k <- clust_ind == k
    if (sum(index_k) > 1){
      Y_k <- Y[index_k,]
      W_k <- W[index_k,]
      n_k <- sum(index_k)
      weights_k_unnorm <- W_k - outer(rep(1,n_k),colMeans(W_k)) - outer(rowMeans(W_k),rep(1,T)) + mean(W_k)
      weights_k <- weights_k_unnorm/(norm(weights_k_unnorm,'f')^2)
      tau_k <- sum(Y_k*weights_k)
      results_full[k,] <- c(tau_k,n_k)
    }
    else {
      results_full[k,] <- c(NA,NA)
    }
  }
  
  agg_tau <- weighted.mean(results_full[,1],results_full[,2],na.rm = TRUE)
  
  
  return(agg_tau)
}

# reproducibility
set.seed(1212)

# Data dimensions
n <- 200
T <- 40
unobserved_types <- 5
unobserved_effect_on_W <- 0.01
# unobserved confounder U_i
U_i <- round(sample(1:unobserved_types, n, replace = T), 0)
U <- rep(U_i, T)
U_mat <- matrix(U, nrow = n, ncol = T)
# Assignment
W <- round(runif(n*T) + rep(U_i, T) * unobserved_effect_on_W, 0)
W_mat <- matrix(W, nrow = n, ncol = T)
# Constant for S_i creation
Z <- cbind(rep(1, T))
# S <- W%*%Z/T
# S_i <- rowSums(matrix(S, nrow = n, ncol = T))/T
# year-time time effect
lambda_t <- rep(1:T, n) + rnorm(n*T)
# random errors
E_y <- rnorm(n*T)
E_mat <- matrix(E_y, nrow = n, ncol = T)
# unit-fixed effect (arbitrary function choice of U_i)
alpha_i <- U_i^2/100
# treatment effect
tau <- 1
# dependent variable
Y <- rep(alpha_i, T) + tau * W + E_y
Y_mat <- matrix(Y, nrow = n, ncol = T)


cor(U, Y)
cor(U, W)


## Number of clusters
K_5 <- floor(n/12)
K_6 <- floor(n/16)
K_7 <- floor(n/24)
K_8 <- floor(n/32)
K_9 <- floor(n/75)
K_10 <- 1
# K_vec <- c(K_1, K_2, K_3, K_4, K_5, K_6, K_7, K_8)
K_vec <- c(K_5, K_6, K_7, K_8, K_9, K_10)

# number of simulations
B <- 10000




# Simulations

# Design 1 Only - No Over-Selection



## Design 1 - No Selection Bias

tau_calc <- pryr::partial(est_function, W=W_mat,Y=Y_mat,Z=Z)

sim_res_1 <- do.call(rbind,lapply(1:B, function(b) {
  
  
  # tau_calc <- pryr::partial(est_function, W=W_b,Y=Y_b,Z=Z)
  a_tau_coll <- sapply(K_vec, FUN=tau_calc)
  
  return(a_tau_coll)
}))

predicted_treatmens <- colMeans(sim_res_1)
rmse_1 <- sqrt(colMeans((sim_res_1-tau)^2))
bias_1 <- colMeans(sim_res_1-tau)

# rmse_1
# bias_1




K_names <- paste0(K_vec, "cluster/s")
table_1 <- round(rbind(predicted_treatmens, rmse_1, bias_1),6)
colnames(table_1) <- c(K_names)
rownames(table_1) <- c('treatment_prediction','rmse_design_1','bias_design_1')

table_1 %>%
  kbl() %>%
  kable_styling()


# Plots
## Predictions
plot(K_vec, predicted_treatmens, col=("black"), main = "Prediction Results", pch = 1, type = "l", ylim = c(0, 1.2))
## RMSE
plot(K_vec, rmse_1, col=("black"), main = "RMSE Results", pch = 1, type = "l")
## Bias
plot(K_vec, bias_1, col=("black"), main = "Bias Results", pch = 1, type = "l")
