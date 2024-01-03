# ------------
# simulation 3 Conformal causal inference for cluster-level and individual-level treatment effect
# ------------
rm(list = ls())
set.seed(123)
expit <- function(x){1/(1+exp(-x))}
library(tidyverse)
library(SuperLearner) # TMLE
library(foreach)
library(doSNOW)
cl <- makeCluster(4)
registerDoSNOW(cl)

# simulation ----------------
m <- 30 # number of clusters
m_test <- 1000 # number of testing data points
# tr_prop <- 0.5 # proportion of data used for training f
pi <- 0.5 # randomization ratio
alpha <- 0.2 # target non-coverage probability
sim_size <- 1000
package_list <- c("tidyverse", "SuperLearner")
# SLmethods <- c("SL.glm", "SL.randomForest")
SLmethods <- c("SL.glm")


tictoc::tic()
results <- foreach(iter = 1:sim_size, .combine = cbind, .packages = package_list) %dopar% {
  # generate observed data --------------------
  N <- sample(10:50, size = m, replace = T)
  C1 <- rnorm(m, mean = 0.1 * N, sd = 1)
  C2 <- rbinom(m, size = 1, prob = expit(C1/2))
  # A <- rbinom(m, size = 1, prob = pi)
  A <- sample(c(rep(1,m/2),rep(0,m/2)), size = m, replace = F)
  # while(sum(A)>=20 || sum(A) <=10){A <- rbinom(m, size = 1, prob = pi)}
  sim_data <- map_dfr(1:m, function(j) {
    X1 <- rbinom(N[j], size = 1, prob = 0.3 + 0.4 * C2[j])
    X2 <- rnorm(N[j], mean = ifelse(C1[j]>0, 1, -1) * mean(X1), sd = 1)
    Y1 <- sin(C1[j]) * (2 *C2[j] - 1) + abs(X1*X2) + rnorm(N[j], sd = 1) + N[j]/50
    Y0 <- sin(C1[j]) * (2 *C2[j] - 1) + abs(X1*X2) + rnorm(N[j], sd = 1) + rnorm(1, mean = 0, sd = 0.5)
    Y <- A[j] * Y1 + (1-A[j]) * Y0
    cluster_id <- j
    data.frame(cbind(Y, X1, X2, A=A[j], N=N[j], C1=C1[j], C2=C2[j]), cluster_id)
  })
  cl_data <- sim_data %>% group_by(cluster_id) %>% summarise_all(mean)
  
  # generate test data
  N_test <- sample(10:50, size = m_test, replace = T)
  C1_test <- rnorm(m_test, mean = 0.1 * N_test, sd = 1)
  C2_test <- rbinom(m_test, size = 1, prob = expit(C1_test/2))
  test_data <- map_dfr(1:m_test, function(j){
    X1_test <- rbinom(N_test[j], size = 1, prob = 0.3 + 0.4 * C2_test[j])
    X2_test <- rnorm(N_test[j], mean = ifelse(C1_test[j]>0, 1, -1) * mean(X1_test), sd = 1)
    Y1_test <- sin(C1_test[j]) * (2 *C2_test[j] - 1) + abs(X1_test*X2_test) + rnorm(N_test[j], sd = 1) + N_test[j]/50
    Y0_test <- sin(C1_test[j]) * (2 *C2_test[j] - 1) + abs(X1_test*X2_test) + rnorm(N_test[j], sd = 1) + rnorm(1, mean = 0, sd = 0.5)
    # A_test <- rbinom(1, size = 1, prob = expit(-C1_test))
    A_test <- rbinom(1, size = 1, prob = 0.9)
    Y_obs <- A_test * Y1_test + (1-A_test) * Y0_test
    trt_effect <- Y1_test - Y0_test
    cluster_id <- j+m
    data.frame(cbind(trt_effect, Y=Y_obs, X1=X1_test, X2=X2_test, A=A_test, N=N_test[j], C1=C1_test[j], C2=C2_test[j]), cluster_id)
  })
  cl_test_data <- test_data %>% group_by(cluster_id) %>% summarise_all(mean)
  
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_ca_A1 <- sample(index_A1, size = 11)
  index_ca_A0 <- sample(index_A0, size = 11)
  index_tr <- setdiff(cl_data$cluster_id, c(index_ca_A1, index_ca_A0))
  
  # Algorithm 1 for marginal cluster-TE -------------
  fit_Ybar <- SuperLearner(cl_data$Y[cl_data$cluster_id %in% index_tr], 
                           X = cl_data[cl_data$cluster_id %in% index_tr, c("A", "X1", "X2", "N", "C1", "C2")], 
                           family = "gaussian", 
                           SL.library = SLmethods)
  S_A1 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_ca_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_ca_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  q_bar_A1 <- quantile(c(S_A1, 1e2), probs = 1-alpha)
  q_bar_A0 <- quantile(c(S_A0, 1e2), probs = 1-alpha)
  
  pred_test_A1 <- predict(fit_Ybar, newdata = cl_test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Ybar, newdata = cl_test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  
  CI1 <- map_dfr(1:nrow(cl_test_data), function(i){
    if(cl_test_data$A[i] == 1){
      data.frame(ci.l = cl_test_data$Y[i] - pred_test_A0[i] - q_bar_A0,
                 ci.r = cl_test_data$Y[i] - pred_test_A0[i] + q_bar_A0)
    } else {
      data.frame(ci.l = pred_test_A1[i] - q_bar_A1 - cl_test_data$Y[i],
                 ci.r = pred_test_A1[i] + q_bar_A1 - cl_test_data$Y[i])
    }
  })
  
  coverage1 <- (cl_test_data$trt_effect >= CI1$ci.l) & (cl_test_data$trt_effect <= CI1$ci.r)
  cov_prob1 <- mean(coverage1)
  avg_length1 <-  mean(CI1$ci.r -  CI1$ci.l)
  
  # Algorithm 2 for marginal cluster-TE ------------------
  # a naive approach (using the result from Algorithm 1)
  CI2_naive <- data.frame(ci.l = pred_test_A1 - pred_test_A0 - q_bar_A0 - q_bar_A1,
                          ci.r = pred_test_A1 - pred_test_A0 + q_bar_A0 + q_bar_A1)
  coverage2_naive <- (cl_test_data$trt_effect >= CI2_naive$ci.l) & (cl_test_data$trt_effect <= CI2_naive$ci.r)
  cov_prob2_naive <- mean(coverage2_naive)
  avg_length2_naive <- mean(CI2_naive$ci.r -  CI2_naive$ci.l)
  
  # a nested approach 
  # index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  # index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  # index_tr2_A1 <- sample(index_A1, size = 10)
  # index_tr2_A0 <- sample(index_A0, size = 10)
  # index_tr1_A1 <- sample(setdiff(index_A1, index_tr2_A1), size = ceiling(0.5 * (length(index_A1)-10)))
  # index_tr1_A0 <- sample(setdiff(index_A0, index_tr2_A0), size = ceiling(0.5 * (length(index_A0)-10)))
  # index_tr1 <- c(index_tr1_A1, index_tr1_A0)
  # index_ca <- setdiff(cl_data$cluster_id, c(index_tr2_A1, index_tr2_A0, index_tr1))
  
  
  # Algorithm 3 for marginal individual-TE--------------
  fit_Y <- SuperLearner(sim_data$Y[sim_data$cluster_id %in% index_tr], 
                        X = sim_data[sim_data$cluster_id %in% index_tr, c("A", "X1", "X2", "N", "C1", "C2")], 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(sim_data$Y[sim_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_ca_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(sim_data$Y[sim_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_ca_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  
  q_A1 <- data.frame(S = c(S_A1,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca_A1], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_ca_A1], 1)) %>%
    mutate(weights = 1/(length(index_ca_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  q_A0 <- data.frame(S = c(S_A0,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca_A0], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_ca_A0], 1)) %>%
    mutate(weights = 1/(length(index_ca_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  pred_test_A1 <- predict(fit_Y, newdata = test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Y, newdata = test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  
  CI3 <- data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0) + (1-test_data$A)*(pred_test_A1 - q_A1 - test_data$Y),
                    ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0) + (1-test_data$A)*(pred_test_A1 + q_A1 - test_data$Y))
  
  coverage3 <- (test_data$trt_effect >= CI3$ci.l) & (test_data$trt_effect <= CI3$ci.r)
  # coverage3 <- aggregate(coverage3, by = list(test_data$cluster_id), FUN = function(x) {x[1]})$x
  cov_prob3 <- mean(coverage3)
  avg_length3 <-  mean(CI3$ci.r -  CI3$ci.l)
  
  # Algorithm 4 for marginal individual-TE ------------------
  # a naive approach (using the result from Algorithm 3)
  CI4_naive <- data.frame(ci.l = pred_test_A1 - pred_test_A0 - q_A0 - q_A1,
                          ci.r = pred_test_A1 - pred_test_A0 + q_A0 + q_A1)
  coverage4_naive <- (test_data$trt_effect >= CI4_naive$ci.l) & (test_data$trt_effect <= CI4_naive$ci.r)
  cov_prob4_naive <- mean(coverage4_naive)
  avg_length4_naive <- mean(CI4_naive$ci.r -  CI4_naive$ci.l)
  
  # Algorithm 3 for local individual-TE--------------
  ss <- dplyr::filter(sim_data, abs(X2) <= 0.5)
  sim_data <- left_join(ss, aggregate(ss$N, by=list(cluster_id = ss$cluster_id), length), by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
  cl_data <- sim_data %>% group_by(cluster_id) %>% summarise_all(mean)
  tt <- dplyr::filter(test_data, abs(X2) <= 0.5)
  test_data <- left_join(tt, aggregate(tt$N, by=list(cluster_id = tt$cluster_id), length), by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_ca_A1 <- sample(index_A1, size = 10)
  index_ca_A0 <- sample(index_A0, size = 10)
  index_tr <- setdiff(cl_data$cluster_id, c(index_ca_A1, index_ca_A0))
  
  fit_Y <- SuperLearner(sim_data$Y[sim_data$cluster_id %in% index_tr], 
                        X = sim_data[sim_data$cluster_id %in% index_tr, c("A", "X1", "X2", "N", "C1", "C2")], 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(sim_data$Y[sim_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_ca_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(sim_data$Y[sim_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_ca_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  
  q_A1 <- data.frame(S = c(S_A1,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca_A1], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_ca_A1], 1)) %>%
    mutate(weights = 1/(length(index_ca_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  q_A0 <- data.frame(S = c(S_A0,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca_A0], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_ca_A0], 1)) %>%
    mutate(weights = 1/(length(index_ca_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  pred_test_A1 <- predict(fit_Y, newdata = test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Y, newdata = test_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  
  CI3_local <-  data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0) + (1-test_data$A)*(pred_test_A1 - q_A1 - test_data$Y),
                           ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0) + (1-test_data$A)*(pred_test_A1 + q_A1 - test_data$Y))
  
  coverage3_local <- (test_data$trt_effect >= CI3_local$ci.l) & (test_data$trt_effect <= CI3_local$ci.r)
  # coverage3 <- aggregate(coverage3, by = list(test_data$cluster_id), FUN = function(x) {x[1]})
  cov_prob3_local <- mean(coverage3_local)
  avg_length3_local <-  mean(CI3_local$ci.r -  CI3_local$ci.l)
  
  # Algorithm 4 for marginal individual-TE ------------------
  # a naive approach (using the result from Algorithm 3)
  CI4_naive_local <- data.frame(ci.l = pred_test_A1 - pred_test_A0 - q_A0 - q_A1,
                                ci.r = pred_test_A1 - pred_test_A0 + q_A0 + q_A1)
  coverage4_naive_local <- (test_data$trt_effect >= CI4_naive_local$ci.l) & (test_data$trt_effect <= CI4_naive_local$ci.r)
  cov_prob4_naive_local <- mean(coverage4_naive_local)
  avg_length4_naive_local <- mean(CI4_naive_local$ci.r -  CI4_naive_local$ci.l)
  
  c(cov_prob1, cov_prob2_naive, cov_prob3, cov_prob4_naive, cov_prob3_local, cov_prob4_naive_local,
    avg_length1, avg_length2_naive, avg_length3, avg_length4_naive, avg_length3_local, avg_length4_naive_local)
}
tictoc::toc()
stopCluster(cl)


summary_table <- matrix(NA, nrow = 6, ncol = 4)
rownames(summary_table) <- c("O-cl", "X-naive-cl", "O-in", "X-naive-in", "O-in-local","X-naive-in-local")
colnames(summary_table) <- c("coverage_avg", "coverage_sd", "length_avg", "length_sd")
summary_table[,1] <- apply(results[1:6,],1,mean)
summary_table[,3] <- apply(results[7:12,],1,mean)
summary_table[,2] <- apply(results[1:6,],1,sd)
summary_table[,4] <- apply(results[7:12,],1,sd)
round(summary_table,3)

xtable::xtable(summary_table, digits = 3)

