# ------------
# simulation 2 Conformal causal inference for individual-level treatment effect
# ------------
rm(list = ls())
set.seed(123)
library(tidyverse)
library(SuperLearner) # TMLE
library(foreach)
library(doSNOW)
cl <- makeCluster(8)
registerDoSNOW(cl)
expit <- function(x){1/(1+exp(-x))}

# simulation ----------------
m <- 100 # number of clusters
m_test <- 1000 # number of testing data points
tr_prop <- 1/2 # proportion of data used for training f
pi <- 0.5 # randomization ratio
alpha <- 0.1 # target non-coverage probability
gamma <- 0.5 # target non-coverage prob for interval boundaries
sim_size <- 1000
package_list <- c("tidyverse", "SuperLearner")
# SLmethods <- c("SL.glm", "SL.rpart", "SL.nnet")
SLmethods <- c("SL.glm", "SL.randomForest")

tictoc::tic()
results <- foreach(iter = 1:sim_size, .combine = cbind, .packages = package_list) %dopar% {
  # generate observed data --------------------
  N <- sample(10:50, size = m, replace = T)
  C1 <- rnorm(m, mean = 0.1 * N, sd = 1)
  C2 <- rbinom(m, size = 1, prob = expit(C1/2))
  A <- rbinom(m, size = 1, prob = pi)
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
  
  oracle_in_te <- quantile(test_data$trt_effect, 1-alpha/2) - quantile(test_data$trt_effect, alpha/2)
  
  # Algorithm 3 for marginal individual-TE--------------
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[index_ca] == 0]
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
  # coverage3 <- aggregate(coverage3, by = list(test_data$cluster_id), FUN = function(x) {x[1]})
  cov_prob3 <- mean(coverage3)
  avg_length3 <-  mean(CI3$ci.r -  CI3$ci.l)
  
  # Algorithm 4 for marginal individual-TE ------------------
  # a naive approach (using the result from Algorithm 3)
  CI4_naive <- data.frame(ci.l = pred_test_A1 - pred_test_A0 - q_A0 - q_A1,
                          ci.r = pred_test_A1 - pred_test_A0 + q_A0 + q_A1)
  coverage4_naive <- (test_data$trt_effect >= CI4_naive$ci.l) & (test_data$trt_effect <= CI4_naive$ci.r)
  cov_prob4_naive <- mean(coverage4_naive)
  avg_length4_naive <- mean(CI4_naive$ci.r -  CI4_naive$ci.l)
  
  # a nested approach 
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr_A1 <- sample(index_A1, size = floor(tr_prop * length(index_A1)))
  index_tr_A0 <- sample(index_A0, size = floor(tr_prop * length(index_A0)))
  index_tr2_A1 <- sample(index_tr_A1, size = floor(tr_prop * length(index_tr_A1)))
  index_tr2_A0 <- sample(index_tr_A0, size = floor(tr_prop * length(index_tr_A0)))
  index_tr <- c(index_tr_A1, index_tr_A0)
  index_tr2 <- c(index_tr2_A1, index_tr2_A0)
  index_tr1 <- setdiff(index_tr, index_tr2)
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  fit_Y <- SuperLearner(sim_data$Y[sim_data$cluster_id %in% index_tr1], 
                        X = sim_data[sim_data$cluster_id %in% index_tr1, c("A", "X1", "X2", "N", "C1", "C2")], 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(sim_data$Y[sim_data$cluster_id %in% index_tr2_A1] - 
                 predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_tr2_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(sim_data$Y[sim_data$cluster_id %in% index_tr2_A0] - 
                 predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_tr2_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  
  q_A1 <- data.frame(S = c(S_A1,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_tr2_A1], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_tr2_A1], 1)) %>%
    mutate(weights = 1/(length(index_tr2_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  q_A0 <- data.frame(S = c(S_A0,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_tr2_A0], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_tr2_A0], 1)) %>%
    mutate(weights = 1/(length(index_tr2_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  pred_A1 <- predict(fit_Y, newdata = sim_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_A0 <- predict(fit_Y, newdata = sim_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  C.l <- sim_data$A * (sim_data$Y -  pred_A0 - q_A0) + (1-sim_data$A) * (pred_A1 - q_A1 - sim_data$Y)
  C.r <- sim_data$A * (sim_data$Y -  pred_A0 + q_A0) + (1-sim_data$A) * (pred_A1 + q_A1 - sim_data$Y)
  fit_C.l <- SuperLearner(Y = C.l[sim_data$cluster_id %in% index_tr], 
                          X = sim_data[sim_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  fit_C.r <- SuperLearner(Y = C.r[sim_data$cluster_id %in% index_tr], 
                          X = sim_data[sim_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  pred_ca_C.l <- predict(fit_C.l, newdata = sim_data[sim_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_ca_C.r <- predict(fit_C.r, newdata = sim_data[sim_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  S_star <- pmax(pred_ca_C.l - C.l[sim_data$cluster_id %in% index_ca], C.r[sim_data$cluster_id %in% index_ca] - pred_ca_C.r)
  q_star <- data.frame(S = c(S_star,1e2), 
                       cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca], m+1),
                       N = c(sim_data$N[sim_data$cluster_id %in% index_ca], 1)) %>%
    mutate(weights = 1/(length(index_ca)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-gamma) %>%
    .[1,1]
  
  pred_test_C.l <-  predict(fit_C.l, newdata = test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_test_C.r <-  predict(fit_C.r, newdata = test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  CI4_nested <- data.frame(ci.l = pred_test_C.l - q_star, 
                           ci.r = pred_test_C.r + q_star)
  coverage4_nested <- (test_data$trt_effect >= CI4_nested$ci.l) & (test_data$trt_effect <= CI4_nested$ci.r)
  cov_prob4_nested <- mean(coverage4_nested)
  avg_length4_nested <- mean(CI4_nested$ci.r -  CI4_nested$ci.l)
  
  # Algorithm 3 for local individual-TE--------------
  ss <- dplyr::filter(sim_data, abs(X2) <= 0.5)
  sim_data <- left_join(ss, aggregate(ss$N, by=list(cluster_id = ss$cluster_id), length), by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
  cl_data <- sim_data %>% group_by(cluster_id) %>% summarise_all(mean)
  tt <- dplyr::filter(test_data, abs(X2) <= 0.5)
  test_data <- left_join(tt, aggregate(tt$N, by=list(cluster_id = tt$cluster_id), length), by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
  oracle_in_te_local <- quantile(test_data$trt_effect, 1-alpha/2) - quantile(test_data$trt_effect, alpha/2)
  
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[index_ca] == 0]
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
  
  CI3_local <- data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0) + (1-test_data$A)*(pred_test_A1 - q_A1 - test_data$Y),
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
  
  # a nested approach 
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr_A1 <- sample(index_A1, size = floor(tr_prop * length(index_A1)))
  index_tr_A0 <- sample(index_A0, size = floor(tr_prop * length(index_A0)))
  index_tr2_A1 <- sample(index_tr_A1, size = floor(tr_prop * length(index_tr_A1)))
  index_tr2_A0 <- sample(index_tr_A0, size = floor(tr_prop * length(index_tr_A0)))
  index_tr <- c(index_tr_A1, index_tr_A0)
  index_tr2 <- c(index_tr2_A1, index_tr2_A0)
  index_tr1 <- setdiff(index_tr, index_tr2)
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  fit_Y <- SuperLearner(sim_data$Y[sim_data$cluster_id %in% index_tr1], 
                        X = sim_data[sim_data$cluster_id %in% index_tr1, c("A", "X1", "X2", "N", "C1", "C2")], 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(sim_data$Y[sim_data$cluster_id %in% index_tr2_A1] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_tr2_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(sim_data$Y[sim_data$cluster_id %in% index_tr2_A0] - 
                predict(fit_Y, newdata = sim_data[sim_data$cluster_id %in% index_tr2_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  
  q_A1 <- data.frame(S = c(S_A1,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_tr2_A1], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_tr2_A1], 1)) %>%
    mutate(weights = 1/(length(index_tr2_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  q_A0 <- data.frame(S = c(S_A0,1e2),
                     cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_tr2_A0], m+1),
                     N = c(sim_data$N[sim_data$cluster_id %in% index_tr2_A0], 1)) %>%
    mutate(weights = 1/(length(index_tr2_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-alpha) %>%
    .[1,1]
  
  pred_A1 <- predict(fit_Y, newdata = sim_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_A0 <- predict(fit_Y, newdata = sim_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  C.l <- sim_data$A * (sim_data$Y -  pred_A0 - q_A0) + (1-sim_data$A) * (pred_A1 - q_A1 - sim_data$Y)
  C.r <- sim_data$A * (sim_data$Y -  pred_A0 + q_A0) + (1-sim_data$A) * (pred_A1 + q_A1 - sim_data$Y)
  fit_C.l <- SuperLearner(Y = C.l[sim_data$cluster_id %in% index_tr], 
                          X = sim_data[sim_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  fit_C.r <- SuperLearner(Y = C.r[sim_data$cluster_id %in% index_tr], 
                          X = sim_data[sim_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  pred_ca_C.l <- predict(fit_C.l, newdata = sim_data[sim_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_ca_C.r <- predict(fit_C.r, newdata = sim_data[sim_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  S_star <- pmax(pred_ca_C.l - C.l[sim_data$cluster_id %in% index_ca], C.r[sim_data$cluster_id %in% index_ca] - pred_ca_C.r)
  q_star <- data.frame(S = c(S_star,1e2), 
                       cluster_id = c(sim_data$cluster_id[sim_data$cluster_id %in% index_ca], m+1),
                       N = c(sim_data$N[sim_data$cluster_id %in% index_ca], 1)) %>%
    mutate(weights = 1/(length(index_ca)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights)) %>%
    filter(cumsum_weights >= 1-gamma) %>%
    .[1,1]
  
  pred_test_C.l <-  predict(fit_C.l, newdata = test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_test_C.r <-  predict(fit_C.r, newdata = test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  CI4_nested_local <- data.frame(ci.l = pred_test_C.l - q_star, 
                           ci.r = pred_test_C.r + q_star)
  coverage4_nested_local <- (test_data$trt_effect >= CI4_nested_local$ci.l) & (test_data$trt_effect <= CI4_nested_local$ci.r)
  cov_prob4_nested_local <- mean(coverage4_nested_local)
  avg_length4_nested_local <- mean(CI4_nested_local$ci.r -  CI4_nested_local$ci.l)
  
  
  c(cov_prob3, cov_prob4_naive, cov_prob4_nested, cov_prob3_local, cov_prob4_naive_local, cov_prob4_nested_local,
    avg_length3, avg_length4_naive, avg_length4_nested, avg_length3_local, avg_length4_naive_local, avg_length4_nested_local,
    oracle_in_te, oracle_in_te_local)

}
tictoc::toc()
stopCluster(cl)

coverage_in_marginal <- data.frame(y = as.vector(t(results[1:3,])),
                                   x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p1 <- ggplot(coverage_in_marginal) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0.7,1)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Coverage probability")

coverage_in_local <- data.frame(y = as.vector(t(results[4:6,])),
                                x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p2 <- ggplot(coverage_in_local) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0.7,1)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Coverage probability")

avg_length_in_marginal <- data.frame(y = as.vector(t(results[7:9,])),
                                     x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p3 <- ggplot(avg_length_in_marginal) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0,13)) +
  geom_hline(yintercept = mean(results[13,]), linetype = "solid") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Length of conformal interval")

avg_length_in_local <- data.frame(y = as.vector(t(results[10:12,])),
                                  x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p4 <- ggplot(avg_length_in_local) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0,13)) +
  geom_hline(yintercept = mean(results[14,]), linetype = "solid") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Length of conformal interval")

library(cowplot)
title1 <- ggdraw() + draw_label("Marginal individual-level treatment effect", fontface = 'bold')
title2 <- ggdraw() + draw_label("Local individual-level treatment effect", fontface = 'bold')
p <- plot_grid(title1, title2, p1, p2, p3, p4, ncol = 2, rel_heights = c(0.1,1,1))
save_plot(paste0("sim-individual-m",m,".png"), p, base_height = 8, base_asp = 1)