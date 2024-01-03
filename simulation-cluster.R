# ------------
# simulation 1 Conformal causal inference for cluster-level treatment effect
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
tr_prop <- 0.5 # proportion of data used for training f
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

  oracle_cl_te <- quantile(cl_test_data$trt_effect, 1-alpha/2) - quantile(cl_test_data$trt_effect, alpha/2)

  # Algorithm 1 for marginal cluster-TE -------------
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[index_ca] == 0]
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
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr_A1 <- sample(index_A1, size = floor(tr_prop * 1.5 * length(index_A1)))
  index_tr_A0 <- sample(index_A0, size = floor(tr_prop * 1.5 * length(index_A0)))
  index_tr <- c(index_tr_A1, index_tr_A0)
  index_tr1 <- c(sample(index_tr_A1, size = floor(tr_prop * length(index_tr_A1))), 
                 sample(index_tr_A0, size = floor(tr_prop * length(index_tr_A0))))
  index_tr2 <- setdiff(index_tr, index_tr1)
  index_tr2_A1 <- index_tr2[cl_data$A[index_tr2] == 1]
  index_tr2_A0 <- index_tr2[cl_data$A[index_tr2] == 0]
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  fit_Ybar <- SuperLearner(cl_data$Y[cl_data$cluster_id %in% index_tr1], 
                           X = cl_data[cl_data$cluster_id %in% index_tr1, c("A", "X1", "X2", "N", "C1", "C2")], 
                           family = "gaussian", 
                           SL.library = SLmethods)
  S_A1 <- abs(cl_data$Y[cl_data$cluster_id %in% index_tr2_A1] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_tr2_A1, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  S_A0 <- abs(cl_data$Y[cl_data$cluster_id %in% index_tr2_A0] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_tr2_A0, c("A", "X1", "X2", "N", "C1", "C2")])$pred)
  q_bar_A1 <- quantile(c(S_A1, 1e2), probs = 1-alpha)
  q_bar_A0 <- quantile(c(S_A0, 1e2), probs = 1-alpha)
  pred_A1 <- predict(fit_Ybar, newdata = cl_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=1))$pred
  pred_A0 <- predict(fit_Ybar, newdata = cl_data[, c("A", "X1", "X2", "N", "C1", "C2")] %>% mutate(A=0))$pred
  C.l <- cl_data$A * (cl_data$Y -  pred_A0 - q_bar_A0) + (1-cl_data$A) * (pred_A1 - q_bar_A1 - cl_data$Y)
  C.r <- cl_data$A * (cl_data$Y -  pred_A0 + q_bar_A0) + (1-cl_data$A) * (pred_A1 + q_bar_A1 - cl_data$Y)
  fit_C.l <- SuperLearner(Y = C.l[cl_data$cluster_id %in% index_tr], 
                          X = cl_data[cl_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  fit_C.r <- SuperLearner(Y = C.r[cl_data$cluster_id %in% index_tr], 
                          X = cl_data[cl_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1", "C2")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  pred_ca_C.l <- predict(fit_C.l, newdata = cl_data[cl_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_ca_C.r <- predict(fit_C.r, newdata = cl_data[cl_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1", "C2")])$pred
  S_star <- pmax(pred_ca_C.l - C.l[cl_data$cluster_id %in% index_ca], C.r[cl_data$cluster_id %in% index_ca] - pred_ca_C.r)
  q_star <- quantile(c(S_star, 1e2), probs = 1-gamma)
  
  pred_test_C.l <-  predict(fit_C.l, newdata = cl_test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  pred_test_C.r <-  predict(fit_C.r, newdata = cl_test_data[, c("X1", "X2", "N", "C1", "C2")])$pred
  CI2_nested <- data.frame(ci.l = pred_test_C.l - q_star, 
                           ci.r = pred_test_C.r + q_star)
  coverage2_nested <- (cl_test_data$trt_effect >= CI2_nested$ci.l) & (cl_test_data$trt_effect <= CI2_nested$ci.r)
  cov_prob2_nested <- mean(coverage2_nested)
  avg_length2_nested <- mean(CI2_nested$ci.r -  CI2_nested$ci.l)
  
  

  # Algorithm 1 for local cluster-TE -------------
  cl_data <- dplyr::filter(cl_data, C2 ==1, C1 >= 2)
  cl_data <- cl_data %>% mutate(cluster_id = 1:nrow(cl_data))
  cl_test_data <- dplyr::filter(cl_test_data, C2 ==1, C1 >= 2)
  oracle_cl_te_local <- quantile(cl_test_data$trt_effect, 1-alpha/2) - quantile(cl_test_data$trt_effect, alpha/2)

  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[index_ca] == 0]
  fit_Ybar <- SuperLearner(cl_data$Y[cl_data$cluster_id %in% index_tr], 
                           X = cl_data[cl_data$cluster_id %in% index_tr, c("A", "X1", "X2", "N", "C1")], 
                           family = "gaussian", 
                           SL.library = SLmethods)
  S_A1 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_ca_A1, c("A", "X1", "X2", "N", "C1")])$pred)
  S_A0 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_ca_A0, c("A", "X1", "X2", "N", "C1")])$pred)
  q_bar_A1 <- quantile(c(S_A1, 1e2), probs = 1-alpha)
  q_bar_A0 <- quantile(c(S_A0, 1e2), probs = 1-alpha)
  
  pred_test_A1 <- predict(fit_Ybar, newdata = cl_test_data[, c("A", "X1", "X2", "N", "C1")] %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Ybar, newdata = cl_test_data[, c("A", "X1", "X2", "N", "C1")] %>% mutate(A=0))$pred
  
  CI1_local <- map_dfr(1:nrow(cl_test_data), function(i){
    if(cl_test_data$A[i] == 1){
      data.frame(ci.l = cl_test_data$Y[i] - pred_test_A0[i] - q_bar_A0,
                 ci.r = cl_test_data$Y[i] - pred_test_A0[i] + q_bar_A0)
    } else {
      data.frame(ci.l = pred_test_A1[i] - q_bar_A1 - cl_test_data$Y[i],
                 ci.r = pred_test_A1[i] + q_bar_A1 - cl_test_data$Y[i])
    }
  })
  
  coverage1_local <- (cl_test_data$trt_effect >= CI1_local$ci.l) & (cl_test_data$trt_effect <= CI1_local$ci.r)
  cov_prob1_local <- mean(coverage1_local)
  avg_length1_local <-  mean(CI1_local$ci.r -  CI1_local$ci.l)
  
  # Algorithm 2 for local cluster-TE ------------------
  # a naive approach (using the result from Algorithm 1)
  CI2_naive_local <- data.frame(ci.l = pred_test_A1 - pred_test_A0 - q_bar_A0 - q_bar_A1,
                          ci.r = pred_test_A1 - pred_test_A0 + q_bar_A0 + q_bar_A1)
  coverage2_naive_local <- (cl_test_data$trt_effect >= CI2_naive_local$ci.l) & (cl_test_data$trt_effect <= CI2_naive_local$ci.r)
  cov_prob2_naive_local <- mean(coverage2_naive_local)
  avg_length2_naive_local <- mean(CI2_naive_local$ci.r -  CI2_naive_local$ci.l)
  
  # a nested approach 
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr_A1 <- sample(index_A1, size = floor(tr_prop * 1.5 * length(index_A1)))
  index_tr_A0 <- sample(index_A0, size = floor(tr_prop * 1.5 * length(index_A0)))
  index_tr <- c(index_tr_A1, index_tr_A0)
  index_tr1 <- c(sample(index_tr_A1, size = floor(tr_prop * length(index_tr_A1))), 
                 sample(index_tr_A0, size = floor(tr_prop * length(index_tr_A0))))
  index_tr2 <- setdiff(index_tr, index_tr1)
  index_tr2_A1 <- index_tr2[cl_data$A[index_tr2] == 1]
  index_tr2_A0 <- index_tr2[cl_data$A[index_tr2] == 0]
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  fit_Ybar <- SuperLearner(cl_data$Y[cl_data$cluster_id %in% index_tr1], 
                           X = cl_data[cl_data$cluster_id %in% index_tr1, c("A", "X1", "X2", "N", "C1")], 
                           family = "gaussian", 
                           SL.library = SLmethods)
  S_A1 <- abs(cl_data$Y[cl_data$cluster_id %in% index_tr2_A1] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_tr2_A1, c("A", "X1", "X2", "N", "C1")])$pred)
  S_A0 <- abs(cl_data$Y[cl_data$cluster_id %in% index_tr2_A0] - 
                predict(fit_Ybar, newdata = cl_data[cl_data$cluster_id %in% index_tr2_A0, c("A", "X1", "X2", "N", "C1")])$pred)
  q_bar_A1 <- quantile(c(S_A1, 1e2), probs = 1-alpha)
  q_bar_A0 <- quantile(c(S_A0, 1e2), probs = 1-alpha)
  pred_A1 <- predict(fit_Ybar, newdata = cl_data[, c("A", "X1", "X2", "N", "C1")] %>% mutate(A=1))$pred
  pred_A0 <- predict(fit_Ybar, newdata = cl_data[, c("A", "X1", "X2", "N", "C1")] %>% mutate(A=0))$pred
  C.l <- cl_data$A * (cl_data$Y -  pred_A0 - q_bar_A0) + (1-cl_data$A) * (pred_A1 - q_bar_A1 - cl_data$Y)
  C.r <- cl_data$A * (cl_data$Y -  pred_A0 + q_bar_A0) + (1-cl_data$A) * (pred_A1 + q_bar_A1 - cl_data$Y)
  fit_C.l <- SuperLearner(Y = C.l[cl_data$cluster_id %in% index_tr], 
                          X = cl_data[cl_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  fit_C.r <- SuperLearner(Y = C.r[cl_data$cluster_id %in% index_tr], 
                          X = cl_data[cl_data$cluster_id %in% index_tr, c("X1", "X2", "N", "C1")], 
                          family = "gaussian", 
                          SL.library = SLmethods)
  pred_ca_C.l <- predict(fit_C.l, newdata = cl_data[cl_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1")])$pred
  pred_ca_C.r <- predict(fit_C.r, newdata = cl_data[cl_data$cluster_id %in% index_ca, c("X1", "X2", "N", "C1")])$pred
  S_star <- pmax(pred_ca_C.l - C.l[cl_data$cluster_id %in% index_ca], C.r[cl_data$cluster_id %in% index_ca] - pred_ca_C.r)
  q_star <- quantile(c(S_star, 1e2), probs = 1-gamma)
  
  pred_test_C.l <-  predict(fit_C.l, newdata = cl_test_data[, c("X1", "X2", "N", "C1")])$pred
  pred_test_C.r <-  predict(fit_C.r, newdata = cl_test_data[, c("X1", "X2", "N", "C1")])$pred
  CI2_nested_local <- data.frame(ci.l = pred_test_C.l - q_star, 
                                 ci.r = pred_test_C.r + q_star)
  coverage2_nested_local <- (cl_test_data$trt_effect >= CI2_nested_local$ci.l) & (cl_test_data$trt_effect <= CI2_nested_local$ci.r)
  cov_prob2_nested_local <- mean(coverage2_nested_local)
  avg_length2_nested_local <- mean(CI2_nested_local$ci.r -  CI2_nested_local$ci.l)
  
  
  c(cov_prob1, cov_prob2_naive, cov_prob2_nested, cov_prob1_local, cov_prob2_naive_local, cov_prob2_nested_local,
    avg_length1, avg_length2_naive, avg_length2_nested, avg_length1_local, avg_length2_naive_local, avg_length2_nested_local,
    oracle_cl_te, oracle_cl_te_local)
}
tictoc::toc()
stopCluster(cl)

library(ggplot2)

coverage_cl_marginal <- data.frame(y = as.vector(t(results[1:3,])),
                                   x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p1 <- ggplot(coverage_cl_marginal) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0.7,1)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Coverage probability")

coverage_cl_local <- data.frame(y = as.vector(t(results[4:6,])),
                                x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p2 <- ggplot(coverage_cl_local) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0.7,1)) +
  geom_hline(yintercept = 0.9, linetype = "dashed") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Coverage probability")

avg_length_cl_marginal <- data.frame(y = as.vector(t(results[7:9,])),
                                     x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p3 <- ggplot(avg_length_cl_marginal) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0,10)) +
  geom_hline(yintercept = mean(results[13,]), linetype = "solid") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Length of conformal interval")

avg_length_cl_local <- data.frame(y = as.vector(t(results[10:12,])),
                                  x = factor(rep(c("O", "B-direct", "B-nested"), each = sim_size), levels = c("O", "B-direct", "B-nested")))
p4 <- ggplot(avg_length_cl_local) +
  stat_boxplot(aes(x=x, y=y), geom='errorbar', linetype=1, width=0.2) +
  geom_boxplot(aes(x=x, y=y, fill=x), outlier.shape = NA) + 
  theme_bw() + ylim(c(0,10)) +
  geom_hline(yintercept = mean(results[14,]), linetype = "solid") +
  theme(legend.position="none", text = element_text(size = 16)) + xlab("") +ylab("Length of conformal interval")

library(cowplot)
title1 <- ggdraw() + draw_label("Marginal cluster-level treatment effect", fontface = 'bold')
title2 <- ggdraw() + draw_label("Local cluster-level treatment effect", fontface = 'bold')
p <- plot_grid(title1, title2, p1, p2, p3, p4, ncol = 2, rel_heights = c(0.1,1,1))
save_plot(paste0("sim-cluster-m",m,".png"), p, base_height = 8, base_asp = 1)
