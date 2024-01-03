rm(list = ls())
set.seed(123)
library(tidyverse)
library(SuperLearner)
SLmethods <- c("SL.glm", "SL.randomForest")

d <- haven::read_sas("data-application/PPACT/ppact_public_bpi.sas7bdat")
dd <- haven::read_sas("data-application/PPACT/ppact_public_use.sas7bdat")

observed_data <- dplyr::select(d, SID, CLUST, INTERVENTION, PEGS, TIMEPOINT,
                               AGE, FEMALE, disable, Current_Smoke, BMI, Alcohol_Abuse, Drug_Abuse, comorbid, Dep_OR_Anx,
                               pain_count, BL_avg_daily, COT_flag) %>%
  filter(TIMEPOINT %in% c(0,12)) %>%
  pivot_wider(names_from = TIMEPOINT, names_glue = "Y{TIMEPOINT}", values_from = PEGS) %>%
  dplyr::select(-SID) 
observed_data <- observed_data[complete.cases(observed_data),]
# observed_data$CLUST <- as.factor(observed_data$CLUST)
colnames(observed_data)[1] <- "cluster_id"
colnames(observed_data)[2] <- "A"
colnames(observed_data)[16] <- "Y"
observed_data_all <- observed_data %>% left_join(group_by(observed_data, cluster_id) %>% summarise(N = n()), by = "cluster_id") %>% as.data.frame()
cl_data_all <- observed_data_all %>% group_by(cluster_id) %>% summarise_all(mean) %>% as.data.frame
alpha <- c(0.1, 0.2, 0.3, 0.4, 0.5)
tr_prop <- 1/2 # proportion of data used for training f
sim_size <- 100



marginal_results <- map(1: sim_size, function(j){
  
  index_test <- sample(cl_data_all$cluster_id, size = 20, replace = F)
  test_data <- filter(observed_data_all, cluster_id %in% index_test)
  cl_test_data <- filter(cl_data_all, cluster_id %in% index_test)
  observed_data <- filter(observed_data_all, !(cluster_id %in% index_test))
  cl_data <- filter(cl_data_all, !(cluster_id %in% index_test))
  
  m <- length(unique(cl_data$cluster_id))
  
  # marginal cluster-level effects -------
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 0]
  fit_Ybar <- SuperLearner(cl_data$Y[cl_data$cluster_id %in% index_tr], 
                           X = dplyr::select(cl_data[cl_data$cluster_id %in% index_tr,], -c("Y", "cluster_id")), 
                           family = "gaussian", 
                           SL.library = SLmethods)
  S_A1 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Ybar, newdata = dplyr::select(cl_data[cl_data$cluster_id %in% index_ca_A1,], -c("Y", "cluster_id")))$pred)
  S_A0 <- abs(cl_data$Y[cl_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Ybar, newdata = dplyr::select(cl_data[cl_data$cluster_id %in% index_ca_A0,], -c("Y", "cluster_id")))$pred)
  q_bar_A1 <- quantile(c(S_A1, 1e2), probs = 1-alpha)
  q_bar_A0 <- quantile(c(S_A0, 1e2), probs = 1-alpha)
  
  pred_test_A1 <- predict(fit_Ybar, newdata = dplyr::select(cl_test_data, -c("Y", "cluster_id")) %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Ybar, newdata = dplyr::select(cl_test_data, -c("Y", "cluster_id")) %>% mutate(A=0))$pred
  
  cl_marginal <- map_dfr(1:length(alpha), function(j){
    CI1 <- data.frame(ci.l = cl_test_data$A*(cl_test_data$Y - pred_test_A0 - q_bar_A0[j]) + (1-cl_test_data$A)*(pred_test_A1 - q_bar_A1[j] - cl_test_data$Y),
                      ci.r = cl_test_data$A*(cl_test_data$Y - pred_test_A0 + q_bar_A0[j]) + (1-cl_test_data$A)*(pred_test_A1 + q_bar_A1[j] - cl_test_data$Y))
    data.frame(avg_length1 = mean(CI1$ci.r-CI1$ci.l), exclude01 = mean(CI1$ci.r <0))
  })
  
  # marginal individual-level effects -------
  fit_Y <- SuperLearner(observed_data$Y[observed_data$cluster_id %in% index_tr], 
                        X = dplyr::select(observed_data[observed_data$cluster_id %in% index_tr,], -c("Y", "cluster_id")), 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A1,], -c("Y", "cluster_id")))$pred)
  S_A0 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A0,], -c("Y", "cluster_id")))$pred)
  
  q_A1d <- data.frame(S = c(S_A1,1e2),
                      cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A1], m+1),
                      N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A1], 1)) %>%
    mutate(weights = 1/(length(index_ca_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights))
  q_A1 <- map_dbl(alpha, ~filter(q_A1d, cumsum_weights >= 1-.)[1,1])
  
  q_A0d <- data.frame(S = c(S_A0,1e2),
                      cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A0], m+1000),
                      N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A0], 1)) %>%
    mutate(weights = 1/(length(index_ca_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights))
  q_A0 <- map_dbl(alpha, ~filter(q_A0d, cumsum_weights >= 1-.)[1,1])
  
  pred_test_A1 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=0))$pred
  
  in_marginal <- map_dfr(1:length(alpha), function(j){
    CI3 <- data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0[j]) + (1-test_data$A)*(pred_test_A1 - q_A1[j] - test_data$Y),
                      ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0[j]) + (1-test_data$A)*(pred_test_A1 + q_A1[j] - test_data$Y))
    data.frame(avg_length3 = mean(CI3$ci.r-CI3$ci.l), exclude03 = mean(CI3$ci.r <0))
  })
  cbind(cl_marginal, in_marginal)
})
mean_marginal_results <- round(Reduce('+',marginal_results)/sim_size,3)
sd_marginal_results <- sqrt((Reduce('+', map(marginal_results, ~.^2)) - Reduce('+',marginal_results)^2/sim_size)/(sim_size-1)) %>% round(digits = 3)

xtable::xtable(data.frame(length1 = paste0(mean_marginal_results[,1], "(", sd_marginal_results[,1], ")"),
           exclude1 = paste0(mean_marginal_results[,2], "(", sd_marginal_results[,2], ")"),
           length2 = paste0(mean_marginal_results[,3], "(", sd_marginal_results[,3], ")"),
           exclude2 = paste0(mean_marginal_results[,4], "(", sd_marginal_results[,4], ")")))

##### conditinal analysis 1 ---------
tictoc::tic()
# observed_data_local <- filter(observed_data_all, Y0 >=7)
# observed_data_local <- filter(observed_data_all, Y0 >= 4 & Y0 < 7)
observed_data_local <- filter(observed_data_all, Y0 < 4)
# observed_data_local <- filter(observed_data_all, Dep_OR_Anx==1)
# observed_data_local <- filter(observed_data_all, FEMALE==1) 
# observed_data_local <- filter(observed_data_all, AGE>=60)
local_N <- aggregate(observed_data_local$N, by=list(cluster_id = observed_data_local$cluster_id), length)
observed_data_local <- left_join(observed_data_local, local_N, by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
cl_data_local <- observed_data_local %>% group_by(cluster_id) %>% summarise_all(mean) %>% as.data.frame
local_indi_results <- map(1:sim_size, function(j){
  index_test <- sample(cl_data_local$cluster_id, size = 14, replace = F)
  test_data <- filter(observed_data_local, cluster_id %in% index_test)
  cl_test_data <- filter(cl_data_local, cluster_id %in% index_test)
  observed_data <- filter(observed_data_local, !(cluster_id %in% index_test))
  cl_data <- filter(cl_data_local, !(cluster_id %in% index_test))
  m <- length(unique(cl_data$cluster_id))
  index_A1 <- cl_data$cluster_id[cl_data$A == 1]
  index_A0 <- cl_data$cluster_id[cl_data$A == 0]
  index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
  index_ca <- setdiff(cl_data$cluster_id, index_tr)
  index_ca_A1 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 1]
  index_ca_A0 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 0]
  
  fit_Y <- SuperLearner(observed_data$Y[observed_data$cluster_id %in% index_tr], 
                        X = dplyr::select(observed_data[observed_data$cluster_id %in% index_tr,], -c("Y", "cluster_id")), 
                        family = "gaussian", 
                        SL.library = SLmethods)
  S_A1 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A1] - 
                predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A1,], -c("Y", "cluster_id")))$pred)
  S_A0 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A0] - 
                predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A0,], -c("Y", "cluster_id")))$pred)
  
  q_A1d <- data.frame(S = c(S_A1,1e2),
                      cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A1], m+1000),
                      N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A1], 1)) %>%
    mutate(weights = 1/(length(index_ca_A1)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights))
  q_A1 <- map_dbl(alpha, ~filter(q_A1d, cumsum_weights >= 1-.)[1,1])
  
  q_A0d <- data.frame(S = c(S_A0,1e2),
                      cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A0], m+1),
                      N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A0], 1)) %>%
    mutate(weights = 1/(length(index_ca_A0)+1) * 1/N) %>%
    arrange(S) %>%
    mutate(cumsum_weights = cumsum(weights))
  q_A0 <- map_dbl(alpha, ~filter(q_A0d, cumsum_weights >= 1-.)[1,1])
  
  pred_test_A1 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=1))$pred
  pred_test_A0 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=0))$pred
  
  in_local <- map_dfr(1:length(alpha), function(j){
    CI3 <- data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0[j]) + (1-test_data$A)*(pred_test_A1 - q_A1[j] - test_data$Y),
                      ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0[j]) + (1-test_data$A)*(pred_test_A1 + q_A1[j] - test_data$Y))
    data.frame(avg_length3 = mean(CI3$ci.r-CI3$ci.l), exclude03 = mean(CI3$ci.r <0))
  })
})
mean_local_results <- (Reduce('+',local_indi_results)/sim_size) %>% round(3)
sd_local_results <- sqrt((Reduce('+', map(local_indi_results, ~.^2)) - Reduce('+',local_indi_results)^2/sim_size)/(sim_size-1)) %>% round(digits = 3)
xtable::xtable(data.frame(length1 = paste0(mean_local_results[,1], "(", sd_local_results[,1], ")"),
                          exclude1 = paste0(mean_local_results[,2], "(", sd_local_results[,2], ")")))
tictoc::toc()




##### conditional analysis 2 ---------
tictoc::tic()
alpha <- 0.2
subgroups <- expression(Y0<=4, Y0<=5 & Y0>4, Y0<=6 & Y0>5, Y0<7 & Y0>6, Y0<8 & Y0>=7, Y0<9 & Y0>=8, Y0>=9)
condi_result_2 <- map_dfr(1: length(subgroups), function(iter){
  observed_data_local <- filter(observed_data_all, eval(subgroups[[iter]]))
  local_N <- aggregate(observed_data_local$N, by=list(cluster_id = observed_data_local$cluster_id), length)
  observed_data_local <- left_join(observed_data_local, local_N, by = join_by(cluster_id)) %>% mutate(N=x) %>% select(-x)
  cl_data_local <- observed_data_local %>% group_by(cluster_id) %>% summarise_all(mean) %>% as.data.frame
  local_indi_results <- map_dfr(1:sim_size, function(j){
    teset_size <- floor(nrow(cl_data_local)/10)
    index_test <- sample(cl_data_local$cluster_id, size = teset_size, replace = F)
    test_data <- filter(observed_data_local, cluster_id %in% index_test)
    cl_test_data <- filter(cl_data_local, cluster_id %in% index_test)
    observed_data <- filter(observed_data_local, !(cluster_id %in% index_test))
    cl_data <- filter(cl_data_local, !(cluster_id %in% index_test))
    m <- length(unique(cl_data$cluster_id))
    index_A1 <- cl_data$cluster_id[cl_data$A == 1]
    index_A0 <- cl_data$cluster_id[cl_data$A == 0]
    index_tr <- c(sample(index_A1, size = floor(tr_prop * length(index_A1))), sample(index_A0, size = floor(tr_prop * length(index_A0))))
    index_ca <- setdiff(cl_data$cluster_id, index_tr)
    index_ca_A1 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 1]
    index_ca_A0 <- index_ca[cl_data$A[cl_data$cluster_id %in% index_ca] == 0]
    
    fit_Y <- SuperLearner(observed_data$Y[observed_data$cluster_id %in% index_tr], 
                          X = dplyr::select(observed_data[observed_data$cluster_id %in% index_tr,], -c("Y", "cluster_id")), 
                          family = "gaussian", 
                          SL.library = SLmethods)
    S_A1 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A1] - 
                  predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A1,], -c("Y", "cluster_id")))$pred)
    S_A0 <- abs(observed_data$Y[observed_data$cluster_id %in% index_ca_A0] - 
                  predict(fit_Y, newdata = dplyr::select(observed_data[observed_data$cluster_id %in% index_ca_A0,], -c("Y", "cluster_id")))$pred)
    
    q_A1d <- data.frame(S = c(S_A1,1e2),
                        cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A1], m+1000),
                        N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A1], 1)) %>%
      mutate(weights = 1/(length(index_ca_A1)+1) * 1/N) %>%
      arrange(S) %>%
      mutate(cumsum_weights = cumsum(weights))
    q_A1 <- filter(q_A1d, cumsum_weights >= 1-alpha)[1,1]
    
    q_A0d <- data.frame(S = c(S_A0,1e2),
                        cluster_id = c(observed_data$cluster_id[observed_data$cluster_id %in% index_ca_A0], m+1),
                        N = c(observed_data$N[observed_data$cluster_id %in% index_ca_A0], 1)) %>%
      mutate(weights = 1/(length(index_ca_A0)+1) * 1/N) %>%
      arrange(S) %>%
      mutate(cumsum_weights = cumsum(weights))
    q_A0 <-filter(q_A0d, cumsum_weights >= 1-alpha)[1,1]
    
    pred_test_A1 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=1))$pred
    pred_test_A0 <- predict(fit_Y, newdata = dplyr::select(test_data, -c("Y", "cluster_id")) %>% mutate(A=0))$pred
    
    data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0) + (1-test_data$A)*(pred_test_A1 - q_A1 - test_data$Y),
               ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0) + (1-test_data$A)*(pred_test_A1 + q_A1 - test_data$Y))
    
    # in_local <- map_dfr(1:length(alpha), function(j){
    #   CI3 <- data.frame(ci.l = test_data$A*(test_data$Y - pred_test_A0 - q_A0[j]) + (1-test_data$A)*(pred_test_A1 - q_A1[j] - test_data$Y),
    #                     ci.r = test_data$A*(test_data$Y - pred_test_A0 + q_A0[j]) + (1-test_data$A)*(pred_test_A1 + q_A1[j] - test_data$Y))
    #   data.frame(avg_length3 = mean(CI3$ci.r-CI3$ci.l), exclude03 = mean(CI3$ci.r <0))
  })
  mutate(local_indi_results) %>% mutate(subgroup = as.character(subgroups[iter]))
})
tictoc::toc()
plot_d2 <- condi_result_2 %>% group_by(subgroup) %>% summarise_all(mean)
plot_d2$subgroup <- c("(6,7)", "[7,8)", "[8,9)", "(1,4]", "(4,5]", "(5,6]", "[9,10)")
ggplot(plot_d2) +
  geom_errorbar(aes(x=subgroup, ymin = ci.l, ymax = ci.r), width = 0.3) +
  ylim(-5,5) + labs(x= "Baseline pain score", y = "80% conformal interval for individual-level treatment effect\n averaged across 100 runs") +
  theme_bw() +
  theme(text= element_text(size = 16)) +
  coord_flip()
ggsave("data-application/caterpillar.png")  

