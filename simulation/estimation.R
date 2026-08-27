## Description: simulation results on estimation
## reproduce Tables S1, S2, and S3. 

library(movMF)
library(trust)
library(parallel)
library(foreach)
library(doSNOW)
library(xtable)
library(ggplot2)
library(patchwork)

source('../main function/basic_function.R')


#############################################################################################
#############################################################################################
## Table S1
#############################################################################################
#############################################################################################

sim_setting = 'full'

if (sim_setting == 'full'){
  sample_size_set <- c(120, 300, 600, 1200, 2400)
  period_upper_bound_set <- c(40, 60, 80, 100, 120)
  true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4)) 
} else{
  sample_size_set <- c(120)
  period_upper_bound_set <- c(40)
  true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4)) 
}


sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  for (i in 1:length(sample_size_set)) {
    for (true_AR_coef_current in true_AR_coef_set) {
      results <- full_workflow_MSE_function(
        sample.size = sample_size_set[i],
        rep = rep,
        true_AR_coef = true_AR_coef_current,
        period_upper = period_upper_bound_set[i],
        d_set = 7
      )
      result_list[[list_ind]] <- results
      list_ind <- list_ind+1
    }
  }
  return(result_list)
}

# Setup Cluster
nCores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(nCores)
registerDoSNOW(cl)

iterations <- 200

# Progress Bar
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# Parallel Loop
sim_results <- foreach(
  i = 1:iterations,
  .combine = 'c',               
  .multicombine = TRUE,           
  .options.snow = opts,
  .packages = c("movMF")          
) %dopar% {
  source('../main function/basic_function.R')
  list(sim_function(rep = i))
}

parallel::stopCluster(cl)
close(pb)

save(sim_results,file = 'est_sim_result_save.RData')

results_list <- get(load("est_sim_result_save.RData"))



result_matrix <- matrix(NA,nrow = 21,ncol = 5)
row.names(result_matrix) <- c('trend component','periodic component','residual',
                              'AR(1) parameter bias','AR(1) parameter MSE',
                              'trend component','periodic component','residual',
                              'AR(2) 1 parameter bias','AR(2) 1 parameter MSE',
                              'AR(2) 2 parameter bias','AR(2) 2 parameter MSE',
                              'trend component','periodic component','residual',
                              'AR(3) 1 parameter bias','AR(3) 1 parameter MSE',
                              'AR(3) 2 parameter bias','AR(3) 2 parameter MSE',
                              'AR(3) 3 parameter bias','AR(3) 3 parameter MSE')

AR2_coef <- c(0.6,-0.3)
AR3_coef <- c(0.3,-0.1,0.4)

for (i in 1:5) {
  ## AR(1)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-2)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-2)]][[2]])
  }  
  result_matrix[1:3,i] <- colMeans(SE_current)
  result_matrix[4,i] <- mean(AR_coef-0.5)
  result_matrix[5,i] <- mean((AR_coef - 0.5)^2)
  ## AR(2)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-1)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-1)]][[2]])
  }  
  result_matrix[6:8,i] <- colMeans(SE_current)
  result_matrix[9,i] <- mean(AR_coef[1,]-AR2_coef[1])
  result_matrix[10,i] <- mean((AR_coef[1,] - AR2_coef[1])^2)
  result_matrix[11,i] <- mean(AR_coef[2,]-AR2_coef[2])
  result_matrix[12,i] <- mean((AR_coef[2,] - AR2_coef[2])^2)
  ## AR(3)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3)]][[2]])
  }  
  result_matrix[13:15,i] <- colMeans(SE_current)
  result_matrix[16,i] <- mean(AR_coef[1,]-AR3_coef[1])
  result_matrix[17,i] <- mean((AR_coef[1,] - AR3_coef[1])^2)
  result_matrix[18,i] <- mean(AR_coef[2,]-AR3_coef[2])
  result_matrix[19,i] <- mean((AR_coef[2,] - AR3_coef[2])^2)
  result_matrix[20,i] <- mean(AR_coef[3,]-AR3_coef[3])
  result_matrix[21,i] <- mean((AR_coef[3,] - AR3_coef[3])^2)
}


formatted_mat <- apply(result_matrix * 100, 2, function(x) sprintf("%.4f", x))
formatted_df <- as.data.frame(formatted_mat)
print(xtable(formatted_df, align = "lccccc"), 
      include.rownames = TRUE, 
      sanitize.text.function = identity)


## period estimation
list_ind <- 1
period_est_list <- list()

for (i in 1:5) {
  ## AR(1)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3-2)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3-2)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3-2)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  ## AR(2)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3-1)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3-1)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3-1)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  ## AR(3)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  
}

# print results for period estimation

for (i in c(1:5)) {
  ## AR(1)
  current_result <- period_est_list[[(i*3-2)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
  ## AR(2)
  current_result <- period_est_list[[(i*3-1)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
  ## AR(3)
  current_result <- period_est_list[[(i*3)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
}


#############################################################################################
#############################################################################################
## Table S2
#############################################################################################
#############################################################################################

sample_size_set <- c(120, 300, 600, 1200, 2400)
period_upper_bound_set <- c(40, 60, 80, 100, 120)
true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4))


sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  for (i in 1:length(sample_size_set)) {
    for (true_AR_coef_current in true_AR_coef_set) {
      results <- full_workflow_MSE_function(
        sample.size = sample_size_set[i],
        rep = rep,
        true_AR_coef = true_AR_coef_current,
        period_upper = period_upper_bound_set[i],
        d_set = 48
      )
      result_list[[list_ind]] <- results
      list_ind <- list_ind+1
    }
  }
  return(result_list)
}

# Setup Cluster
nCores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(nCores)
registerDoSNOW(cl)

iterations <- 200

# Progress Bar
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# Parallel Loop
sim_results <- foreach(
  i = 1:iterations,
  .combine = 'c',               
  .multicombine = TRUE,           
  .options.snow = opts,
  .packages = c("movMF")          
) %dopar% {
  source('../main function/basic_function.R')
  list(sim_function(rep = i))
}

parallel::stopCluster(cl)
close(pb)

save(sim_results,file = 'est_sim_result_save.RData')

results_list <- get(load("est_sim_result_save.RData"))



result_matrix <- matrix(NA,nrow = 21,ncol = 5)
row.names(result_matrix) <- c('trend component','periodic component','residual',
                              'AR(1) parameter bias','AR(1) parameter MSE',
                              'trend component','periodic component','residual',
                              'AR(2) 1 parameter bias','AR(2) 1 parameter MSE',
                              'AR(2) 2 parameter bias','AR(2) 2 parameter MSE',
                              'trend component','periodic component','residual',
                              'AR(3) 1 parameter bias','AR(3) 1 parameter MSE',
                              'AR(3) 2 parameter bias','AR(3) 2 parameter MSE',
                              'AR(3) 3 parameter bias','AR(3) 3 parameter MSE')

AR2_coef <- c(0.6,-0.3)
AR3_coef <- c(0.3,-0.1,0.4)

for (i in 1:5) {
  ## AR(1)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-2)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-2)]][[2]])
  }  
  result_matrix[1:3,i] <- colMeans(SE_current)
  result_matrix[4,i] <- mean(AR_coef-0.5)
  result_matrix[5,i] <- mean((AR_coef - 0.5)^2)
  ## AR(2)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-1)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-1)]][[2]])
  }  
  result_matrix[6:8,i] <- colMeans(SE_current)
  result_matrix[9,i] <- mean(AR_coef[1,]-AR2_coef[1])
  result_matrix[10,i] <- mean((AR_coef[1,] - AR2_coef[1])^2)
  result_matrix[11,i] <- mean(AR_coef[2,]-AR2_coef[2])
  result_matrix[12,i] <- mean((AR_coef[2,] - AR2_coef[2])^2)
  ## AR(3)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3)]][[2]])
  }  
  result_matrix[13:15,i] <- colMeans(SE_current)
  result_matrix[16,i] <- mean(AR_coef[1,]-AR3_coef[1])
  result_matrix[17,i] <- mean((AR_coef[1,] - AR3_coef[1])^2)
  result_matrix[18,i] <- mean(AR_coef[2,]-AR3_coef[2])
  result_matrix[19,i] <- mean((AR_coef[2,] - AR3_coef[2])^2)
  result_matrix[20,i] <- mean(AR_coef[3,]-AR3_coef[3])
  result_matrix[21,i] <- mean((AR_coef[3,] - AR3_coef[3])^2)
}


formatted_mat <- apply(result_matrix * 100, 2, function(x) sprintf("%.4f", x))
formatted_df <- as.data.frame(formatted_mat)
print(xtable(formatted_df, align = "lccccc"), 
      include.rownames = TRUE, 
      sanitize.text.function = identity)


## period estimation
sample_size_set <- c(120, 300, 600, 1200, 2400)
period_upper_bound_set <- c(40, 60, 80, 100, 120)
list_ind <- 1
period_est_list <- list()

for (i in 1:5) {
  ## AR(1)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3-2)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3-2)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    ##BIC_vec = BIC_function_trandition(loss_vec_save_np,lambda_seq,sample_size,100,period_max)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3-2)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    ##BIC_vec = BIC_function_trandition(loss_vec_save_np,lambda_seq,sample_size,100,period_max)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  ## AR(2)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3-1)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3-1)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    ##BIC_vec = BIC_function_trandition(loss_vec_save_np,lambda_seq,sample_size,100,period_max)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3-1)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  ## AR(3)
  period_est_vec1 <- NULL
  period_est_vec2 <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    
    bw_current <- result_current[[(i*3)]][[5]]
    sample_size <- sample_size_set[i]
    period_max <- period_upper_bound_set[i]
    
    ## period est after trend removal
    loss_vec_save_np <- result_current[[(i*3)]][[3]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current1 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec1 <- c(period_est_vec1,period_est_current1)
    
    
    ## period est after trend and periodic component removal
    loss_vec_save_np <- result_current[[(i*3)]][[4]]
    lambda_max = 1
    lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
    BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,1)
    BIC_ind = which(BIC_vec==min(BIC_vec))[14]
    
    lambda_choose = lambda_seq[BIC_ind]
    
    loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose
    ##plot(loss_vec_save_wp2)
    period_est_current2 <- which(loss_vec_save_wp2==min(loss_vec_save_wp2))
    period_est_vec2 <- c(period_est_vec2,period_est_current2)
    
    print(c(t,period_est_current1,period_est_current2))
  }  
  
  period_est_list[[list_ind]] <- rbind(period_est_vec1,period_est_vec2)
  list_ind <- list_ind + 1
  
  
}

# print results for period estimation

for (i in c(1:5)) {
  ## AR(1)
  current_result <- period_est_list[[(i*3-2)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
  ## AR(2)
  current_result <- period_est_list[[(i*3-1)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
  ## AR(3)
  current_result <- period_est_list[[(i*3)]]
  ## before
  period_est_df <- current_result[1,]
  print(mean(period_est_df==12))
  ## after
  period_est_df <- current_result[2,]
  print(mean(period_est_df==1))
}




#############################################################################################
#############################################################################################
## Table S3
#############################################################################################
#############################################################################################


## supplementary extrinsic approach
sample_size_set <- c(120, 300, 600, 1200, 2400)
true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4))

sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  for (i in 1:length(sample_size_set)) {
    for (true_AR_coef_current in true_AR_coef_set) {
      results <- extrinsic_full_workflow_MSE_function(
        sample.size = sample_size_set[i],
        rep = rep,
        true_AR_coef = true_AR_coef_current
      )
      result_list[[list_ind]] <- results
      list_ind <- list_ind+1
    }
  }
  return(result_list)
}

# Setup Cluster
nCores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(nCores)
registerDoSNOW(cl)

iterations <- 200

# Progress Bar
pb <- txtProgressBar(max = iterations, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# Parallel Loop
sim_results <- foreach(
  i = 1:iterations,
  .combine = 'c',               
  .multicombine = TRUE,           
  .options.snow = opts,
  .packages = c("movMF")          
) %dopar% {
  source('../main function/basic_function.R')
  list(sim_function(rep = i))
}

parallel::stopCluster(cl)
close(pb)

save(sim_results,file = 'extrinsic_est_sim_result_save.RData')

results_list <- get(load("extrinsic_est_sim_result_save.RData"))



result_matrix <- matrix(NA,nrow = 18,ncol = 5)
row.names(result_matrix) <- c('periodic component','residual',
                              'AR(1) parameter bias','AR(1) parameter MSE',
                              'periodic component','residual',
                              'AR(2) 1 parameter bias','AR(2) 1 parameter MSE',
                              'AR(2) 2 parameter bias','AR(2) 2 parameter MSE',
                              'periodic component','residual',
                              'AR(3) 1 parameter bias','AR(3) 1 parameter MSE',
                              'AR(3) 2 parameter bias','AR(3) 2 parameter MSE',
                              'AR(3) 3 parameter bias','AR(3) 3 parameter MSE')

AR2_coef <- c(0.6,-0.3)
AR3_coef <- c(0.3,-0.1,0.4)

for (i in 1:5) {
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-2)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-2)]][[2]])
  }  
  result_matrix[1:2,i] <- colMeans(SE_current)
  result_matrix[3,i] <- mean(AR_coef-0.5)
  result_matrix[4,i] <- mean((AR_coef - 0.5)^2)
  ## AR(2)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3-1)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3-1)]][[2]])
  }  
  result_matrix[5:6,i] <- colMeans(SE_current)
  result_matrix[7,i] <- mean(AR_coef[1,]-AR2_coef[1])
  result_matrix[8,i] <- mean((AR_coef[1,] - AR2_coef[1])^2)
  result_matrix[9,i] <- mean(AR_coef[2,]-AR2_coef[2])
  result_matrix[10,i] <- mean((AR_coef[2,] - AR2_coef[2])^2)
  ## AR(3)
  SE_current <- NULL
  AR_coef <- NULL
  for (t in 1:length(results_list)) {
    result_current <- results_list[[t]]
    SE_current <- rbind(SE_current,result_current[[(i*3)]][[1]])
    AR_coef <- cbind(AR_coef,result_current[[(i*3)]][[2]])
  }  
  result_matrix[11:12,i] <- colMeans(SE_current)
  result_matrix[13,i] <- mean(AR_coef[1,]-AR3_coef[1])
  result_matrix[14,i] <- mean((AR_coef[1,] - AR3_coef[1])^2)
  result_matrix[15,i] <- mean(AR_coef[2,]-AR3_coef[2])
  result_matrix[16,i] <- mean((AR_coef[2,] - AR3_coef[2])^2)
  result_matrix[17,i] <- mean(AR_coef[3,]-AR3_coef[3])
  result_matrix[18,i] <- mean((AR_coef[3,] - AR3_coef[3])^2)
}


formatted_mat <- apply(result_matrix * 100, 2, function(x) sprintf("%.4f", x))
formatted_df <- as.data.frame(formatted_mat)
print(xtable(formatted_df, align = "lccccc"), 
      include.rownames = TRUE, 
      sanitize.text.function = identity)

