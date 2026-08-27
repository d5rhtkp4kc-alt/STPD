## Description: simulation results on prediction
## get required statistical model results needed for Figures 4, S3, and S4.

library(movMF)
library(parallel)
library(foreach)
library(doSNOW)


source('../main function/basic_function.R')


#############################################################################################
#############################################################################################
## Parametric model prediction results needed for Figures 4 and S3
#############################################################################################
#############################################################################################

sample_size_set <- c(120, 300, 600)
true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4))

sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  for (i in 1:length(sample_size_set)) {
    for (true_AR_coef_current in true_AR_coef_set) {
      results <- full_workflow_prediction_function(
        sample.size = sample_size_set[i],
        rep = rep,
        true_AR_coef = true_AR_coef_current,
        d_set=7
      )
      results <- results[, ncol(results):1]
      result_list[[list_ind]] <- results
      list_ind <- list_ind + 1
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

save(sim_results,file = 'prediction_sim_result_save_d7.RData')

#############################################################################################
#############################################################################################
## Parametric model prediction results needed for Figures S4
#############################################################################################
#############################################################################################

sample_size_set <- c(120, 300, 600)
true_AR_coef_set <- list(c(0.5),c(0.6,-0.3),c(0.3,-0.1,0.4))


sim_function <- function(rep) {
  result_list <- list()
  list_ind <- 1
  for (i in 1:length(sample_size_set)) {
    for (true_AR_coef_current in true_AR_coef_set) {
      results <- full_workflow_prediction_function(
        sample.size = sample_size_set[i],
        rep = rep,
        true_AR_coef = true_AR_coef_current,
        d_set=48
      )
      results <- results[, ncol(results):1]
      result_list[[list_ind]] <- results
      list_ind <- list_ind + 1
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

save(sim_results,file = 'prediction_sim_result_save_d48.RData')

