## Description: simulation results on prediction
## generate required simulated data for neural network results needed for Figures 4, S3, and S4.

library(movMF)
library(parallel)
library(foreach)
library(doSNOW)

source('../main function/basic_function.R')

#############################################################################################
#############################################################################################
## Data for neural networks needed for Figures 4 and S3
#############################################################################################
#############################################################################################
sample_size_set <- c(120, 300, 600)
true_AR_coef_set <- list(c(0.5), c(0.6, -0.3), c(0.3, -0.1, 0.4))
d_set <- 7
num_replicates <- 200

# Base directory where all folders will be stored
base_dir <- "simulated_datasets_7"

if (!dir.exists(base_dir)) {
  dir.create(base_dir)
}

for (sample.size in sample_size_set) {
  for (i in seq_along(true_AR_coef_set)) {
    
    current_coef <- true_AR_coef_set[[i]]
    
    # Create folder name using AR1, AR2, AR3 based on the index 'i'
    folder_name <- sprintf("n_%d_AR%d", sample.size, i)
    folder_path <- file.path(base_dir, folder_name)
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    
    cat(sprintf("Generating %d datasets for %s...\n", num_replicates, folder_name))
    

    for (rep in 1:num_replicates) {
      AR_list <- AR_generate_function(
        seed_set = rep, 
        sample_size = sample.size, 
        d = d_set, 
        coef = current_coef
      )
      
      sphere_matrix <- matrix(NA, nrow = sample.size, ncol = d_set)
      for (t in 1:sample.size) {
        mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
        angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
        sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]], angle_n) %*% mu
      }
      
      period_matrix <- period_generate_funtion(AR_list)
      y_matrix <- trend_generate_funtion(period_matrix)
      
      file_name <- sprintf("sim_%03d.csv", rep)
      file_path_full <- file.path(folder_path, file_name)
      
      # Write to CSV
      write.table(y_matrix, file = file_path_full, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
    }
  }
}

cat("All datasets generated successfully!\n")


#############################################################################################
#############################################################################################
## Data for neural networks needed for  Figures S4
#############################################################################################
#############################################################################################
sample_size_set <- c(120, 300, 600)
true_AR_coef_set <- list(c(0.5), c(0.6, -0.3), c(0.3, -0.1, 0.4))
d_set <- 48
num_replicates <- 200

# Base directory where all folders will be stored
base_dir <- "simulated_datasets_48"

if (!dir.exists(base_dir)) {
  dir.create(base_dir)
}

for (sample.size in sample_size_set) {
  for (i in seq_along(true_AR_coef_set)) {
    
    current_coef <- true_AR_coef_set[[i]]
    
    # Create folder name using AR1, AR2, AR3 based on the index 'i'
    folder_name <- sprintf("n_%d_AR%d", sample.size, i)
    folder_path <- file.path(base_dir, folder_name)
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    
    cat(sprintf("Generating %d datasets for %s...\n", num_replicates, folder_name))
    
    
    for (rep in 1:num_replicates) {
      AR_list <- AR_generate_function(
        seed_set = rep, 
        sample_size = sample.size, 
        d = d_set, 
        coef = current_coef
      )
      
      sphere_matrix <- matrix(NA, nrow = sample.size, ncol = d_set)
      for (t in 1:sample.size) {
        mu <- c(1:d_set) / sqrt(sum((c(1:d_set))^2))
        angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
        sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]], angle_n) %*% mu
      }
      
      period_matrix <- period_generate_funtion(AR_list)
      y_matrix <- trend_generate_funtion(period_matrix)
      
      file_name <- sprintf("sim_%03d.csv", rep)
      file_path_full <- file.path(folder_path, file_name)
      
      # Write to CSV
      write.table(y_matrix, file = file_path_full, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
    }
  }
}

cat("All datasets generated successfully!\n")