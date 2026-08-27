## Description: real application for U.S. energy generation compositions
## reproduce Figures 5, 6, and S2, 
## and to generate required statistical model results needed for Figure S5.

library(ggplot2)
library(patchwork)
library(psych)
library(osqp)
library(transport)
library(easyCODA)
library(tidyr)

source('../main function/basic_function.R')


#############################################################################################
#############################################################################################
## Preprocessing
#############################################################################################
#############################################################################################

df <- read.csv('dataset/US_energy.csv')
df[df == "--"] <- 0
period_list <- list()
for (i in 2:dim(df)[2]) {
  current_vec <- as.numeric(df[,i])
  compositional_vec <- rep(0,7)
  #coal
  compositional_vec[1] <- current_vec[1]/sum(current_vec)
  #petroleum
  compositional_vec[2] <- sum(current_vec[2:3])/sum(current_vec)
  #gas
  compositional_vec[3] <- sum(current_vec[4:5])/sum(current_vec)
  #nuclear
  compositional_vec[4] <- current_vec[6]/sum(current_vec)
  #Conventional hydroelectric
  compositional_vec[5] <- current_vec[7]/sum(current_vec)
  #Renewables (wind, geothermal, biomass (total))
  compositional_vec[6] <- sum(current_vec[8:10])/sum(current_vec)
  #solar
  compositional_vec[7] <- sum(current_vec[11:12])/sum(current_vec)
  
  period_list[[i-1]] <- sqrt(compositional_vec)
}


y_matrix <- matrix(0,nrow = length(period_list), ncol = length(period_list[[1]]))
for (i in 1:length(period_list)) {
  y_matrix[i,] <- period_list[[i]]
}

y_matrix <- y_matrix[1:276,] ## use 2001-2023

write.table(y_matrix, file = 'energy.csv', sep = ",", 
            row.names = FALSE, col.names = FALSE)

sample_size = nrow(y_matrix)



#############################################################################################
#############################################################################################
## Model fitting
#############################################################################################
#############################################################################################

#### Trend Estimation ####
x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)
trend_component <- LocSpheReg(xin = as.vector(x_input), yin = y_matrix)
bw_current <- trend_component$optns$bw ##0.1657616
trend_component <- trend_component$yout


mu <- intrinsic_mean_sphere(y_matrix)
remove_trend_result <- get_spherical_residuals(y_matrix, trend_component,mu = mu)
remove_trend_data <- remove_trend_result$residuals


#### Period Estimation ####

loss_vec_save_np <- period_est_sphere_function(remove_trend_data,40,sample_size,0)

period_max=40
lambda_max = 1
lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,10)
BIC_ind = which(BIC_vec==min(BIC_vec))[14]
lambda_choose = lambda_seq[BIC_ind]
loss_vec_save_wp = loss_vec_save_np + seq(1,period_max) * lambda_choose


#### Periodic Component Estimation ####
p=12
x <- covariate_function(p,sample_size)
periodic_component <- GloSpheGeoReg(xin = x, yin = remove_trend_data)
final_result <- get_spherical_residuals(remove_trend_data, periodic_component)
final_resid_data <- final_result$residuals

## check periodicity
loss_vec_save_np <- period_est_sphere_function(final_resid_data,40,sample_size,0)
period_max=40
lambda_max = 1
lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,10)
BIC_ind = which(BIC_vec==min(BIC_vec))[14]
lambda_choose = lambda_seq[BIC_ind]
loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose


#############################################################################################
#############################################################################################
## Figure 5
#############################################################################################
#############################################################################################
month_dates <- seq(as.Date("2001-01-01"), as.Date("2023-12-01"), by = "month")

prepare_plot_data <- function(mat, x_axis_vals) {
  df <- as.data.frame(100 * mat^2)
  colnames(df) <- c('Coal','Petroleum','Gas','Nuclear','Conventional hydroelectric','Renewables','Solar')
  df$x <- x_axis_vals
  return(tidyr::pivot_longer(df, 
                             cols = -x, 
                             names_to = "Energy_Source", 
                             values_to = "value"))
}

custom_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")

create_styled_plot_4year <- function(data, title) {
  
  start_date <- as.Date("2001-01-01")
  end_limit  <- as.Date("2023-12-01")
  month_breaks  <- seq(start_date, end_limit, by = "3 years")
  
  ggplot(data, aes(x = x, y = value, color = Energy_Source, linetype = Energy_Source)) +
    geom_line(linewidth = 0.7) +
    scale_linetype_manual(values = custom_linetypes) +
    scale_x_date(
      breaks = month_breaks,
      date_labels = "%Y %b",
      limits = c(start_date, end_limit)
    ) +
    labs(y = 'Percentage(%)', title = title) +
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text=element_text(size = 12),
          legend.text=element_text(size = 14),
          plot.margin=margin(5,15,5,5),
          panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))
}

## fig b
loss_df <- data.frame(loss_vec_save_wp)
period_max=40
colnames(loss_df) <- 'loss_value'
p2 <- ggplot(loss_df, aes(x = seq(1,period_max), y = loss_value)) + 
  geom_line() +
  labs(
    title = '(b) Penalized RSS for de-trended time series',
    x = expression(vartheta),  
    y = 'Penalized RSS'
  ) +
  scale_x_continuous(limits = c(1, period_max),breaks = c(1, 12, 24, 36)) +  # Custom x-axis tick marks
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 14),
        plot.margin=margin(5,15,5,5),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))

## fig d
loss_df <- data.frame(loss_vec_save_wp2)
period_max=40
colnames(loss_df) <- 'loss_value'
p4 <- ggplot(loss_df, aes(x = seq(1,period_max), y = loss_value)) + 
  geom_line() +
  labs(
    title = '(d) Penalized RSS for de-trended and de-seasonalized time series',
    x = expression(vartheta),  
    y = 'Penalized RSS'
  ) +
  scale_x_continuous(limits = c(1, period_max),breaks = c(1, 12, 24, 36)) +  # Custom x-axis tick marks
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 14),
        plot.margin=margin(5,15,5,5),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))
# fig a: trend component est
p1 <- create_styled_plot_4year(prepare_plot_data(trend_component, month_dates), "(a) Estimated trend component")

# fig c: periodic component est
p3 <- create_styled_plot_4year(prepare_plot_data(periodic_component, month_dates), "(c) Estimated periodic component")

combined_plot <- p1 + p2 + p3 + p4 +
  plot_layout(nrow = 2, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot


#############################################################################################
#############################################################################################
## Figure 6
#############################################################################################
#############################################################################################
create_styled_plot <- function(data, title) {
  
  start_date <- as.Date("2001-01-01")
  end_limit  <- as.Date("2023-12-01")
  month_breaks  <- seq(start_date, end_limit, by = "3 years")
  
  ggplot(data, aes(x = x, y = value, color = Energy_Source, linetype = Energy_Source)) +
    geom_line(linewidth = 0.7) +
    scale_linetype_manual(values = custom_linetypes) +
    scale_x_date(
      breaks = month_breaks,
      date_labels = "%Y %b",
      limits = c(start_date, end_limit)
    ) +
    labs(y = 'Percentage(%)', title = title) +
    theme_classic()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size = 14),
          axis.text=element_text(size = 12),
          legend.text=element_text(size = 14),
          plot.margin=margin(5,15,5,5),
          panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))
}
# Plot (1): original data
p1 <- create_styled_plot(prepare_plot_data(y_matrix, month_dates), "(a) Original time series")

# Plot (2): remove trend
p2 <- create_styled_plot(prepare_plot_data(remove_trend_data, month_dates), "(b) De-trended time series")

# Plot (3): remove period
p3 <- create_styled_plot(prepare_plot_data(final_resid_data, month_dates), "(c) De-trended and de-seasonalized time series")

combined_plot <- p1 + p2 + p3 +
  plot_layout(nrow = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot

#############################################################################################
#############################################################################################
## Figure S2
#############################################################################################
#############################################################################################


clr_trans_function <- function(dat){
  dat_square <- dat^2
  ## deal with zero
  eps <- 1e-6
  y_pos <- dat_square
  y_pos[y_pos == 0] <- eps
  y_pos <- y_pos / rowSums(y_pos)
  y_clr_mat <- CLR(y_pos)$LR
  return(y_clr_mat)
}


inv_clr <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    return(t(apply(x, 1, function(row) {
      exp_row <- exp(row)
      exp_row / sum(exp_row)
    })))
  } else if (is.numeric(x)) {
    exp_x <- exp(x)
    return(exp_x / sum(exp_x))
  } else {
    stop("Input must be a numeric vector, matrix, or data frame.")
  }
}

y_clr_mat <- clr_trans_function(y_matrix)
matplot(y_clr_mat, type = "l", lty = 1, col = 1:nrow(y_clr_mat))

trend_component_clr <- clr_trans_function(trend_component)
period_component_clr <- clr_trans_function(periodic_component)
remove_trend_result <- y_clr_mat-trend_component_clr
remove_trend_data <- sqrt(inv_clr(remove_trend_result))
final_result <- remove_trend_result-period_component_clr
matplot(final_result, type = "l", lty = 1, col = 1:nrow(final_result))
final_resid_data <- sqrt(inv_clr(final_result))

# Plot (1): original data
p1 <- create_styled_plot(prepare_plot_data(y_matrix, month_dates), "(a) Original time series")

# Plot (2): remove trend
p2 <- create_styled_plot(prepare_plot_data(remove_trend_data, month_dates), "(b) De-trended time series")

# Plot (3): remove period
p3 <- create_styled_plot(prepare_plot_data(final_resid_data, month_dates), "(c) De-trended and de-seasonalized time series")

combined_plot <- p1 + p2 + p3 +
  plot_layout(nrow = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot

##########################################################################################
##########################################################################################
## Statistical model prediction
##########################################################################################
##########################################################################################
## rolling from B-steps ahead prediction to 1-step ahead prediction

order_upper_set <- c(1:20)
test_num <- 30
prediction_error_matrix <- matrix(NA,nrow = 3,ncol = test_num)

for (i in 1:test_num) {
  y_matrix_current <- y_matrix[i:nrow(y_matrix),]
  zhu_AR_predict <- stationary_AR_predict_function(
    data = y_matrix_current,
    candidate.order = order_upper_set,
    test.num = (test_num+1-i),
    predict.step.num = (test_num+1-i),
    eval.method = 'moving'
  )
  zhu_DAR_predict <- difference_AR_predict_function(
    data = y_matrix_current,
    candidate.order = order_upper_set,
    test.num = (test_num+1-i),
    predict.step.num = (test_num+1-i),
    eval.method = 'moving'
  )
  PTAR_predict <- PTAR_predict_function(
    data = y_matrix_current,
    period_est = 12,
    candidate.order = order_upper_set,
    test.num = (test_num+1-i),
    predict.step.num = (test_num+1-i),
    eval.method = 'moving'
  )
  prediction_error_matrix[1,i] <- mean(zhu_AR_predict[[1]]) 
  prediction_error_matrix[2,i] <- mean(zhu_DAR_predict[[1]]) 
  prediction_error_matrix[3,i] <- mean(PTAR_predict[[1]]) 
  print(c(mean(zhu_AR_predict[[1]]),mean(zhu_DAR_predict[[1]]),mean(PTAR_predict[[1]])))
}

prediction_error_matrix <- prediction_error_matrix[, ncol(prediction_error_matrix):1]
save(prediction_error_matrix,file = 'US_energy_pred_30.RData')


