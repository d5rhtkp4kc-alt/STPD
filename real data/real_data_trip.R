## Description: real application for New York City Citi bike sharing system
## reproduce Figures 7, 8, S6, S7, S8, S9, and S10, 
## and to generate required statistical model results needed for Figure S11.


library(ggplot2)
library(patchwork)
library(psych)
library(osqp)
library(transport)
library(tidyr)
library(plot3D)

source('../main function/basic_function.R')

#############################################################################################
#############################################################################################
## Model fitting
#############################################################################################
#############################################################################################

y_matrix <- read.csv("dataset/trip.csv", header = FALSE)
y_matrix <- as.matrix(y_matrix)
sample_size = nrow(y_matrix)

#### Trend Estimation ####
x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1)
trend_component <- LocSpheReg(xin = as.vector(x_input), yin = y_matrix)
bw_current <- trend_component$optns$bw ##0.3317327
trend_component <- trend_component$yout
mu <- intrinsic_mean_sphere(y_matrix)
remove_trend_result <- get_spherical_residuals(y_matrix, trend_component,mu = mu)
remove_trend_data <- remove_trend_result$residuals

#### Period Estimation ####

loss_vec_save_np <- period_est_sphere_function(remove_trend_data,25,sample_size,0)
period_max=25
lambda_max = 1
lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,10)
BIC_ind = which(BIC_vec==min(BIC_vec))[14]
lambda_choose = lambda_seq[BIC_ind]
loss_vec_save_wp = loss_vec_save_np + seq(1,period_max) * lambda_choose

#### Periodic Component Estimation ####
p=7
x <- covariate_function(p,sample_size)
periodic_component <- GloSpheGeoReg(xin = x, yin = remove_trend_data)
final_result <- get_spherical_residuals(remove_trend_data, periodic_component)
final_resid_data <- final_result$residuals


## check periodicity
loss_vec_save_np <- period_est_sphere_function(final_resid_data,25,sample_size,0)
period_max=25
lambda_max = 1
lambda_seq = seq(0.0001,lambda_max,lambda_max/200)
BIC_vec = BIC_function_new(loss_vec_save_np,lambda_seq,sample_size,bw_current,10)
BIC_ind = which(BIC_vec==min(BIC_vec))[14]
lambda_choose = lambda_seq[BIC_ind]
loss_vec_save_wp2 = loss_vec_save_np + seq(1,period_max) * lambda_choose


#############################################################################################
#############################################################################################
## Figure 7
#############################################################################################
#############################################################################################
loss_df <- data.frame(loss_vec_save_wp)
period_max=25
colnames(loss_df) <- 'loss_value'
p1 <- ggplot(loss_df, aes(x = seq(1,period_max), y = loss_value)) + 
  geom_line() +
  labs(
    title = '(a) Penalized RSS for de-trended time series',
    x = expression(vartheta),  # Add LaTeX-style expression
    y = 'Penalized RSS'
  ) +
  scale_x_continuous(limits = c(1, period_max),breaks = c(1, 7, 14, 21)) +  # Custom x-axis tick marks
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 14),
        plot.margin=margin(5,15,5,5),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))

loss_df <- data.frame(loss_vec_save_wp2)
period_max=25
colnames(loss_df) <- 'loss_value'
p2 <- ggplot(loss_df, aes(x = seq(1,period_max), y = loss_value)) + 
  geom_line() +
  labs(
    title = '(b) Penalized RSS for de-trended and de-seasonalized time series',
    x = expression(vartheta),  # Add LaTeX-style expression
    y = 'Penalized RSS'
  ) +
  scale_x_continuous(limits = c(1, period_max),breaks =c(1, 7, 14, 21)) +  # Custom x-axis tick marks
  theme_classic()+
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text=element_text(size = 12),
        legend.text=element_text(size = 14),
        plot.margin=margin(5,15,5,5),
        panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
        legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))

p1 + p2

#############################################################################################
#############################################################################################
## Figure S6
#############################################################################################
#############################################################################################

plot_function <- function(matrix_input, select.ind, title.input='', y.title='',line.adj=-2) {
  # 1. Subset Data
  matrix_input <- matrix_input[select.ind, ]
  
  # 2. Setup Grids
  grp_grid <- seq(0, 24, length.out = ncol(matrix_input))
  hours_grid <- 1:nrow(matrix_input)
  z_final <- t(matrix_input^2)
  
  # 3. Colors
  my_cols <- colorRampPalette(c("#E0F7FA", "#80DEEA", "#FFF59D", "#FFEB3B"))(100)
  
  # 4. Plotting
  par(mar = c(1.5, 1.5, 1, 0.5))
  
  res <- persp3D(x = grp_grid, 
                 y = hours_grid, 
                 z = z_final,
                 col = my_cols, 
                 colkey = FALSE, 
                 facets = TRUE,
                 border = "#4682B4", 
                 lwd = 0.1, 
                 alpha = 0.6,
                 theta = 150, phi = 25, 
                 expand = 0.5, 
                 bty = "b2",
                 ticktype = "detailed",
                 zlim = c(0, 0.08),
                 clim = c(0, 0.08),
                 # REMOVED main = title.input from here
                 xlab = "Time (Hour of Day)", 
                 ylab = "Observation", 
                 zlab = "")
  
  mtext(title.input, side = 3, line = line.adj, font = 1, cex = 1)
  title(ylab = y.title, font.lab = 1, line = -1, cex.lab = 1.5)
}

par(mfrow = c(1, 3)) 

plot_function(y_matrix,c(1:nrow(y_matrix)),'(a) Original time series','')
plot_function(remove_trend_data,c(1:nrow(remove_trend_data)),'(b) De-trended time series','')
plot_function(final_resid_data,c(1:nrow(final_resid_data)),'(c) De-trended and de-seasonalized time series','')

#############################################################################################
#############################################################################################
## Figure 8
#############################################################################################
#############################################################################################
dates <- seq(as.Date("2019-03-01"), length.out = nrow(y_matrix), by = "day")
day_names <- weekdays(dates)
weekday_ind <- !(day_names %in% c("Saturday", "Sunday"))
weekend_ind <- day_names %in% c("Saturday", "Sunday")

par(mfrow = c(2, 3)) 

plot_function(y_matrix,weekday_ind,'Original time series','Weekday',line.adj=-0.5)
plot_function(remove_trend_data,weekday_ind,'De-trended time series','',line.adj=-0.5)
plot_function(final_resid_data,weekday_ind,'De-trended and de-seasonalized time series','',line.adj=-0.5)
plot_function(y_matrix,weekend_ind,'','Weekend')
plot_function(remove_trend_data,weekend_ind,'','')
plot_function(final_resid_data,weekend_ind,'','')


#############################################################################################
#############################################################################################
## Figure S7
#############################################################################################
#############################################################################################
dates <- seq(as.Date("2019-03-01"), length.out = nrow(y_matrix), by = "day")
day_names <- weekdays(dates)


par(mfrow = c(1, 3)) 

plot_function(y_matrix,day_names %in% c("Monday"),'Original time series','Monday',line.adj=-0.5)
plot_function(remove_trend_data,day_names %in% c("Monday"),'De-trended time series','',line.adj=-0.5)
plot_function(final_resid_data,day_names %in% c("Monday"),'De-trended and de-seasonalized time series','',line.adj=-0.5)

plot_function(y_matrix,day_names %in% c("Tuesday"),'','Tuesday')
plot_function(remove_trend_data,day_names %in% c("Tuesday"),'','')
plot_function(final_resid_data,day_names %in% c("Tuesday"),'','')

plot_function(y_matrix,day_names %in% c("Wednesday"),'','Wednesday')
plot_function(remove_trend_data,day_names %in% c("Wednesday"),'','')
plot_function(final_resid_data,day_names %in% c("Wednesday"),'','')

plot_function(y_matrix,day_names %in% c("Thursday"),'','Thursday')
plot_function(remove_trend_data,day_names %in% c("Thursday"),'','')
plot_function(final_resid_data,day_names %in% c("Thursday"),'','')

#############################################################################################
#############################################################################################
## Figure S8
#############################################################################################
#############################################################################################

plot_function(y_matrix,day_names %in% c("Friday"),'Original time series','Friday',line.adj=-0.5)
plot_function(remove_trend_data,day_names %in% c("Friday"),'De-trended time series','',line.adj=-0.5)
plot_function(final_resid_data,day_names %in% c("Friday"),'De-trended and de-seasonalized time series','',line.adj=-0.5)

plot_function(y_matrix,day_names %in% c("Saturday"),'','Saturday')
plot_function(remove_trend_data,day_names %in% c("Saturday"),'','')
plot_function(final_resid_data,day_names %in% c("Saturday"),'','')


plot_function(y_matrix,day_names %in% c("Sunday"),'','Sunday')
plot_function(remove_trend_data,day_names %in% c("Sunday"),'','')
plot_function(final_resid_data,day_names %in% c("Sunday"),'','')


#############################################################################################
#############################################################################################
## Figure S9
#############################################################################################
#############################################################################################

par(mfrow = c(1, 2)) 
plot_function(trend_component,c(1:nrow(trend_component)),'(a) Estimated trend component')
plot_function(periodic_component,c(1:nrow(periodic_component)),'(b) Estimated periodic component')


par(mfrow = c(2, 2)) 

plot_function(trend_component,weekday_ind,'Estimated trend component','Weekday',line.adj=-0.5)
plot_function(periodic_component,weekday_ind,'Estimated periodic component','',line.adj=-0.5)
plot_function(trend_component,weekend_ind,'','Weekend')
plot_function(periodic_component,weekend_ind,'','')



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
    eval.method = 'moving',
    error.measure='JCD'
  )
  zhu_DAR_predict <- difference_AR_predict_function(
    data = y_matrix_current,
    candidate.order = order_upper_set,
    test.num = (test_num+1-i),
    predict.step.num = (test_num+1-i),
    eval.method = 'moving',
    error.measure='JCD'
  )
  PTAR_predict <- PTAR_predict_function(
    data = y_matrix_current,
    period_est = 7,
    candidate.order = order_upper_set,
    test.num = (test_num+1-i),
    predict.step.num = (test_num+1-i),
    eval.method = 'moving',
    error.measure='JCD'
  )
  prediction_error_matrix[1,i] <- mean(zhu_AR_predict[[1]]) 
  prediction_error_matrix[2,i] <- mean(zhu_DAR_predict[[1]]) 
  prediction_error_matrix[3,i] <- mean(PTAR_predict[[1]]) 
  print(c(mean(zhu_AR_predict[[1]]),mean(zhu_DAR_predict[[1]]),mean(PTAR_predict[[1]])))
}

prediction_error_matrix <- prediction_error_matrix[, ncol(prediction_error_matrix):1]
save(prediction_error_matrix,file = 'transport_pred_30.RData')

