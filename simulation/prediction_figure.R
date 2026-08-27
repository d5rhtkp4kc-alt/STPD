## Description: simulation results on prediction
## reproduce Figures 4, S3, and S4.

library(ggplot2)
library(patchwork)
library(scales)

#############################################################################################
#############################################################################################
## Figures 4
#############################################################################################
#############################################################################################

## LSTM results
files <- list.files("LSTM_7", full.names = TRUE)
LSTM_list <- list()
for (i in 1:length(files)) {
  LSTM_list[[i]] <- t(as.matrix(read.csv(files[i])[,2:3]))
}

## TPFN results
files <- list.files("TPFN_7", full.names = TRUE)
TPFN_list <- list()
for (i in 1:length(files)) {
  TPFN_list[[i]] <- read.csv(files[i])[,3]
}

# statistical model results
sim_results <- get(load('prediction_sim_result_save_d7.RData'))


title_list <- list('Spherical AR(1)', 'Spherical AR(2)', 'Spherical AR(3)')


plot_list_fig1 <- list()


for (k in 1:3) {
  result_matrix <- 0
  for (j in 1:200) {
    result_matrix <- result_matrix + sim_results[[j]][[k]]
  }
  result_matrix <- result_matrix / 200
  
  df <- data.frame(
    index = rep(1:ncol(result_matrix), 6),
    angle = c(result_matrix[1,], result_matrix[2,], result_matrix[3,], 
              LSTM_list[[k]][1,], LSTM_list[[k]][2,], TPFN_list[[k]]),
    group = rep(c("SAR", "DSAR", "TPSAR", 'LSTM-S', 'LSTM-P', 'TimePFN'), each = ncol(result_matrix))
  )
  
  p <- ggplot(df, aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
    geom_point(size = 1.5, fill = NA, stroke = 1) +
    geom_line() +
    scale_color_manual(values = c("SAR" = "#009E73", "DSAR" = "#56B4E9", "TPSAR" = "#D55E00", "LSTM-S" = "#CC79A7", "LSTM-P" = "#0072B2","TimePFN" = "#E69F00")) +
    scale_linetype_manual(values = c("SAR" = "dashed", "DSAR" = "dotted", "TPSAR" = "solid", "LSTM-S" = "dotdash", "LSTM-P" = "longdash","TimePFN" = "twodash")) +
    scale_shape_manual(values = c("SAR" = 24, "DSAR" = 22, "TPSAR" = 21, "LSTM-S" = 23, "LSTM-P" = 25,"TimePFN" = 8)) +
    scale_x_continuous(limits = c(1, ncol(result_matrix)), breaks = scales::breaks_pretty(n = 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.6, "cm"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.margin = margin(5, 15, 5, 5)
    )
  
  if (k == 1) { 
    p <- p + scale_y_continuous(limits = c(0.2, 0.44), labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_text(size = 14))
  } else if (k == 3) { 
    p <- p + scale_y_continuous(
      limits = c(0.2, 0.44), labels = label_number(accuracy = 0.01),
      sec.axis = sec_axis(~., name = "T=120", breaks = NULL)
    ) + theme(axis.title.y.left = element_blank(), axis.title.y.right = element_text(size = 14, angle = 270))
  } else { 
    p <- p + scale_y_continuous(limits = c(0.2, 0.44), labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_blank())
  }
  
  p <- p + labs(
    title = title_list[[k]],
    x = "Number of Steps Ahead",
    y = if (k == 1) "Prediction Error" else ""
  )
  
  plot_list_fig1[[k]] <- p
}


fig1 <- plot_list_fig1[[1]] + plot_list_fig1[[2]] + plot_list_fig1[[3]] +
  plot_layout(nrow = 1, guides = 'collect') & 
  theme(legend.position = 'bottom', legend.key.width = unit(3, "cm"), legend.spacing.x = unit(0.5, 'cm'))

fig1

#############################################################################################
#############################################################################################
## Figures S3
#############################################################################################
#############################################################################################
plot_list_fig2 <- list()
for (k in 4:9) {
  result_matrix <- 0
  for (j in 1:200) {
    result_matrix <- result_matrix + sim_results[[j]][[k]]
  }
  result_matrix <- result_matrix / 200
  
  df <- data.frame(
    index = rep(1:ncol(result_matrix), 6),
    angle = c(result_matrix[1,], result_matrix[2,], result_matrix[3,], 
              LSTM_list[[k]][1,], LSTM_list[[k]][2,], TPFN_list[[k]]),
    group = rep(c("SAR", "DSAR", "TPSAR", 'LSTM-S', 'LSTM-P', 'TimePFN'), each = ncol(result_matrix))
  )
  
  i_row <- if (k <= 6) 1 else 2        
  i_ar <- (k - 1) %% 3 + 1             
  current_limits <- if (i_row == 1) c(0.2, 0.42) else c(0.2, 0.48)
  current_T <- if (i_row == 1) "300" else "600"
  

  p <- ggplot(df, aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
    geom_point(size = 1.5, fill = NA, stroke = 1) +
    geom_line() +
    scale_color_manual(values = c("SAR" = "#009E73", "DSAR" = "#56B4E9", "TPSAR" = "#D55E00", "LSTM-S" = "#CC79A7", "LSTM-P" = "#0072B2","TimePFN" = "#E69F00")) +
    scale_linetype_manual(values = c("SAR" = "dashed", "DSAR" = "dotted", "TPSAR" = "solid", "LSTM-S" = "dotdash", "LSTM-P" = "longdash","TimePFN" = "twodash")) +
    scale_shape_manual(values = c("SAR" = 24, "DSAR" = 22, "TPSAR" = 21, "LSTM-S" = 23, "LSTM-P" = 25,"TimePFN" = 8)) +
    scale_x_continuous(limits = c(1, ncol(result_matrix)), breaks = scales::breaks_pretty(n = 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.6, "cm"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.margin = margin(5, 15, 5, 5)
    )
  
  if (i_ar == 1) { 
    p <- p + scale_y_continuous(limits = current_limits, labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_text(size = 14))
  } else if (i_ar == 3) { 
    p <- p + scale_y_continuous(
      limits = current_limits, labels = label_number(accuracy = 0.01),
      sec.axis = sec_axis(~., name = paste0("T=", current_T), breaks = NULL)
    ) + theme(axis.title.y.left = element_blank(), axis.title.y.right = element_text(size = 14, angle = 270))
  } else { 
    p <- p + scale_y_continuous(limits = current_limits, labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_blank())
  }
  
  p <- p + labs(
    title = if (i_row == 1) title_list[[i_ar]] else NULL,
    x = if (i_row == 2) "Number of Steps Ahead" else "",
    y = if (i_ar == 1) "Prediction Error" else ""
  )
  
  plot_list_fig2[[k - 3]] <- p
}

fig2 <- plot_list_fig2[[1]] + plot_list_fig2[[2]] + plot_list_fig2[[3]] +
  plot_list_fig2[[4]] + plot_list_fig2[[5]] + plot_list_fig2[[6]] +
  plot_layout(nrow = 2, guides = 'collect') & 
  theme(legend.position = 'bottom', legend.key.width = unit(3, "cm"), legend.spacing.x = unit(0.5, 'cm'))

fig2

#############################################################################################
#############################################################################################
## Figures S4
#############################################################################################
#############################################################################################

## LSTM results
files <- list.files("LSTM_48", full.names = TRUE)
LSTM_list <- list()
for (i in 1:length(files)) {
  LSTM_list[[i]] <- t(as.matrix(read.csv(files[i])[,2:3]))
}

## TPFN results
files <- list.files("TPFN_48", full.names = TRUE)
TPFN_list <- list()
for (i in 1:length(files)) {
  TPFN_list[[i]] <- read.csv(files[i])[,3]
}

## statistical model results
sim_results <- get(load('prediction_sim_result_save_d48.RData'))


sample_size_set <- c(120, 300, 600)
title_list <- list('Spherical AR(1)', 'Spherical AR(2)', 'Spherical AR(3)')


y_limits_list <- list(c(0.2, 0.44), c(0.2, 0.42), c(0.2, 0.48)) 

data_frame_list <- list()
plot_list <- list()
total_plots <- 9


for (k in 1:total_plots) {
  
  
  i_sample <- (k - 1) %/% 3 + 1  
  i_ar <- (k - 1) %% 3 + 1       
  
  result_matrix <- 0
  for (j in 1:200) {
    result_matrix <- result_matrix + sim_results[[j]][[k]]
  }
  result_matrix <- result_matrix / 200
  
  # Build DataFrame
  data_frame_list[[k]] <- data.frame(
    index = rep(1:ncol(result_matrix), 6),
    angle = c(result_matrix[1,], result_matrix[2,], result_matrix[3,], 
              LSTM_list[[k]][1,], LSTM_list[[k]][2,], TPFN_list[[k]]),
    group = rep(c("SAR", "DSAR", "TPSAR", 'LSTM-S', 'LSTM-P', 'TimePFN'), each = ncol(result_matrix))
  )
  
  # Base ggplot object
  p <- ggplot(data_frame_list[[k]], aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
    geom_point(size = 1.5, fill = NA, stroke = 1) +
    geom_line() +
    scale_color_manual(values = c("SAR" = "#009E73", "DSAR" = "#56B4E9", "TPSAR" = "#D55E00", "LSTM-S" = "#CC79A7", "LSTM-P" = "#0072B2","TimePFN" = "#E69F00")) +
    scale_linetype_manual(values = c("SAR" = "dashed", "DSAR" = "dotted", "TPSAR" = "solid", "LSTM-S" = "dotdash", "LSTM-P" = "longdash","TimePFN" = "twodash")) +
    scale_shape_manual(values = c("SAR" = 24, "DSAR" = 22, "TPSAR" = 21, "LSTM-S" = 23, "LSTM-P" = 25,"TimePFN" = 8)) +
    scale_x_continuous(limits = c(1, ncol(result_matrix)), breaks = scales::breaks_pretty(n = 5)) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.6, "cm"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_rect(fill = NA, color = "black"),
      plot.margin = margin(5, 15, 5, 5)
    )
  
  # Apply dynamic Y-axis scales, Secondary axes, and Labels based on position
  current_limits <- y_limits_list[[i_sample]]
  
  if (i_ar == 1) { 
    # FIRST COLUMN
    p <- p + scale_y_continuous(limits = current_limits, labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_text(size = 14))
    
  } else if (i_ar == 3) { 
    # LAST COLUMN 
    p <- p + scale_y_continuous(
      limits = current_limits, 
      labels = label_number(accuracy = 0.01),
      sec.axis = sec_axis(~., name = paste0("T=", sample_size_set[i_sample]), breaks = NULL)
    ) + theme(
      axis.title.y.left = element_blank(),
      axis.title.y.right = element_text(size = 14, angle = 270)
    )
    
  } else { 
    # MIDDLE COLUMN
    p <- p + scale_y_continuous(limits = current_limits, labels = label_number(accuracy = 0.01)) +
      theme(axis.title.y = element_blank())
  }
  
  p <- p + labs(
    title = if (i_sample == 1) title_list[[i_ar]] else NULL,
    x = if (i_sample == 3) "Number of Steps Ahead" else "",
    y = if (i_ar == 1) "Prediction Error" else ""
  )
  
  plot_list[[k]] <- p
}


combined_plot <- plot_list[[1]] + plot_list[[2]] + plot_list[[3]] +
  plot_list[[4]] + plot_list[[5]] + plot_list[[6]] +
  plot_list[[7]] + plot_list[[8]] + plot_list[[9]] +
  plot_layout(nrow = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(3, "cm"), 
    legend.spacing.x = unit(0.5, 'cm')
  )

combined_plot