## reproduce Figures S5 and S11.

library(ggplot2)
library(patchwork)
library(scales)

#############################################################################################
#############################################################################################
## Figures S5
#############################################################################################
#############################################################################################

prediction_error_matrix <- get(load('US_energy_pred_30.RData'))
LSTM_result <- as.matrix(read.csv('energy_LSTM.csv')[,2:31])
timePFN_result <- read.csv('timepfn_energy_results.csv')[,2]

df <- data.frame(
  index = rep(1:ncol(prediction_error_matrix), 6),
  angle = c(prediction_error_matrix[1,], prediction_error_matrix[2,],prediction_error_matrix[3,],LSTM_result[1,],LSTM_result[2,],timePFN_result),
  group = rep(c("SAR", "DSAR","TPSAR",'LSTM-S','LSTM-P','TimePFN'), each = ncol(prediction_error_matrix))
)



df$method_type <- ifelse(
  df$group %in% c("SAR", "DSAR", "TPSAR", "TimePFN"), 
  "AR Models and TimePFN", 
  "LSTM Models"
)


all_groups <- c("SAR", "DSAR", "TPSAR", "LSTM-S", "LSTM-P", "TimePFN")

shared_layers <- list(
  scale_color_manual(
    values = c("SAR" = "#009E73", "DSAR" = "#56B4E9", "TPSAR" = "#D55E00", "LSTM-S" = "#CC79A7", "LSTM-P" = "#0072B2","TimePFN" = "#E69F00"), 
    limits = all_groups
  ),
  scale_linetype_manual(
    values = c("SAR" = "dashed", "DSAR" = "dotted", "TPSAR" = "solid", "LSTM-S" = "dotdash", "LSTM-P" = "longdash","TimePFN" = "twodash"), 
    limits = all_groups
  ),
  scale_shape_manual(
    values = c("SAR" = 24, "DSAR" = 22, "TPSAR" = 21, "LSTM-S" = 23, "LSTM-P" = 25,"TimePFN" = 8), 
    limits = all_groups
  ),
  scale_x_continuous(limits = c(1, ncol(prediction_error_matrix))),
  theme_classic(),
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.title = element_text(size = 14, hjust = 0.5)
  )
)


df_ar <- df
df_ar$angle[df_ar$method_type == "LSTM Models"] <- NA

df_lstm <- df
df_lstm$angle[df_lstm$method_type == "AR Models and TimePFN"] <- NA


# Plot 1 (AR Methods)
p1 <- ggplot(df_ar, aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
  # Added na.rm = TRUE so R doesn't warn you about the missing NA values we just created
  geom_point(size = 1.5, fill = NA, stroke = 1, na.rm = TRUE) + 
  geom_line(na.rm = TRUE) +
  labs(title = "AR Models and TimePFN", x = "Number of Steps Ahead", y = "Prediction Error") +
  shared_layers

# Plot 2 (LSTM Methods)
p2 <- ggplot(df_lstm, aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
  geom_point(size = 1.5, fill = NA, stroke = 1, na.rm = TRUE) + 
  geom_line(na.rm = TRUE) +
  labs(title = "LSTM Models", x = "Number of Steps Ahead", y = "Prediction Error") +
  shared_layers


p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    legend.key.width = unit(3, "cm")
  )


#############################################################################################
#############################################################################################
## Figures S11
#############################################################################################
#############################################################################################

prediction_error_matrix <- get(load('transport_pred_30.RData'))
LSTM_result <- as.matrix(read.csv('trip_LSTM.csv')[,2:31])
timePFN_result <- read.csv('timepfn_trip_results.csv')[,2]

df <- data.frame(
  index = rep(1:ncol(prediction_error_matrix), 6),
  angle = c(prediction_error_matrix[1,], prediction_error_matrix[2,],prediction_error_matrix[3,],LSTM_result[1,],LSTM_result[2,],timePFN_result),
  group = rep(c("SAR", "DSAR","TPSAR",'LSTM-S','LSTM-P','TimePFN'), each = ncol(prediction_error_matrix))
)



ggplot(df, aes(x = index, y = angle, color = group, shape = group, linetype = group)) +
  geom_point(size = 1.5, fill = NA, stroke = 1) +
  geom_line() +
  scale_color_manual(values = c(
    "SAR"   = "#009E73",  
    "DSAR"  = "#56B4E9", 
    "TPSAR" = "#D55E00",
    "LSTM-S" = "#CC79A7", # Added: Purplish pink
    "LSTM-P" = "#0072B2",  # Added: Dark blue
    "TimePFN" = "#E69F00"
  )) +
  scale_linetype_manual(values = c(
    "SAR"   = "dashed",  
    "DSAR"  = "dotted", 
    "TPSAR" = "solid",
    "LSTM-S" = "dotdash", # Added: Dot-dash line
    "LSTM-P" = "longdash", # Added: Long dashed line
    "TimePFN" = "twodash"
  )) +
  scale_shape_manual(values = c(
    "SAR"   = 24, # Triangle point up
    "DSAR"  = 22, # Square
    "TPSAR" = 21, # Circle
    "LSTM-S" = 23, # Added: Diamond
    "LSTM-P" = 25,  # Added: Triangle point down
    "TimePFN" = 8
  )) +
  scale_x_continuous(limits = c(1, ncol(prediction_error_matrix))) + 
  labs(
    x = "Number of Steps Ahead",
    y = "Prediction Error"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 14),
    legend.key.width = unit(3, "cm"), 
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.margin = margin(5, 15, 5, 5)
  )

