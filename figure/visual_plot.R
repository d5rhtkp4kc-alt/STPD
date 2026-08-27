## Description: figures not based on simulation results or real data analysis results
## reproduce Figures 4, S3, and S4.

library(movMF)
library(ggplot2)
library(ggforce)
library(patchwork)

source('../main function/basic_function.R')

#############################################################################################
#############################################################################################
## Figure 1 
#############################################################################################
#############################################################################################

AR_list <- AR_generate_function_plot(seed_set=1,sample_size=200,d=7,coef=c(0.3,-0.1,0.4))


sphere_matrix <- matrix(NA,nrow = 200,ncol = 7)
for (t in 1:200) {
  mu <- c(1:7) / sqrt(sum((c(1:7))^2))
  angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
  sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
}
period_matrix <- period_generate_funtion(AR_list)
periodic_function_matrix <- periodic_component_function(sample.size,d=7,period=12)
periodic_function_matrix_full <- matrix(NA,nrow = 200,ncol = 7)
for (k in 1:200) {
  current_ind <- k + 12 - floor((k + 12 - 1)/12) * 12
  periodic_function_matrix_full[k,] <- periodic_function_matrix[current_ind,]
}
trend_matrix <- trend_generate_funtion(period_matrix)



prepare_plot_data <- function(mat, x_axis_vals) {
  df <- as.data.frame(mat)
  colnames(df) <- paste0("Coordinate ", 1:7)
  df$x <- x_axis_vals
  return(tidyr::pivot_longer(df, cols = starts_with("Coord"), names_to = "Coordinate", values_to = "value"))
}


custom_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "11")

create_styled_plot <- function(data, x_lab, y_lab, title) {
  ggplot(data, aes(x = x, y = value, color = Coordinate, linetype = Coordinate)) +
    geom_line(linewidth = 0.7) +
    scale_linetype_manual(values = custom_linetypes) +
    labs(y = y_lab, x = x_lab, title = title) +
    theme_classic()+
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          legend.text=element_text(size = 14),
          plot.margin=margin(5,15,5,5),
          panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom",plot.title = element_text(hjust = 0,size=14))
}


y_matrix <- trend_matrix
sample_size <- 200
x_input <- matrix(c(1:sample_size)/sample_size,nrow = sample_size, ncol = 1) 
trend_component <- LocSpheReg(xin = as.vector(x_input), yin = y_matrix)
trend_component <- trend_component$yout



mu <- intrinsic_mean_sphere(y_matrix)
remove_trend_result <- get_spherical_residuals(y_matrix, trend_component,mu = mu)
remove_trend_data <- remove_trend_result$residuals


p=12
x <- covariate_function(p,sample_size)
periodic_component <- GloSpheGeoReg(xin = x, yin = remove_trend_data)
final_result <- get_spherical_residuals(remove_trend_data, periodic_component)
final_resid_data <- final_result$residuals

t_vals <- 1:200
u_vals <- t_vals / 200
# Plot (1): 
p1 <- create_styled_plot(prepare_plot_data(y_matrix, t_vals), 
                         "", 'Value', "(a) Original simulated time series")

# Plot (2): 
p2 <- create_styled_plot(prepare_plot_data(remove_trend_data, t_vals), 
                         "", 'Value', "(b) De-trended time series")

# Plot (3):
p3 <- create_styled_plot(prepare_plot_data(final_resid_data, t_vals), 
                         "Observation index", 'Value', "(c) De-trended and de-seasonalized time series")


combined_plot <- p1 + p2 + p3 +
  plot_layout(nrow = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    # Increase the width of the legend lines
    legend.key.width = unit(3, "line"), 
    # Optional: add a bit of space between legend items
    legend.spacing.x = unit(0.5, 'cm')
  )


print(combined_plot)

#############################################################################################
#############################################################################################
## Figure S1
#############################################################################################
#############################################################################################
remove_trend_result <- y_matrix-trend_component
row_norms <- sqrt(rowSums(remove_trend_result^2))
remove_trend_data <- remove_trend_result / row_norms

final_result <- remove_trend_data-periodic_component
row_norms <- sqrt(rowSums(final_result^2))
final_resid_data <- final_result/row_norms

# Plot (1): 
p1 <- create_styled_plot(prepare_plot_data(y_matrix, t_vals), 
                         "", 'Value', "(a) Original simulated time series")

# Plot (2): 
p2 <- create_styled_plot(prepare_plot_data(remove_trend_data, t_vals), 
                         "", 'Value', "(b) De-trended time series")

# Plot (3):
p3 <- create_styled_plot(prepare_plot_data(final_resid_data, t_vals), 
                         "Observation index", 'Value', "(c) De-trended and de-seasonalized time series")


combined_plot <- p1 + p2 + p3 +
  plot_layout(nrow = 3, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(3, "line"), 
    legend.spacing.x = unit(0.5, 'cm')
  )

print(combined_plot)

#############################################################################################
#############################################################################################
## Figure 2
#############################################################################################
#############################################################################################


title_style <- element_text(hjust = 0.5, size = 14, margin = margin(t = 10, b = 20))
common_theme <- theme_void() + 
  theme(
    plot.title = title_style,
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 40, b = 10, l = 40)
  )

# --- Euclidean Space Plot ---
df_euclid <- data.frame(
  x = c(0, 3, 1, 2), 
  y = c(0, 2, 3, -1),
  label = c("italic(c)", "italic(a)", "italic(b)", "italic(d)")
)

p1 <- ggplot(df_euclid) +
  geom_vline(xintercept = 0, color = "gray80", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "gray80", linewidth = 0.5) +
  geom_segment(aes(x = 1, y = 3, xend = 3, yend = 2), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               linewidth = 1.2, color = "firebrick") +
  geom_segment(aes(x = 0, y = 0, xend = 2, yend = -1), 
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"), 
               linewidth = 1.2, color = "firebrick") +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 3), linetype = "dashed", color = "blue", linewidth = 1) +
  geom_segment(aes(x = 2, y = -1, xend = 3, yend = 2), linetype = "dashed", color = "blue", linewidth = 1) +
  geom_point(aes(x, y), size = 3, color = "firebrick") +
  geom_text(aes(x, y, label = label), vjust = -1.2, size = 6, parse = TRUE, fontface = "bold.italic") +
  coord_fixed() +
  expand_limits(y = c(-1.5, 4), x = c(-0.5, 3.5)) +
  labs(title = "(a) Optimal transport in Euclidean space") +
  common_theme

# --- Helper Functions ---
norm_vec <- function(x) { sqrt(sum(x^2)) }
geo_dist <- function(x, y) { acos(max(min(sum(x * y), 1), -1)) }

get_rotation_matrix <- function(g_from, g_to) {
  p <- length(g_from)
  theta <- geo_dist(g_from, g_to)
  if (theta < 1e-10) return(diag(p))
  u1 <- g_from
  proj <- g_to - sum(g_to * u1) * u1
  u2 <- proj / norm_vec(proj)
  Q <- (u2 %*% t(u1)) - (u1 %*% t(u2))
  R <- diag(p) + sin(theta) * Q + (1 - cos(theta)) * (Q %*% Q)
  return(R)
}

get_smooth_geodesic <- function(p1, p2, n = 500) {
  p1 <- as.numeric(p1 / norm_vec(p1))
  p2 <- as.numeric(p2 / norm_vec(p2))
  theta <- acos(sum(p1 * p2))
  if (theta < 1e-10) return(data.frame(V1=rep(p1[1],n), V2=rep(p1[2],n), V3=rep(p1[3],n)))
  t <- seq(0, 1, length.out = n)
  path <- sapply(t, function(ti) { (sin((1-ti)*theta)*p1 + sin(ti*theta)*p2)/sin(theta) })
  return(as.data.frame(t(path)))
}

# --- Spherical Plot  ---
g1_3d <- c(cos(2.2), 0, sin(2.2)) 
g3_3d <- c(cos(3.1), 0, sin(3.1))  
q2 <- c(-0.1, -0.6, 0.5)
g2_3d <- as.numeric(q2 / norm_vec(q2)) 
R_mat <- get_rotation_matrix(g1_3d, g3_3d)
g4_3d <- drop(R_mat %*% g2_3d)



path_red_top <- get_smooth_geodesic(g1_3d, g2_3d) # b to a
path_red_bot <- get_smooth_geodesic(g3_3d, g4_3d) # c to d
path_blue1   <- get_smooth_geodesic(g1_3d, g3_3d) # b to c
path_blue2   <- get_smooth_geodesic(g2_3d, g4_3d) # a to d

t_eq <- seq(0, pi, length.out = 1000)
eq_front <- data.frame(x = cos(t_eq + pi), y = 0.3 * sin(t_eq + pi))
eq_back  <- data.frame(x = cos(t_eq), y = 0.3 * sin(t_eq))


p2 <- ggplot() +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), color = "black", linewidth = 0.5) +
  geom_path(data = eq_back, aes(x = x, y = y), linetype = "dashed", color = "gray75", linewidth = 0.4) +
  geom_path(data = eq_front, aes(x = x, y = y), color = "gray40", linewidth = 0.4) +
  geom_path(data = path_blue1, aes(x = V1, y = V3), color = "blue", linewidth = 1, linetype = "dashed") +
  geom_path(data = path_blue2, aes(x = V1, y = V3), color = "blue", linewidth = 1, linetype = "dashed") +
  geom_path(data = path_red_top, aes(x = V1, y = V3), color = "firebrick", linewidth = 1.2,
            arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  geom_path(data = path_red_bot, aes(x = V1, y = V3), color = "firebrick", linewidth = 1.2,
            arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  geom_point(aes(x = c(g1_3d[1], g2_3d[1], g3_3d[1], g4_3d[1]), 
                 y = c(g1_3d[3], g2_3d[3], g3_3d[3], g4_3d[3])), color = "firebrick", size = 2) +
  annotate("text", x = g1_3d[1]-0.15, y = g1_3d[3]+0.05, label = "italic(b)", parse=T, size=6) +
  annotate("text", x = g2_3d[1]+0.15, y = g2_3d[3], label = "italic(a)", parse=T, size=6) +
  annotate("text", x = g3_3d[1]-0.15, y = g3_3d[3]-0.05, label = "italic(c)", parse=T, size=6) +
  annotate("text", x = g4_3d[1]+0.15, y = g4_3d[3], label = "italic(d)", parse=T, size=6) +
  labs(title = "(b) Optimal transport in sphere") +
  coord_fixed(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1)) +
  common_theme


p1 + p2

#############################################################################################
#############################################################################################
## Figure 3
#############################################################################################
#############################################################################################

AR_list <- AR_generate_function_plot(seed_set=1,sample_size=200,d=7,coef=c(0.3,-0.1,0.4))


sphere_matrix <- matrix(NA,nrow = 200,ncol = 7)
for (t in 1:200) {
  mu <- c(1:7) / sqrt(sum((c(1:7))^2))
  angle_n <- sqrt(sum((AR_list[[t]] %*% mu)^2))
  sphere_matrix[t,] <- exp_from_log_rotation_predict(AR_list[[t]] ,angle_n) %*% mu
}
period_matrix <- period_generate_funtion(AR_list)
periodic_function_matrix <- periodic_component_function(sample.size,d=7,period=12)
periodic_function_matrix_full <- matrix(NA,nrow = 200,ncol = 7)
for (k in 1:200) {
  current_ind <- k + 12 - floor((k + 12 - 1)/12) * 12
  periodic_function_matrix_full[k,] <- periodic_function_matrix[current_ind,]
}
trend_matrix <- trend_generate_funtion(period_matrix)


u_grid <- (1:200) / 200

trend_function_matrix <- t(sapply(u_grid, trend_component_function,
                                  d=7,
                                  inc_idx = c(1,2,3),
                                  A = 0.6))

## ggplot
t_vals <- 1:200
u_vals <- t_vals / 200

# Plot (1): g(t) vs t
p1 <- create_styled_plot(prepare_plot_data(periodic_function_matrix_full, t_vals), 
                         "t", expression(g(t)), "(a) Periodic component")

# Plot (2): f(u) vs u
p2 <- create_styled_plot(prepare_plot_data(trend_function_matrix, u_vals), 
                         "u", expression(f(u)), "(b) Trend component")

# Plot (3): R_t^(2) vs t
p3 <- create_styled_plot(prepare_plot_data(sphere_matrix, t_vals), 
                         "t", expression(R[t]^{(2)}), "(c) Spherical AR")

# Plot (4): R_t^(1) vs t
p4 <- create_styled_plot(prepare_plot_data(period_matrix, t_vals), 
                         "t", expression(R[t]^{(1)}), "(d) Spherical AR + periodic component")

# Plot (5): y_t vs t
p5 <- create_styled_plot(prepare_plot_data(trend_matrix, t_vals), 
                         "t", expression(y[t]), "(e)  Spherical AR + periodic and trend components")


layout <- "
AAABBB
CCDDEE
"


combined_plot <- p1 + p2 + p3 + p4 + p5 + 
  plot_layout(design = layout, guides = 'collect') & 
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(3, "line"), 
    legend.spacing.x = unit(0.5, 'cm')
  )

# Display the result
print(combined_plot)

