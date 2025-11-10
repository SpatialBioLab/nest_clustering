# Load required libraries
library(ggplot2)
library(dplyr)

# 2D STATIC PLOT

plot_static_kernel_effect_gam <- function(model, data) {
  library(ggplot2)
  library(viridis)
  
  kernel_mean_val <- mean(data$kernel_mean_z, na.rm = TRUE)
  ref_praia <- unique(data$Praia)[1]
  
  mean_kernel <- mean(data$kernel_within_z, na.rm = TRUE)
  sd_kernel   <- sd(data$kernel_within_z, na.rm = TRUE)
  
  kernel_levels <- c(q[1], q[2], q[3], q[4])
  names(kernel_levels) <- c("Low", "Medium", "High", "Hotspot")
  
  dist_seq <- seq(min(data$dist_z, na.rm = TRUE), max(data$dist_z, na.rm = TRUE), length.out = 100)
  veg_seq  <- seq(min(data$dist_veget_z, na.rm = TRUE), max(data$dist_veget_z, na.rm = TRUE), length.out = 100)
  
  plot_data_all <- do.call(rbind, lapply(names(kernel_levels), function(level) {
    new_grid <- expand.grid(dist_z = dist_seq, dist_veget_z = veg_seq)
    new_grid$kernel_within_z <- kernel_levels[[level]]
    new_grid$kernel_mean_z   <- kernel_mean_val
    new_grid$surveyID <- ref_praia
    new_grid$kernel_level <- factor(level, levels = c("Low", "Medium", "High", "Hotspot"))
    
    new_grid$predicted_probability <- predict(model, newdata = new_grid, type = "response")
    new_grid
  }))
  
  ggplot(plot_data_all, aes(x = dist_z, y = dist_veget_z, fill = predicted_probability)) +
    geom_tile() +
    facet_wrap(~ kernel_level, nrow = 1) +
    scale_fill_viridis_c(name = "Predicted\nprobability") +
    labs(
      title = "GAM: Nests vs. random points",
      x = "Distance to water (z-score)",
      y = "Distance to vegetation (z-score)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_static_kernel_effect_gam(mod_gam, dados)