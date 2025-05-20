# Load required libraries
library(ggplot2)
library(lme4)
library(ggeffects)
library(dplyr)
library(plotly)


# --- Data Loading and Preparation ---
data <- read.csv2("path_to_csv.csv")
data$log_KDE <- log1p(as.numeric(data$KDE))   # log(1 + KDE) 


# --- Logistic Regression Models ---
# Nest vs. Random
glm_nr_epa <- glm(factor(Point == "nest") ~ dist_water * log_KDE + dist_veg,
                  data = data %>% filter(Point != "regular"),
                  family = binomial)
summary(glm_nr_epa)


# Nest vs. Regular
glm_nreg_epa <- glm(factor(Point == "nest") ~ dist_water * log_KDE + dist_veg,
                    data = data %>% filter(Point != "random"),
                    family = binomial)
summary(glm_nreg_epa)


# 3D PLOTS

calculate_and_plot_three_planes <- function(model, data) { 
  # Define the three levels of log_epa
  log_KDE_levels <- quantile(data$log_KDE, probs = c(0.5, 0.75, 1), na.rm = TRUE)
  names(log_KDE_levels) <- c("Low", "Medium", "High")
  
  # Create grid of dist_water and dist_veg
  dist_water_seq <- seq(min(data$dist_water, na.rm = TRUE), max(data$dist_water, na.rm = TRUE), length.out = 100)
  dist_veg_seq <- seq(min(data$dist_veg, na.rm = TRUE), max(data$dist_veg, na.rm = TRUE), length.out = 100)
  
  # Create base grid
  base_grid <- expand.grid(dist_water = dist_water_seq, dist_veg = dist_veg_seq)
  
  # Start plotly figure
  fig <- plot_ly()
  
  # Loop through log_epa levels and add surfaces
  for (level in names(log_KDE_levels)) {
    new_data <- base_grid %>%
      mutate(log_epa = log_KDE_levels[[level]])
    
    # Predict
    preds <- predict(model, newdata = new_data, type = "response")
    new_data$predicted_probability <- preds
    
    # Reshape predictions into matrix for surface
    z_matrix <- matrix(preds, 
                       nrow = length(dist_water_seq), 
                       ncol = length(dist_veg_seq))
    
    # Add surface to plot
    fig <- fig %>% add_surface(
      x = dist_veg_seq,
      y = dist_water_seq,
      z = z_matrix,
      showscale = FALSE,
      name = paste("log_KDE:", level),
      colorscale = switch(level,
                          "Low" = list(c(0, '#1b9e77'), c(1, '#1b9e77')),
                          "Medium" = list(c(0, '#d95f02'), c(1, '#d95f02')),
                          "High" = list(c(0, '#7570b3'), c(1, '#7570b3')))
    )
  }
  
  # Final layout
  fig <- fig %>% layout(
    title = "Predicted Nesting Probability",
    scene = list(
      xaxis = list(title = "Distance to veg"),
      yaxis = list(title = "Distance to water"),
      zaxis = list(title = "Predicted Probability")
    )
  )
  
  return(fig)
}
overall_plot <- calculate_and_plot_three_planes(model, data)
overall_plot
htmlwidgets::saveWidget(overall_plot, "3D_PLOT.html")



# 2D STATIC PLOT

plot_static_all_beaches_probability <- function(model, data) {
  # Create prediction grid for all beaches combined
  log_KDE_levels <- quantile(data$log_KDE, probs = c(0.5, 0.75, 1), na.rm = TRUE)
  names(log_KDE_levels) <- c("Low", "Medium", "High")
  
  # Create prediction grid for all combinations of distance to water, vegetation, and log_epa
  new_data_all <- expand.grid(
    dist_water = seq(min(data$dist_water), max(data$dist_water), length.out = 100),
    dist_veg = seq(min(data$dist_veg), max(data$dist_veg), length.out = 100),
    log_KDE = log_KDE_levels
  )
  
  # Predict probabilities for all combinations of variables
  predictions_all <- predict(model, newdata = new_data_all, type = "response")
  plot_data_all <- cbind(new_data_all, predicted_probability = predictions_all)
  
  # Plot
  ggplot(plot_data_all, aes(x = dist_water, y = dist_veg, fill = predicted_probability)) +
    geom_tile() +
    facet_wrap(~ log_KDE, labeller = label_bquote(log_KDE == .(log_KDE))) +
    scale_fill_gradientn(
      colors = my_color_ramp(100),   # Generate 100 interpolated colors
      name = "Predicted\nProbability",
      limits = c(0, 1),              # Ensures scale goes from 0 to 1
      oob = scales::squish           # Prevents out-of-bounds values from throwing errors
    ) +
    labs(
      title = "Predicted Nesting Probability Across All Beaches",
      x = "Distance to Water",
      y = "Distance to Vegetation"
    ) +
    theme_minimal()
}

plot_static_all_beaches_probability(model, data)



