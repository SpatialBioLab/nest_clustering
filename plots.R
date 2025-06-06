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
  log_KDE_levels <- quantile(data$log_KDE, probs = c(0.5, 0.7, 0.9), na.rm = TRUE)
  
  p50 <- log_epa_levels[[1]]
  p70 <- log_epa_levels[[2]]
  p90 <- log_epa_levels[[3]]
  
  # Step 2: Classify each point into a category
  data$accum_cat <- cut(data$log_epa,
                        breaks = c(-Inf, p50, p70, p90, Inf),
                        labels = c("Low", "Medium", "High", "Hotspot"),
                        include.lowest = TRUE)
   
  mean_log_epa <- tapply(data$log_epa, data$accum_cat, mean, na.rm = TRUE)
  
  # Create grid of predictors
  dist_water_seq <- seq(min(data$dist_water, na.rm = TRUE),
                        max(data$dist_water, na.rm = TRUE), length.out = 100)
  dist_veg_seq <- seq(min(data$dist_veg, na.rm = TRUE),
                      max(data$dist_veg, na.rm = TRUE), length.out = 100)
  base_grid <- expand.grid(dist_water = dist_water_seq, dist_veg = dist_veg_seq)
  
  # Start plotly figure
  fig <- plot_ly()
  
  # Loop through log_epa levels and add surfaces
  for (cat in names(mean_log_epa)) {
    new_data <- base_grid %>%
      mutate(log_epa = mean_log_epa[[cat]])
    
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
      name = cat,
      colorscale = switch(cat,
                          "Low" = list(c(0, '#1b9e77'), c(1, '#1b9e77')),
                          "Medium" = list(c(0, '#d95f02'), c(1, '#d95f02')),
                          "High" = list(c(0, '#7570b3'), c(1, '#7570b3')),
                          "Hotspot" = list(c(0, '#2d004b'), c(1, '#2d004b)))
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
  log_epa_levels <- quantile(data$log_epa, probs = c(0.5, 0.7, 0.9), na.rm = TRUE)
  
  p50 <- log_epa_levels[[1]]
  p70 <- log_epa_levels[[2]]
  p90 <- log_epa_levels[[3]]

  data$accum_cat <- cut(data$log_epa,
    breaks = c(-Inf, p50, p70, p90, Inf),
    labels = c("Low (<50%)", "Medium (50–70%)", "High (70–90%)", "Hotspot (>90%)"),
    include.lowest = TRUE)

  log_epa_means <- data %>%
    group_by(accum_cat) %>%
    summarize(log_epa = mean(log_epa, na.rm = TRUE)) %>%
    na.omit()

  # Create prediction grid for all combinations of distance to water, vegetation, and log_epa
  dist_water_seq <- seq(min(data$dist_water), max(data$dist_water), length.out = 100)
  dist_veg_seq <- seq(min(data$dist_veg), max(data$dist_veg), length.out = 100)
  
  plot_data_all <- do.call(rbind, lapply(1:nrow(log_epa_means), function(i) {
    level <- log_epa_means$accum_cat[i]
    value <- log_epa_means$log_epa[i]
    
    new_grid <- expand.grid(
      dist_water = dist_water_seq,
      dist_veg = dist_veg_seq
    ) %>%
      mutate(log_epa = value, accum_cat = level)
    
    preds <- predict(model, newdata = new_grid, type = "response")
    new_grid$predicted_probability <- preds
    return(new_grid)
  }))

  # Plot
    ggplot(plot_data_all, aes(x = dist_water, y = dist_veg, fill = predicted_probability)) +
    geom_tile() +
    facet_wrap(~ accum_cat, nrow = 1) +  # <-- This makes 4 plots in one row
    scale_fill_viridis_c(name = "Predicted\nProbability") +
    labs(
      title = "Predicted Nesting Probability by Accumulation Range",
      x = "Distance to Water",
      y = "Distance to Vegetation"
    ) +
    theme_minimal()
  
}

plot_static_all_beaches_probability(model, data)



