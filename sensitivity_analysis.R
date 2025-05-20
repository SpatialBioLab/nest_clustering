# SENSITIVITY ANALYSIS

# Load necessary libraries
library(dplyr)      # For data manipulation
library(glmnet)     # For GLM models (if you used this)
library(MASS)       # For mvrnorm (if needed for simulation)
library(tidyverse)

# load the points of nests/random/regular with associated distance to waterline, distance to vegetation, and kernel of marine debris
data <- read.csv2("path_to_points.csv")
data$log_KDE <- log1p(as.numeric(data$KDE))  # equivalent to log(1 + x)

  
# 1.  Define the function to add random error to the nest coordinates
add_error <- function(data, error_sd = 3) {
  #' Adds random error to nest coordinates (dist_water, dist_veg).
  #'
  #' @param data A data frame containing nest location data, including
  #'             columns for 'dist_water' and 'dist_veg'.
  #' @param error_sd The standard deviation of the random error
  #'                 (default: 3 meters).
  #' @return A data frame with added random error to 'dist_water' and
  #'          'dist_veg'.
  
  # Create new columns with simulated data
  dist_water_sim <- as.numeric(data$dist_water) + rnorm(nrow(data), mean = 0, sd = error_sd)
  dist_veg_sim <- as.numeric(data$dist_veg) + rnorm(nrow(data), mean = 0, sd = error_sd)
  
  # Ensure distances don't go negative
  dist_water_sim <- pmax(dist_water_sim, 0)
  dist_veg_sim <- pmax(dist_veg_sim, 0)
  
  # Add the new columns to the data frame
  data$dist_water_sim <- dist_water_sim
  data$dist_veg_sim <- dist_veg_sim
  return(data)
}


# 2.  Function to fit the GLM model
fit_glm_model <- function(data, comparison_group = "random") {
  #' Fits the GLM model predicting nest presence.
  #'
  #' @param data A data frame containing the data, including original and simulated
  #'             distances ('dist_water', 'dist_veg', 'dist_water_sim', 'dist_veg_sim'),
  #'             log transformed marine debris density ('log_epa' or 'log_tri')
  #'             and nest presence ('Point').  Point should have levels "nest",
  #'             "regular", and "random".
  #' @param comparison_group  The group to compare "nest" against.  Should be
  #'                          either "random" or "regular".
  #' @return The fitted GLM model object.
  if (!comparison_group %in% c("random", "regular")) {
    stop("comparison_group must be either 'random' or 'regular'")
  }
  
  # Create binary variable
  is_nest <- ifelse(data$Point == "nest", 1, 0)
  
  # Add the binary variable to the data frame
  data$is_nest <- is_nest
  
  model <- glm(
    is_nest ~ dist_water_sim * log_KDE + dist_veg_sim,
    family = binomial,
    data = data %>% filter(Point %in% c("nest", comparison_group))
  )
  return(model)
}

# 3. Function to extract and compare coefficients
extract_coefficients <- function(model) {
  #' Extracts coefficients from a GLM model.
  #'
  #' @param model A fitted GLM model object.
  #' @return A data frame with the coefficients, their standard errors,
  #'         and p-values.
  coef_summary <- summary(model)$coefficients
  coef_df <- as.data.frame(coef_summary)
  coef_df$term <- rownames(coef_df)
  return(coef_df)
}

# 4.  Perform the sensitivity analysis
sensitivity_analysis <- function(original_data, n_simulations = 100, comparison_group = "random") {
  
  # 1. Adiciona erro nos dados originais
  original_data_sim <- add_error(original_data)
  
  # 2. Ajusta o modelo com os dados simulados originais
  original_model <- fit_glm_model(original_data_sim, comparison_group = comparison_group)
  original_coefs <- extract_coefficients(original_model)
  
  simulation_coefs <- list()
  
  for (i in 1:n_simulations) {
    simulated_data <- add_error(original_data)
    simulated_model <- fit_glm_model(simulated_data, comparison_group = comparison_group)
    simulation_coefs[[i]] <- extract_coefficients(simulated_model)
  }
  
  all_sim_coefs <- bind_rows(simulation_coefs, .id = "simulation")
  
  return(list(
    original_model = original_model,
    original_coefs = original_coefs,
    all_sim_coefs = all_sim_coefs
  ))
}


# 5.  Run the sensitivity analysis (example)

results_vs_random <- sensitivity_analysis(data, n_simulations = 100, comparison_group = "random")
results_vs_regular <- sensitivity_analysis(data, n_simulations = 100, comparison_group = "regular")


# 6.  Examine the results
# 6.1  Compare original coefficients
print("Original Coefficients (vs random):")
print(results_vs_random$original_coefs)
print("Original Coefficients (vs regular):")
print(results_vs_regular$original_coefs)

# 6.2  Summarize simulation coefficients
print("\nSummary of Simulated Coefficients (vs random):")
print(results_vs_random$all_sim_coefs %>%
        group_by(term) %>%
        summarise(
          mean_estimate = mean(Estimate),
          mean_SE = mean(`Std. Error`),
          mean_pvalue   = mean(`Pr(>|z|)`)
        ))

print("\nSummary of Simulated Coefficients (vs regular):")
print(results_vs_regular$all_sim_coefs %>%
        group_by(term) %>%
        summarise(
          mean_estimate = mean(Estimate),
          mean_SE = mean(`Std. Error`),
          mean_pvalue   = mean(`Pr(>|z|)`)
        ))
