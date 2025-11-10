# SENSITIVITY ANALYSIS

# Load necessary libraries
library(dplyr)      
library(MASS)      
library(tidyverse)
library(mgcv)

# data loading
data <- read.csv2("path_to_csv.csv")

# log-transformation of kernel # add a small offset to kernel to avoid log(0)
min_positive <- min(data$KDE[data$KDE > 0], na.rm = TRUE)
pseudo <- min_positive / 2   # a small offset (half the smallest positive)
data$log_KDE <- log(data$KDE + pseudo)

# between-within kernel and scaling
data <- data %>%
  group_by(surveyID) %>%
  mutate(KDE_mean = mean(log_KDE, na.rm = TRUE),
         KDE_within = log_KDE - KDE_mean) %>%
  ungroup() %>%
  mutate(dist_z = as.numeric(scale(as.numeric(dist_water))),
         KDE_mean_z = as.numeric(scale(KDE_mean)),
         KDE_within_z = as.numeric(scale(KDE_within)),
         dist_veget_z = as.numeric(scale(as.numeric(dist_veget))),
         width_z = as.numeric(scale(as.numeric(width))),
         lat = as.numeric(lat),
         lon = as.numeric(lon))


# function to add random error to the nest coordinates
add_error <- function(data, error_sd = 3) {
  
  # Create new columns with simulated data
  dist_water_sim <- as.numeric(data$dist_z) + rnorm(nrow(data), mean = 0, sd = error_sd)
  dist_veg_sim <- as.numeric(data$dist_veget_z) + rnorm(nrow(data), mean = 0, sd = error_sd)
  
  # Ensure distances don't go negative
  dist_water_sim <- pmax(dist_water_sim, 0)
  dist_veg_sim <- pmax(dist_veg_sim, 0)
  
  # Add the new columns to the data frame
  data$dist_water_sim <- dist_water_sim
  data$dist_veg_sim <- dist_veg_sim
  return(data)
}


# fit the model
fit_gam_model <- function(data, comparison_group = "random") {
  
  # Create binary variable
  is_nest <- ifelse(data$Point == "nest", 1, 0)
  
  # Add the binary variable to the data frame
  data$is_nest <- is_nest
  
  model <- gam(
    is_nest ~ s(dist_water_sim, KDE_within_z) + 
    KDE_mean_z +
    dist_veget_sim + 
    s(surveyID, bs = "re") + (1),
    family = binomial,
    data = data %>% filter(Point %in% c("nest", comparison_group)),
    method = "REML")
  return(model)
}

# extract and compare coefficients
extract_coefficients <- function(model) {
  coef_summary <- summary(model)$coefficients
  coef_df <- as.data.frame(coef_summary)
  coef_df$term <- rownames(coef_df)
  return(coef_df)
}

# sensitivity analysis function

sensitivity_analysis <- function(original_data, n_simulations = 100, comparison_group = "random") {
  original_data_sim <- add_error(original_data)
  original_model <- fit_gam_model(original_data_sim, comparison_group = comparison_group)
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


# run the sensitivity analysis

results_vs_random <- sensitivity_analysis(data, n_simulations = 100, comparison_group = "random")
results_vs_regular <- sensitivity_analysis(data, n_simulations = 100, comparison_group = "regular")

# summarize simulation coefficients
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
