
### TESTING MORPHOLOGY EFFECT 
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


# REGRESSION BY POINT (width per surveyID)
m_litter <- lm(KDE_within_z ~ width_z + I(dist_z) + lon + lat, data = data)
data$litter_resid <- resid(m_litter)


# GAMM with litter_resid

# Nest vs. Random 
gam_random_morph <- gam(
  factor(Point == "nest") ~ s(dist_z, by = litter_resid) +
  s(surveyID, bs = "re"),
  family = binomial,
  data = data %>% filter(Point != "random"),
  method = "REML")
summary(gam_random_morph)

summary(gam_random_morph)
plot.gam(gam_random_morph, pages=1, shade=TRUE)
vis.gam(gam_random_morph, type="response")


# Nest vs. Regular 
gam_regular_morph <- gam(
  factor(Point == "nest") ~ s(dist_z, by = litter_resid) +
  s(surveyID, bs = "re"),
  family = binomial,
  data = data %>% filter(Point != "regular"),
  method = "REML")
summary(gam_regular_morph)

summary(gam_regular_morph)
plot.gam(gam_regular_morph, pages=1, shade=TRUE)
vis.gam(gam_regular_morph, type="response")
