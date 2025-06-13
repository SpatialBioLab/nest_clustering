#### -- AVERAGE ELEVATION PROFILES OF BEACHES

library(terra)
library(sf)
library(dplyr)

# --- Data ---
# perpendicular transects to the shoreline
transects <- st_read("path_to_transects_shapefile") %>% st_transform(utm_crs) %>% filter(code == "code")

# nests points locations
nests <- read.csv2("path_to_csv") %>%
  mutate(across(c(Code, Point), as.factor), across(c(lon, lat), as.numeric)) %>%
  filter(Point == "nest") %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>% 
  st_transform(utm_crs) %>% 
  filter(Code == "code")

# marine debris points locations
debris <- st_read("path_to_debris_shapefile") %>% st_transform(utm_crs)

# digital surface model
dsm <- rast("path_to_dsm") %>% project(utm_crs)

# --- Extract Elevation Profiles Along Transects ---
profiles <- list()
global_min <- Inf; global_max <- -Inf

for (i in seq_len(nrow(transects))) {
  t_vect <- terra::vect(transects[i, ])
  elevs <- terra::extract(dsm, t_vect, xy = TRUE, points = TRUE)
  if (nrow(elevs) > 0) {
    df <- dplyr::select(elevs, x, y, elevation = 2)
    start <- st_coordinates(transects[i, ])[1, 1:2]
    df$cross_shore_distance <- sqrt((df$x - start[1])^2 + (df$y - start[2])^2)
    profiles[[i]] <- df
    global_min <- min(global_min, df$cross_shore_distance, na.rm = TRUE)
    global_max <- max(global_max, df$cross_shore_distance, na.rm = TRUE)
  }
}

# --- Average elevation profile ---
profile_avg <- data.frame()
if (length(profiles) > 0) {
  combined <- bind_rows(profiles)
  breaks <- seq(floor(global_min), ceiling(global_max), by = 2)
  if (length(breaks) < 2) breaks <- c(global_min, global_max + 1)
  
  profile_avg <- combined %>%
    group_by(dist_bin = cut(cross_shore_distance, breaks = breaks, labels = FALSE)) %>%
    summarise(
      cross_shore_center = median(cross_shore_distance, na.rm = TRUE),
      mean_elevation = mean(elevation, na.rm = TRUE),
      sd_elevation = sd(elevation, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cross_shore_center)
}

# --- Project nest points to transects ---
nests_proj <- data.frame()
if (nrow(nests) > 0) {
  nearest_idx <- st_nearest_feature(nests, transects)
  nests_proj <- purrr::map2_dfr(seq_len(nrow(nests)), nearest_idx, function(j, idx) {
    line <- transects[idx, ]
    proj <- st_nearest_points(nests[j, ], line) %>% st_cast("POINT")
    dist <- sqrt(sum((st_coordinates(proj)[1, ] - st_coordinates(line)[1, ])^2))
    elev <- terra::extract(dsm, vect(nests[j, ]))[, 2]
    data.frame(cross_shore_distance = dist, elevation = elev)
  })
}

# --- Sample and project plastic points ---

if (nrow(debris) > 0) {
  nearest_idx <- st_nearest_feature(debris, transects)
  plastic_proj <- purrr::map2_dfr(seq_len(nrow(debris)), nearest_idx, function(k, idx) {
    line <- transects[idx, ]
    proj <- st_nearest_points(debris[k, ], line) %>% st_cast("POINT")
    dist <- sqrt(sum((st_coordinates(proj)[1, ] - st_coordinates(line)[1, ])^2))
    elev <- terra::extract(dsm, vect(debris[k, ]))[, 2]
    data.frame(cross_shore_distance = dist, elevation = elev)
  })
}
