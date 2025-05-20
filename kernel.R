# ------------------------------------------------------------------------------
# Script: KDE of Plastic Debris Lines and Extraction to Points
# Description: Performs KDE on polyline debris data and extracts values to points
# ------------------------------------------------------------------------------
# Load libraries
library(sf)
library(terra)
library(tmap)
library(spatstat)
library(SpatialKDE)
library(dplyr)
library(stringr)

# ------------------------------------------------------------------------------
# 1. KDE on Debris Line Centroids
# ------------------------------------------------------------------------------

# Load debris lines already split to grid of 1m
debris <- st_read("path_to_shape.shp")
debris <- st_centroid(debris)
debris <- st_transform(debris, crs = 32627)
coords <- st_coordinates(debris)

# Create spatstat window
bbox <- st_bbox(debris)
win <- owin(xrange = c(bbox["xmin"], bbox["xmax"]),
            yrange = c(bbox["ymin"], bbox["ymax"]))
debris_ppp <- ppp(x = coords[, 1], y = coords[, 2], window = win)

# Estimate optimal bandwidth (optional)
bw <- bw.ppl(debris_ppp)  # Poisson likelihood cross-validation
band_width <- 3  # manually selected based on bw.ppl result

# KDE with different kernels
kde <- kde(debris, band_width = band_width, kernel = "epanechnikov", cell_size = 1)

# Save KDE shapefiles
st_write(kde, "kernel.shp", delete_dsn = TRUE)

# ------------------------------------------------------------------------------
# 2. Convert KDE Shapefiles to Rasters (1m Resolution)
# ------------------------------------------------------------------------------

# Create output folder
kde_folder <- "outputs/kernel_kde"
raster_output_folder <- file.path(kde_folder, "rasters")
dir.create(raster_output_folder, showWarnings = FALSE)

# List KDE shapefiles
shp_files <- list.files(kde_folder, pattern = "\\.shp$", full.names = TRUE)

# Loop through and rasterize
for (shp_file in shp_files) {
  kde_sf <- st_read(shp_file, quiet = TRUE)
  kde_vect <- vect(kde_sf)
  
  kde_col <- names(kde_sf)[sapply(kde_sf, is.numeric)][1]  # first numeric column
  
  r_template <- rast(ext(kde_vect), resolution = 1, crs = crs(kde_vect))
  kde_raster <- rasterize(kde_vect, r_template, field = kde_col)
  
  base_name <- tools::file_path_sans_ext(basename(shp_file))
  out_raster_path <- file.path(raster_output_folder, paste0(base_name, ".tif"))
  writeRaster(kde_raster, out_raster_path, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Extract KDE Values to Point Locations
# ------------------------------------------------------------------------------

# Load nest/random/regular point data and project
points_csv <- "path_to_csv.csv"
points_wgs84 <- read.csv(points_csv, sep = ";")
points_sf <- st_as_sf(points_wgs84, coords = c("lon", "lat"), crs = 4326)
points_utm <- st_transform(points_sf, crs = 32627)
points_utm$KDE_value <- NA_real_

# List raster files
raster_files <- list.files(raster_output_folder, pattern = "\\.tif$", full.names = TRUE)

# Extract KDE values by matching code
unique_codes <- unique(points_utm$Code)

for (code in unique_codes) {
  matching_raster <- raster_files[str_detect(basename(raster_files), fixed(code))]
  
  if (length(matching_raster) == 1) {
    r <- rast(matching_raster)
    points_subset <- points_utm %>% filter(Code == code)
    values <- terra::extract(r, vect(points_subset))[, 2]
    points_utm$KDE_value[points_utm$Code == code] <- values
  } else {
    warning(paste("No unique raster found for code:", code))
  }
}

# Reproject back to WGS84 and extract coordinates
points_wgs84 <- st_transform(points_utm, crs = 4326)
coords <- st_coordinates(points_wgs84)
points_wgs84$lon <- coords[, 1]
points_wgs84$lat <- coords[, 2]

# Save final CSV
write.csv(points_wgs84, "points_with_kernel", row.names = FALSE)
