# crop_pop_density.R
# E.Whiting 2022-6-17

# crop population density down to use when finding city center in find_city_center_from_pop_dens in urban_boundaries.R
# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(dplyr)
library(terra)

path_to_data <- "../Data/" # change your path here

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))

# -------------------------------
gridded_data_ghsl <- rast(paste0(path_to_data, "GHSL/", "GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif"))
gridded_data_ghsl <- project(gridded_data_ghsl, "epsg:4326", threads = TRUE)

# List of city codes
cities <- c(return_codes_of_cities())
for (ii in seq_along(cities)) {
  # find city information: given lat, lon
  city <- cities[ii]
  city_latlon <- find_city_latlon(city)
  lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0))
  latlon_diff <- c(-1, 1, -1 * lat_adj, 1 * lat_adj) # delta lat, delta lat, delta lon, delta lon
  urban_domain <- terra::ext(as.vector(cbind(city_latlon[2] + latlon_diff[3], city_latlon[2] + latlon_diff[4],
                                             city_latlon[1] + latlon_diff[1], city_latlon[1] + latlon_diff[2])))
  cropped_ghsl <- terra::crop(gridded_data_ghsl, urban_domain)
  pop_density_df <- as.data.frame(cropped_ghsl, xy = TRUE)
  names(pop_density_df)[3] <- "value"

  write.csv(pop_density_df, paste0(path_to_data, "city_information/pop_density_cropped/", city, "_cropped_pop_density.csv"), row.names = FALSE)
}
