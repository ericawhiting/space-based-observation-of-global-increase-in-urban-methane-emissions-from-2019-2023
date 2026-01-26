# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(sf)
library(sfheaders)
library(dplyr)
library(terra)
library(tidyterra)
library(ggplot2)
library(cowplot)
library(doParallel)
library(foreach)

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))

path_to_data <- "../Data/" # change your path here


# -------------------------------
# Load GHSL population, edgar ch4 and edgar co, and GHSL urbanization
# -------------------------------
path_to_population_data <- "/" # change your path here
pop_grid <- rast(paste0(path_to_data, "GHSL/", "GHS_POP_E2020_GLOBE_R2023A_4326_30ss_V1_0.tif"))
pop_grid <- project(pop_grid, "epsg:4326", threads = TRUE)



# -------------------------------
# Loop through each city
# -------------------------------
# List of city codes
cities <- c(return_codes_of_cities())

# remove parallelization call and change %dopar% to %do% if needed
registerDoParallel(16)
ideal_box_df <- foreach(ii = seq_along(cities), .combine = rbind) %dopar% {
    # find city code and city center (based on peak on population density)
    city <- cities[ii]
    city_latlon <- find_city_center_from_pop_dens(city) # from urban_boundaries.R
    lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0))
    # create extent
    large_urban_domain <- terra::ext(as.vector(cbind(city_latlon[2] -  1 * lat_adj, city_latlon[2] + 1 * lat_adj,
                                                city_latlon[1] - 1, city_latlon[1] + 1)))
    city_grid <- terra::crop(pop_grid, large_urban_domain)

    # set range of boxsizes (degrees lat lon from city center)
    delta_lat <- seq(0.3, 0.9, 0.01)

    delta_lon <- delta_lat * lat_adj

    # range boxsizes to find population density
    population_df <- data.frame(boxsize = numeric(), population = numeric(), area = numeric(), density = numeric())
    for (i in seq_along(delta_lat)) {
        urban_domain <- terra::ext(as.vector(cbind(city_latlon[2] - delta_lon[i], city_latlon[2] + delta_lon[i],
                                                  city_latlon[1] - delta_lat[i], city_latlon[1] + delta_lat[i])))
        cropped_pop <- terra::crop(city_grid, urban_domain)
        area_value <- terra::expanse(cropped_pop, unit = "km")
        sum_value <- terra::global(cropped_pop, "sum")
        population_df <- population_df %>% add_row(boxsize = delta_lat[i] * 2, # box size is the delta_lat size
                        population  = sum_value$sum,
                        area = area_value$area,
                        density = sum_value$sum / area_value$area)
    }

    # -------------------------------
    # normalize population density for each city and find where slope is -0.45
    # -------------------------------
    norm_df <- data.frame(boxsize = population_df$boxsize, norm_density = population_df$density / max(population_df$density))
    rise <- diff(norm_df$norm_density)
    run <- diff(norm_df$boxsize)
    slope <- rise / run
    ideal_box <- population_df[which.min(abs(slope - -0.45)), "boxsize"] # *2 of the delta lon, whole length here

    # ideal box is delta degrees latitude
    ideal_box_longitude_adj <- ideal_box * 1 / cos(city_latlon[1] * (pi / 180.0))

    return(data.frame(code = city, delta_lat = ideal_box / 2, delta_lon = ideal_box_longitude_adj / 2, total_box_length = ideal_box)) # total_box_length is the latitude degrees
}


# -------------------------------
# save boxsizes as csv
# -------------------------------
write.table(ideal_box_df, file = paste0(path_to_data, "ideal_boxsize.csv"), col.names = TRUE, sep = ",", row.names = FALSE)
