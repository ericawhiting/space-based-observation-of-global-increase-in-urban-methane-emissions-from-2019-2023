# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(dplyr)
library(terra)
library(tidyterra)

path_to_analysis_code <- "./" # change your path here

path_to_data <- "../Data/" # change your path here
city_information_path <- paste0(path_to_data, "city_information/city_information.csv")

# -------------------------------
load_city_information_df <- function() {
  df_city_information <- read.csv(city_information_path, header = TRUE)
  return(df_city_information)
}

read_city_info <- function(city) {
  #' find row of city information for city from city_information.csv
  #'
  #' @param city the 2 or 3 letter code unique to that city
  #' @return single row of df of city code, name, country, lat, lon, latlon

  df_city_information <- load_city_information_df()
  city_information <- df_city_information[df_city_information$code == city, ]
  if (dim(city_information)[1] != 1) {
    stop("non-unique city code in city_information.csv")
  }
  return(city_information)
}

read_ideal_boxsize_lat <- function(city) {
  # returns delta_lat
  df_box <- read.csv(paste0(path_to_data, "city_information/ideal_boxsize.csv"), header = TRUE)
  box_info <- df_box[df_box$code == city, ]
  return(box_info$delta_lat)
}

read_ideal_boxsize_lon <- function(city) {
  # returns delta_lon
  df_box <- read.csv(paste0(path_to_data, "city_information/ideal_boxsize.csv"), header = TRUE)
  box_info <- df_box[df_box$code == city, ]
  return(box_info$delta_lon)
}

find_city_latlon <- function(city) {
  #' find lat lon of city based on city code
  #'
  #' @param city the 2 or 3 letter code unique to that city
  #' @return the latitude and longitude in a vector

  city_latlon <- as.numeric(strsplit(read_city_info(city)$latlon, ",")[[1]])
  return(city_latlon)
}

find_city_name <- function(city) {
  #' find name (city, country) based on city code
  #'
  #' @param city the 2 or 3 letter code unique to that city
  #' @return the city name as character

  city_name <- read_city_info(city)$name
  return(city_name)
}
find_city <- function(city) {
  city_name <- read_city_info(city)$city
  return(city_name)
}
find_city_country <- function(city) {
  #' find country based on city code
  #'
  #' @param city the 2 or 3 letter code unique to that city
  #' @return the country name as character

  country_name <- read_city_info(city)$country
  return(country_name)
}

return_codes_of_cities <- function() {
  #' return list of city codes

  df_city_information <- load_city_information_df()
  return(df_city_information$code)
}

return_names_of_cities <- function() {
  #' return list of cities' names, country, must call with ()

  df_city_information <- load_city_information_df()
  return(df_city_information$name)
}

return_cities_in_country <- function(country) {
  #' return df of city info for cities in country
  #'
  #' @param country country name as string, ex: 'Japan' or 'USA'
  #' @return dataframe of city_information for cities in that country

  df_city_information <- load_city_information_df()
  city_information_by_country <- df_city_information[df_city_information$country == country, ]
  return(city_information_by_country)
}
return_cities_in_region <- function(region) {
  df_city_information <- load_city_information_df()
  city_information_by_region <- df_city_information[df_city_information$region == region, ]
  return(city_information_by_region)
}

is_C40 <- function(city) {
  city_C40 <- read_city_info(city)$C40 # 0 if not C40, 1 if C40
  # dmv is not C40 but dc is
  return(city_C40)
}

# return a specific column by name
get_column_by_name <- function(dataframe, column_name) {
  if (column_name %in% colnames(dataframe)) {
    return(dataframe[[column_name]])
  } else {
    stop(paste("Column name", column_name, "not found in the dataframe."))
  }
}

return_cities_in_region_var <- function(region, var) {
  # var is string
  df_city_information_by_region <- return_cities_in_region(region)
  df_var <- get_column_by_name(df_city_information_by_region, var)
  return(df_var)
}

return_list_of_regions <- function() {
  df_city_information <- load_city_information_df()
  return(unique(df_city_information$region))
}

find_city_region <- function(city) {
  city_region <- read_city_info(city)$region
  return(city_region)
}

find_city_coastal_flag <- function(city) {
  coastal_flag <- read_city_info(city)$coastal
  return(coastal_flag)
}

find_city_center_from_pop_dens <- function(city) {
  ### could instead read from city_center_from_pop_density.csv instead of needing to calculate by using find_city_latlon_adjusted below
  # output lat lon
  # original city center above has been used to crop large section of pop data, then we find the max in it
  city_latlon <- find_city_latlon(city) #returns lat, lon  (y,x) of city center
  lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0))
  latlon_diff <- c(-0.2, 0.2, -0.2 * lat_adj, 0.2 * lat_adj) # delta lat, delta lat, delta lon, delta lon
  #                y - delta y (lat),               y + delta y (lat),              x - delta x (lon),              x + delta x (lon)
  latlon_range <- c(city_latlon[1] + latlon_diff[1], city_latlon[1] + latlon_diff[2], city_latlon[2] + latlon_diff[3], city_latlon[2] + latlon_diff[4])

  df_pop_density <- read.csv(paste0(path_to_data, "city_information/pop_density_cropped/", city, "_cropped_pop_density.csv"), header = TRUE)
  df_pop_density <- df_pop_density[(df_pop_density["x"] > latlon_range[3]) &
                                   (df_pop_density["x"] < latlon_range[4]) &
                                   (df_pop_density["y"] > latlon_range[1]) &
                                   (df_pop_density["y"] < latlon_range[2]), ]
  peak_density_row <- df_pop_density[which.max(df_pop_density$value), ]

  center_lonlat <- c(peak_density_row$x[1],  peak_density_row$y[1]) # x,y (lon, lat)
  center_latlon <- rev(center_lonlat)
  return(center_latlon)
}

find_city_latlon_adjusted <- function(city) {
  # pull city center that was based oirignally on pop dens
  df_city_locations_adjusted <- read.csv(paste0(path_to_data, "city_information/city_center_from_pop_density.csv"), header = TRUE)
  df_city <- df_city_locations_adjusted[df_city_locations_adjusted$code == city, ]
  city_latlon <- c(df_city$lat, df_city$lon)
  return(city_latlon)
}

write_df_of_city_centers <- function() {
  cities <- c(return_codes_of_cities())
  city_center_pop_density <- data.frame(city = cities, lat = NA, lon = NA)
  for (i in seq_along(cities)) {
    city <- cities[i]
    latlon <- find_city_center_from_pop_dens(city)
    city_center_pop_density[i, "lat"] <- latlon[1]
    city_center_pop_density[i, "lon"] <- latlon[2]
  }
  write.table(city_center_pop_density, file = paste0(path_to_data, "city_information/city_center_from_pop_density.csv"), col.names = TRUE, sep = ",", row.names = FALSE)
}

find_EDGAR_uncertainty <- function(city) {
  #' pulls uncertainty on EDGAR emissions inventory values based on country or country group
  #' These values come from the SI of M. Crippa et al., Gridded emissions of air pollutants for the period 1970–2012 within Edgar v4.3.2. Earth Syst. Sci. Data 10, 1987–2013 (2018)

  country <- find_city_country(city)
  EU28_list <- c("EU28", "Belgium", "Bulgaria", "Czech Republic", "Denmark", "Germany", "Estonia", "Ireland", "Greece",
                 "Spain", "France", "Croatia", "Italy", "Cyprus", "Latvia", "Lithuania", "Luxembourg", "Hungary", "Malta",
                 "Netherlands", "The Netherlands", "Austria", "Poland", "Portugal", "Romania", "Slovenia", "Slovakia",
                 "Finland", "Sweden", "United Kingdom")
  if (country == "China") {
    uncert <- 94.4
  } else if (country == "India") {
    uncert <- 118.7
  } else if (country == "Brazil") {
    uncert <- 123.3
  } else if (country == "USA") {
    uncert <- 44.2
  } else if (country == "Russia") {
    uncert <- 25.9
  } else if (country %in% EU28_list) {
    uncert <- 64.6
  } else {
    uncert <- 108.1
  }
return(uncert / 100)
}

find_city_hemisphere <- function(city) {
  city_latlon <- find_city_latlon(city)
  if (city_latlon[1] > 0) {
    return("N")
  } else {
     return("S")
  }
}

return_city_boundary_flags <- function(city) {
  # due to neighboring cities, we use this flag to download TROPOMI data under one city code and access it
  # be careful to then remove these cities from analysis so as to not double count emissions!
  if ((city == "dc") || (city == "blt")) {
    return("dmv")
  } else if ((city == "naj") || (city == "zhe")) {
    return("nzh")
  } else if (city == "tok") {
    return("yok")
  } else if (city == "hk") {
    return("shz")
  } else if (city == "rot") {
    return("ams")
  } else if ((city == "eku") || (city == "tsh")) {
    return("joh")
  } else {
    return(city)
  }
}

return_urban_domain_rast <- function(city) {
  city_latlon <- find_city_center_from_pop_dens(city)
  delta_lat <- read_ideal_boxsize_lat(city)
  lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0)) # used to adjust delta longitudes for boxes based on where city is latitude wise
  delta_lon <- delta_lat * lat_adj
  coastal_flag <- find_city_coastal_flag(city) # return 1 if near coast or considered coastal city by eye
  delta_lat <- delta_lat + 0.25 * delta_lat * coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land
  delta_lon <- delta_lon + 0.25 * delta_lon * coastal_flag
  urban_domain_rast <- terra::ext(as.vector(cbind(city_latlon[2] - delta_lon, city_latlon[2] + delta_lon, city_latlon[1] - delta_lat, city_latlon[1] + delta_lat)))
  return(urban_domain_rast)
}