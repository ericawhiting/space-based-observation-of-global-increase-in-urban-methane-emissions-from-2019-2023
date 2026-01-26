# how to run:
# >> cut -d, -f1 /pathtoData/city_information/city_codes.csv | xargs -n1 Rscript TROPOMI_ch4co_city_analysis.R

# -------------------------------
# Command Line Argument: city
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
city <- args[1] # eg: city = 'tok'

# -------------------------------
# Inputs
# -------------------------------
subset_method <- "box"


save_data <- FALSE
pixel_coverage <- 0.2 # 20% area coverage

save_plots <- FALSE
save_spatial_plots <- FALSE

#path for output plots
path_to_plots <- "../plots/tropomi_overpasses/" # change your path here

# check to see if temp_plots is empty before writing pdfs there
if (save_plots) {
  check_files <- list.files(paste0(path_to_plots, "temp_plots/"))
  if (length(check_files) > 0) {
    stop("temp_plots directory is not empty!")
  }
}

# -------------------------------
# for sensitivity studies
# -------------------------------
# Move City Center
center_direction <- "S"
center_move_x <- 0.00 # move longitude, negative value = W direction
center_move_y <- 0.00 # move lat, negative value = S direction
center_move <- center_move_x + center_move_y

# Change City Size
adjustment <- 1 # default is 1
# run between 0.55, .7, 0.85, 1, 1.15, 1.3, 1.45


# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
# removing some that are not used
library(cowsay)
library(ncdf4) # package for netcdf manipulation
library(ggplot2) # package for plotting
library(ggmap)
library(stringr)
library(pracma) # ffor string comparison function, ODR fit
library(cowplot) #for ggplot grid function
library(lmodel2) # for RMA linear fit
library(deming)# for deming fit
library(dplyr)
library(terra)
library(cowplot) #for plot_grid
library(sf)
library(sfheaders)
library(foreach)
library(doParallel)
library(qpdf) # to combine temp pdfs
library(CFtime)


path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))
source(paste0(path_to_utils, "calc_dry_column_mixing_ratio.R"))
source(paste0(path_to_utils, "orbits_to_filter.R"))
source(paste0(path_to_utils, "split_dates.R"))
source(paste0(path_to_utils, "plot_overpass.R"))

path_to_data <- "../Data/" # change your path here

# path to TROPOMI data
tropomi_path <- paste0(path_to_data, "TROPOMI_city_subsets/") # then we stored by year as subdirectories

# city_list comes from utils/write_city_list.R
city_list_path <- paste0(path_to_data, "/city_information/city_list_of_tropomi_files/")
city_list_fname <- paste0(city_list_path, city, "_files.csv")



# -------------------------------
# Set up Parallelization
# -------------------------------
num_cores <-  35
registerDoParallel(cores = num_cores)

# -------------------------------
# Urban Mask
# -------------------------------
# City Information
# Defining city center using pop density peak
city_latlon <- find_city_center_from_pop_dens(city)
city_name <- find_city_name(city)

# Read code of data for given city code, relevant if part of larger region (ex: dc as dmv)
domain_code <- return_city_boundary_flags(city)

# get filelist with paired TROPOMI CH4 and CO filesnames
city_files <- read.csv(city_list_fname) # ch4 and co columns of filenames with paths

# Pull Orbit numbers for co and ch4 files
co_orbit_list <- str_match(city_files$co, "_(\\d{5})_")[, 2]
ch4_orbit_list <- str_match(city_files$ch4, "_(\\d{5})_")[, 2]
city_orbits <- cbind(ch4_orbit_list, co_orbit_list)

# flag:
ch4_orbit_flag <- which(!ch4_orbit_list %in% co_orbit_list)
co_orbit_flag <- which(!co_orbit_list %in% ch4_orbit_list)
if ((length(ch4_orbit_flag) > 0) || (length(co_orbit_flag) > 0)) {
    paste0("non matching orbit detected")
    paste0("CH4 non matching orbit: ", ch4_orbit_flag)
    paste0("CO non matching orbit: ", co_orbit_flag)
}

# -------------------------------
# Assign Filters
# -------------------------------
land_classification_filter <- TRUE # both CO and CH4
co_scatter_filter <- TRUE
qa_filter <- TRUE #both CO and CH4
viirs_filter <- FALSE
viirs_orbit_filter <- TRUE # orbits where VIIRS was down

#pick ONLY one
ch4_hu_amt_2016_filter <- FALSE # filters from H. Hu AMT 2016
ch4_exp_filter <- TRUE # experimental filters from SRON

#CO filter thresholds
co_qa_filter <- 0 # qa > 0
prec_filt_co <- 5 # precision filter
t_threshold  <- 0.5  # clear sky optical thickness threshold
cloud_height_threshold <- 5000 # clear sky cloud height threshold

#CH4 filters thresholds
ch4_qa_filter <- 0 # qa >0
prec_filt_ch4 <- 10 # precision filter
ch4_aero <- 0.07 #SWIR aerosol optical thickness
ch4_albedo_limit <- 0.02

# orbit list to filter out
orbits_2019 <- return_orbits_to_filter(2019)
orbits_2020 <- return_orbits_to_filter(2020)
orbits_2021 <- return_orbits_to_filter(2021)
orbits_2022 <- return_orbits_to_filter(2022)
orbits_2023 <- return_orbits_to_filter(2023)

#Define box sizes
city_latlon <- find_city_center_from_pop_dens(city)
city_latlon <- city_latlon + c(center_move_y, center_move_x)
delta_lat <- read_ideal_boxsize_lat(city) * adjustment
lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0)) # used to adjust delta longitudes for boxes based on where city is latitude wise
delta_lon <- delta_lat * lat_adj

coastal_flag <- find_city_coastal_flag(city) # return 1 if near coast or considered coastal city by eye
delta_lat <- delta_lat + 0.25 * delta_lat * coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land
delta_lon <- delta_lon + 0.25 * delta_lon * coastal_flag


interpret_land_water <- function(value) {
  if (is.na(value)) {
    return(0)
  } else if (value %% 2 == 0) {
    return(1)
  } else {
    return(0)
  }
}

check_bits_15_or_19 <- function(flag) {
  bit_15_mask <- 0x00008000 # cloud_warning
  bit_19_mask <- 0x00080000 # low_cloud_fraction_warning
  bit_15_set <- (bitwAnd(flag, bit_15_mask) != 0)
  bit_19_set <- (bitwAnd(flag, bit_19_mask) != 0) # returns where bits are set
  return(bit_15_set || bit_19_set)
}


# -------------------------------
# Open Data from City List (list of tropomi files to pull), Initiate Vars
# -------------------------------

# -------------------------------
# Iterate through Dates/overpasses
# -------------------------------
#for (i in seq_len(dim(city_orbits)[1])) {
ch4_co_df <- foreach(i = seq_len(dim(city_orbits)[1]), .combine = rbind) %dopar% {
  # download data from data_source tropomi, does not work for tropomi_wfmdfind corresponding tropomi file
  # CH4
  nc_pathname_ch4 <- city_files[i, 1] # has path in it
  nc_filename_ch4 <- basename(city_files[i, 1])
  nc_ch4 <- nc_open(nc_pathname_ch4)
  date_ch4 <- substr(nc_filename_ch4, 21, 28)
  start_time_ch4 <- substr(nc_filename_ch4, 30, 35)
  orbit_ch4 <- substr(nc_filename_ch4, 53, 57) # another way to double check matching orbit numbers
  #CO
  nc_pathname_co <- city_files[i, 2] # has path in it
  nc_filename_co <- basename(city_files[i, 1])
  nc_co <- nc_open(nc_pathname_co)
  date_co <- substr(nc_filename_co, 21, 28)
  start_time_co <- substr(nc_filename_co, 30, 35)
  orbit_co <- substr(nc_filename_co, 53, 57)
  if (!strcmp(date_ch4, date_co)) {
    print("CH4 & CO date not for the same date")
    break
  }
  if (!strcmp(orbit_ch4, orbit_co)) {
    print("CH4 & CO orbit numbers do not match")
    break
  }


  # save date information, ch4 and co file dates verified to match above
  dates_vec <- date_ch4[[1]]
  # -------------------------------
  #read out tropomi data
  # -------------------------------
  #CH4

  #LOCATION DATA
  ch4_time_utc <- ncvar_get(nc_ch4, "PRODUCT/time_utc")
  ch4_lat <- ncvar_get(nc_ch4, "PRODUCT/latitude")
  ch4_lon <- ncvar_get(nc_ch4, "PRODUCT/longitude")
  ch4_lat_bounds <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")
  ch4_lon_bounds <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")

  #GAS DATA
  ch4_mixing_ratio_bc <- ncvar_get(nc_ch4, "PRODUCT/methane_mixing_ratio_bias_corrected")
  ch4_mixing_ratio_prec <- ncvar_get(nc_ch4, "PRODUCT/methane_mixing_ratio_precision")

  #Other data
  ch4_qa <- ncvar_get(nc_ch4, "PRODUCT/qa_value")
  ch4_quality_flag <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/processing_quality_flags")
  ch4_avg_kernel <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel", start = c(12, 1, 1, 1), count = c(1, dim(ch4_lat)[1], dim(ch4_lat)[2], 1))
  ch4_surface_pres <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure")
  ch4_opt_thickness <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_optical_thickness_SWIR")
  ch4_aero_size <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_size")
  ch4_aero_height <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/aerosol_mid_altitude")
  ch4_albedo <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_SWIR") # ch4_swir_albedo, forcing consistent naming across products
  ch4_nir_albedo <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_NIR")
  ch4_chisquare <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/chi_square")
  ch4_prior <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/methane_profile_apriori", start = c(12, 1, 1, 1), count = c(1, dim(ch4_lat)[1], dim(ch4_lat)[2], 1)) # to get prior column enhancement
  ch4_surface_class <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification")
  ch4_surface_alt <- ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude")
  eastward_wind <-  ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/eastward_wind")
  northward_wind <-  ncvar_get(nc_ch4, "PRODUCT/SUPPORT_DATA/INPUT_DATA/northward_wind")
  nc_close(nc_ch4)

  #CO

  #LOCATION DATA
  co_time_utc <- ncvar_get(nc_co, "PRODUCT/time_utc")
  co_lat <- ncvar_get(nc_co, "PRODUCT/latitude")
  co_lon <- ncvar_get(nc_co, "PRODUCT/longitude")
  co_lat_bounds <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds")
  co_lon_bounds <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds")

  #GAS DATA
  co_total_column <- ncvar_get(nc_co, "PRODUCT/carbonmonoxide_total_column") #carbonmonoxide_total_column_corrected
  co_total_column_prec <- ncvar_get(nc_co, "PRODUCT/carbonmonoxide_total_column_precision")
  h2o_total_column <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/water_total_column")

  #Other DATA
  co_qa <- ncvar_get(nc_co, "PRODUCT/qa_value")
  # pressure_levels = ncvar_get(nc_co,'DETAILED_RESULTS/pressure_levels')
  # solar_zenith_ang = ncvar_get(nc_co,'GEOLOCATIONS/solar_zenith_angle')
  # viewing_zenith_ang = ncvar_get(nc_co,'GEOLOCATIONS/viewing_zenith_angle')
  co_quality_flag <- ncvar_get(nc_co,  "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/processing_quality_flags")
  co_scat_opt_thickness <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/scattering_optical_thickness_SWIR")
  co_scat_height <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/height_scattering_layer")
  co_surface_alt <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude")
  co_surface_class <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_classification")
  co_albedo <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_2325") # co_albedo1, forcing consistent naming across products
  co_albedo2 <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/surface_albedo_2335")
  co_avg_kernel <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/DETAILED_RESULTS/column_averaging_kernel", start = c(50, 1, 1, 1), count = c(1, dim(co_lat)[1], dim(co_lat)[2], 1))
  co_surface_pres <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_pressure")
  co_prior <- ncvar_get(nc_co, "PRODUCT/SUPPORT_DATA/INPUT_DATA/carbonmonoxide_profile_apriori", start = c(50, 1, 1, 1), count = c(1, dim(co_lat)[1], dim(co_lat)[2], 1)) # prior column enhancement
  nc_close(nc_co)

  # might want to add back in original cropping so that spatial plots show more than just urban domain
  # -------------------------------
  # Crop to city
  # -------------------------------
  city_coords <- cbind(city_latlon[1] - delta_lat, city_latlon[1] + delta_lat,
                         city_latlon[2] - delta_lon, city_latlon[2] + delta_lon)
  city_mask <- which((as.vector(ch4_lat) > city_coords[1]) & (as.vector(ch4_lat) < city_coords[2]) &
                      (as.vector(ch4_lon) > city_coords[3]) & (as.vector(ch4_lon) < city_coords[4]))
  #ch4: subset data within box boundary
  ch4_lat_city <- as.vector(ch4_lat)[city_mask]
  ch4_lon_city <- as.vector(ch4_lon)[city_mask]
  ch4_lat_b1 <- as.vector(ch4_lat_bounds[1, , ])[city_mask]
  ch4_lat_b2 <- as.vector(ch4_lat_bounds[2, , ])[city_mask]
  ch4_lat_b3 <- as.vector(ch4_lat_bounds[3, , ])[city_mask]
  ch4_lat_b4 <- as.vector(ch4_lat_bounds[4, , ])[city_mask]
  ch4_lon_b1 <- as.vector(ch4_lon_bounds[1, , ])[city_mask]
  ch4_lon_b2 <- as.vector(ch4_lon_bounds[2, , ])[city_mask]
  ch4_lon_b3 <- as.vector(ch4_lon_bounds[3, , ])[city_mask]
  ch4_lon_b4 <- as.vector(ch4_lon_bounds[4, , ])[city_mask]
  ch4_mr_city <- as.vector(ch4_mixing_ratio_bc)[city_mask]
  ch4_avg_kernel_city <- as.vector(ch4_avg_kernel)[city_mask]
  ch4_prec_city <- as.vector(ch4_mixing_ratio_prec)[city_mask]
  ch4_albedo_city <- as.vector(ch4_albedo)[city_mask]
  ch4_nir_albedo_city <- as.vector(ch4_nir_albedo)[city_mask]
  ch4_prior_city <- as.vector(ch4_prior)[city_mask]
  ch4_qa_city <- as.vector(ch4_qa)[city_mask]
  ch4_quality_flag_city <- as.vector(ch4_quality_flag)[city_mask]
  ch4_surface_pres_city <- as.vector(ch4_surface_pres)[city_mask]
  ch4_opt_thickness_city <- as.vector(ch4_opt_thickness)[city_mask]
  ch4_aero_size_city <- as.vector(ch4_aero_size)[city_mask]
  ch4_aero_height_city <- as.vector(ch4_aero_height)[city_mask]
  ch4_chisquare_city <- as.vector(ch4_chisquare)[city_mask]
  ch4_surface_alt_city <- as.vector(ch4_surface_alt)[city_mask]
  ch4_surface_class_city <- as.vector(ch4_surface_class)[city_mask]
  eastward_wind_city <- as.vector(eastward_wind)[city_mask]
  northward_wind_city <- as.vector(northward_wind)[city_mask]

  #co: subset data within box boundary
  co_lat_city <- as.vector(co_lat)[city_mask]
  co_lon_city <- as.vector(co_lon)[city_mask]
  co_lat_b1 <- as.vector(co_lat_bounds[1, , ])[city_mask]
  co_lat_b2 <- as.vector(co_lat_bounds[2, , ])[city_mask]
  co_lat_b3 <- as.vector(co_lat_bounds[3, , ])[city_mask]
  co_lat_b4 <- as.vector(co_lat_bounds[4, , ])[city_mask]
  co_lon_b1 <- as.vector(co_lon_bounds[1, , ])[city_mask]
  co_lon_b2 <- as.vector(co_lon_bounds[2, , ])[city_mask]
  co_lon_b3 <- as.vector(co_lon_bounds[3, , ])[city_mask]
  co_lon_b4 <- as.vector(co_lon_bounds[4, , ])[city_mask]
  co_total_column_city <- as.vector(co_total_column)[city_mask]
  co_total_column_prec_city <- as.vector(co_total_column_prec)[city_mask]
  h2o_total_column_city <- as.vector(h2o_total_column)[city_mask]
  co_avg_kernel_city <- as.vector(co_avg_kernel)[city_mask]
  co_scat_opt_thickness_city <- as.vector(co_scat_opt_thickness)[city_mask]
  co_scat_height_city <- as.vector(co_scat_height)[city_mask]
  co_surface_pres_city <- as.vector(co_surface_pres)[city_mask]
  co_qa_city <- as.vector(co_qa)[city_mask]
  co_albedo_city <- as.vector(co_albedo)[city_mask]
  co_prior_city <- as.vector(co_prior)[city_mask]
  co_surface_alt_city <- as.vector(co_surface_alt)[city_mask]
  co_surface_class_city <- as.vector(co_surface_class)[city_mask]

  # -------------------------------
  # calculate dry column mixing ratio for CO data      # from calc_dry_column_mixing_ratio.R
  # -------------------------------
  dry_temp <- calc_dry_column_mixing_ratio(co_surface_pres_city, co_total_column_city, co_total_column_prec_city, h2o_total_column_city)
  co_mr_city <- dry_temp[[1]]
  co_prec_city <- dry_temp[[2]]

  # -------------------------------
  # Apply data filters ------------
  # -------------------------------
  if (co_scatter_filter) {
    #CO scaling filter and precision
    ind_filt <- which((co_scat_opt_thickness_city >= t_threshold) |
                      (co_scat_height_city >= cloud_height_threshold) |
                      co_prec_city > prec_filt_co |
                      co_albedo_city < 0.02)
    co_mr_city[ind_filt] <- NA
    co_prec_city[ind_filt] <- NA
  }

  if (ch4_hu_amt_2016_filter) {
    ind_filt2 <- which(ch4_opt_thickness_city * ch4_aero_height_city / ch4_aero_size_city > 120) # threshold from H. Hu 2016 AMT
    # left off here
    ch4_mr_city[ind_filt2] <- NA
    ch4_prec_city[ind_filt2] <- NA

    ind_filt3 <- which(ch4_albedo_city < 0.02) # threshold from H. Hu 2016 AMT
    ch4_mr_city[ind_filt3] <- NA
    ch4_prec_city[ind_filt3] <- NA
  }

  if (ch4_exp_filter) {
    #SWIR aerosol optical thickness and precision
    ind_filt6 <- which(ch4_opt_thickness_city >= ch4_aero | ch4_prec_city > prec_filt_ch4)
    ch4_mr_city[ind_filt6] <- NA
    ch4_prec_city[ind_filt6] <- NA
    #Note qa filter is below

    ind_filt3 <- which(ch4_albedo_city < ch4_albedo_limit) # threshold from H. Hu 2016 AMT
    ch4_mr_city[ind_filt3] <- NA
    ch4_prec_city[ind_filt3] <- NA
  }

  if (qa_filter) {
    #qa filters
    #CH4
    ind_filt4 <- which(ch4_qa_city <= ch4_qa_filter)
    ch4_mr_city[ind_filt4] <- NA
    ch4_prec_city[ind_filt4] <- NA
    #CO
    ind_filt5 <- which(co_qa_city <= co_qa_filter)
    co_mr_city[ind_filt5] <- NA
    co_prec_city[ind_filt5] <- NA
  }

  if (viirs_filter) {
    cloud_warning <- sapply(ch4_quality_flag_city, check_bits_15_or_19)
    ind_filt9 <- which(cloud_warning == TRUE)
    ch4_mr_city[ind_filt9] <- NA
    ch4_prec_city[ind_filt9] <- NA
  }

  if (viirs_orbit_filter) {
    if (orbit_co %in% c(orbits_2019, orbits_2020, orbits_2021, orbits_2022, orbits_2023)) {
      ch4_mr_city[] <- NA
      ch4_prec_city[] <- NA
    }
  }

  # land classification
  if (land_classification_filter) {
    # filter ch4 water pixels
    land_water_classification_ch4 <- sapply(ch4_surface_class_city, interpret_land_water)
    ind_filt7 <- which(land_water_classification_ch4 == 0)
    land_water_classification_ch4[ind_filt7] <- NA
    # filter co water pixels
    land_water_classification_co <- sapply(co_surface_class_city, interpret_land_water)
    ind_filt8 <- which(land_water_classification_co == 0)
    land_water_classification_co[ind_filt8] <- NA
  }
  # -------------------------------
  # save data
  # -------------------------------
  filtered_data <- data.frame(ch4 = ch4_mr_city, ch4_avg_kernel = ch4_avg_kernel_city, ch4_surf_prior = ch4_prior_city, ch4_prec = ch4_prec_city,
                              ch4_albedo = ch4_albedo_city, ch4_nir_albedo = ch4_nir_albedo_city,  ch4_qa = ch4_qa_city, ch4_quality_flag = ch4_quality_flag_city,
                              ch4_surface_pres = ch4_surface_pres_city, ch4_opt_thickness = ch4_opt_thickness_city, ch4_aero_size = ch4_aero_size_city,
                              ch4_aero_height = ch4_aero_height_city, ch4_chisquare = ch4_chisquare_city, ch4_surface_alt = ch4_surface_alt_city,
                              ch4_surface_class = ch4_surface_class_city,
                              #land_water_classification_ch4 = land_water_classification_ch4, land_water_classification_co = land_water_classification_co,
                              co = co_mr_city, co_avg_kernel = co_avg_kernel_city, co_surf_prior = co_prior_city,
                              co_prec = co_prec_city, co_scat_opt_thickness = co_scat_opt_thickness_city, co_scat_height = co_scat_height_city,
                              co_surface_pres = co_surface_pres_city, co_qa = co_qa_city, co_albedo = co_albedo_city, co_surface_alt = co_surface_alt_city,
                              co_surface_class = co_surface_class_city, ch4_lat = ch4_lat_city, ch4_lon = ch4_lon_city, co_lat = co_lat_city, co_lon = co_lon_city,
                              ch4_lat_b1 = ch4_lat_b1, ch4_lat_b2 = ch4_lat_b2, ch4_lat_b3 = ch4_lat_b3, ch4_lat_b4 = ch4_lat_b4,
                              ch4_lon_b1 = ch4_lon_b1, ch4_lon_b2 = ch4_lon_b2, ch4_lon_b3 = ch4_lon_b3, ch4_lon_b4 =  ch4_lon_b4,
                              co_lat_b1 = co_lat_b1, co_lat_b2 = co_lat_b2, co_lat_b3 = co_lat_b3, co_lat_b4 = co_lat_b4,
                              co_lon_b1 = co_lon_b1, co_lon_b2 = co_lon_b2, co_lon_b3 = co_lon_b3, co_lon_b4 =  co_lon_b4,
                              eastward_wind = eastward_wind_city, northward_wind = northward_wind_city)
  full_num_pixels <- dim(filtered_data)[1]
  filtered_data <- filtered_data[complete.cases(filtered_data), ]

  fraction_filtered_data <- dim(filtered_data)[1] / full_num_pixels
  if (is.na(fraction_filtered_data)) {
    fraction_filtered_data <- 0
  }

  # Analyze if there are enough co-located ch4, co pixels that pass filtering requirements (pixel threshold)
  if (fraction_filtered_data >= pixel_coverage) {
    ch4_spatial_df <- data.frame(ch4_lat = ch4_lat_city, ch4_lon = ch4_lon_city,
                                 ch4_lat_b1 = ch4_lat_b1, ch4_lat_b2 = ch4_lat_b2, ch4_lat_b3 = ch4_lat_b3, ch4_lat_b4 = ch4_lat_b4,
                                 ch4_lon_b1 = ch4_lon_b1, ch4_lon_b2 = ch4_lon_b2, ch4_lon_b3 = ch4_lon_b3, ch4_lon_b4 = ch4_lon_b4)
    co_spatial_df <- data.frame(co_lat = co_lat_city, co_lon = co_lon_city,
                                co_lat_b1 = co_lat_b1, co_lat_b2 = co_lat_b2, co_lat_b3 = co_lat_b3, co_lat_b4 = co_lat_b4,
                                co_lon_b1 = co_lon_b1, co_lon_b2 = co_lon_b2, co_lon_b3 = co_lon_b3, co_lon_b4 = co_lon_b4)
    # fill back in with  NAs to make spatial plotting easy
    filtered_filled_data_ch4 <- merge(filtered_data %>%
                                         select(ch4, ch4_lat, ch4_lon, ch4_lat_b1, ch4_lat_b2, ch4_lat_b3, ch4_lat_b4, ch4_lon_b1, ch4_lon_b2, ch4_lon_b3, ch4_lon_b4),
                                      ch4_spatial_df, by = c("ch4_lat", "ch4_lon", "ch4_lat_b1", "ch4_lat_b2", "ch4_lat_b3", "ch4_lat_b4", "ch4_lon_b1", "ch4_lon_b2", "ch4_lon_b3", "ch4_lon_b4"),
                                      all.y = TRUE)

    filtered_filled_data_co <- merge(filtered_data %>%
                                         select(co, co_lat, co_lon, co_lat_b1, co_lat_b2, co_lat_b3, co_lat_b4, co_lon_b1, co_lon_b2, co_lon_b3, co_lon_b4),
                                     co_spatial_df, by = c("co_lat", "co_lon", "co_lat_b1", "co_lat_b2", "co_lat_b3", "co_lat_b4", "co_lon_b1", "co_lon_b2", "co_lon_b3", "co_lon_b4"),
                                     all.y = TRUE)
    # -------------------------------
    # calculate background for mr using quantiles
    # -------------------------------
    # initialize background as min
    mr_bkg_percentile <- data.frame("CH4_mr_bkg" = min(filtered_data$ch4), "CO_mr_bkg" = min(filtered_data$co))

    for (percentile in seq(5, 25, by = 1)) {
      ch4_mr_bkg <- quantile(filtered_data$ch4, percentile / 100)
      co_mr_bkg <- quantile(filtered_data$co, percentile / 100)
      #combine
      mr_bkg_percentile <- rbind(mr_bkg_percentile, data.frame("CH4_mr_bkg" = ch4_mr_bkg, "CO_mr_bkg" = co_mr_bkg))
    }
    rownames(mr_bkg_percentile) <- c(c("min", "5%", "6%", "7%", "8%", "9%", "10%", "11%", "12%", "13%", "14%", "15%",
                                      "16%", "17%", "18%", "19%", "20%", "21%", "22%", "23%", "24%", "25%"))


    # save pixel mixing ratio backgrounds calculated with 15 percentile background
    ch4_mr_bkg_percentile <- mr_bkg_percentile[12, 1]
    co_mr_bkg_percentile <- mr_bkg_percentile[12, 2]

    # -------------------------------
    # calculate background or non-urban prior using percentile to find background pixels, listed for ch4, co
    # same method as finding background in mr but could be different pixels
    # listed as method 3 in compare_prior_enhancement_methods.R
    # -------------------------------
     # initialize background as min
    prior_bkg_percentile <- data.frame("CH4_prior_bkg" = min(filtered_data$ch4_surf_prior), "CO_prior_bkg" = min(filtered_data$co_surf_prior))

    for (percentile in seq(5, 25, by = 1)) {
      ch4_prior_bkg <- quantile(filtered_data$ch4_surf_prior, percentile / 100)
      co_prior_bkg <- quantile(filtered_data$co_surf_prior, percentile / 100)
      #combine
      prior_bkg_percentile <- rbind(prior_bkg_percentile, data.frame("CH4_prior_bkg" = ch4_prior_bkg, "CO_prior_bkg" = co_prior_bkg))
    }
    rownames(prior_bkg_percentile) <- c(c("min", "5%", "6%", "7%", "8%", "9%", "10%", "11%", "12%", "13%", "14%", "15%",
                                          "16%", "17%", "18%", "19%", "20%", "21%", "22%", "23%", "24%", "25%"))

    # save pixel prior backgrounds calculated with 15 percentile background
    ch4_prior_bkg_percentile <- prior_bkg_percentile[12, 1]
    co_prior_bkg_percentile <- prior_bkg_percentile[12, 2]
    ch4_prior_mean_enh <- mean(filtered_data$ch4_surf_prior - prior_bkg_percentile[12, 1])
    co_prior_mean_enh <- mean(filtered_data$co_surf_prior - prior_bkg_percentile[12, 2])
    ch4_prior_sd <- sd(filtered_data$ch4_surf_prior)
    ch4_prior_mean <- mean(filtered_data$ch4_surf_prior)
    co_prior_sd <- sd(filtered_data$co_surf_prior)
    co_prior_mean <- mean(filtered_data$co_surf_prior)
    prior_term_ch4 <- mean((1 - filtered_data$ch4_avg_kernel) * (filtered_data$ch4_surf_prior - prior_bkg_percentile[12, 1]) / filtered_data$ch4_avg_kernel)
    prior_term_co <- mean((1 - filtered_data$co_avg_kernel) * (filtered_data$co_surf_prior - prior_bkg_percentile[12, 2]) / filtered_data$co_avg_kernel)
    # -------------------------------
    # Calculate Ratio of Summed Enhancements:
    # subtract off background, scale by surface averaging kernel, adjust for prior bias to urban regions
    # -------------------------------
    # enhancement at each pixel for one overpass # dim length of filtered data, 5 (each type of background calc)

    ch4_enh <- (outer(filtered_data$ch4, mr_bkg_percentile[, 1], "-") / filtered_data$ch4_avg_kernel)
    ch4_prior_adjustment <- ((1 - filtered_data$ch4_avg_kernel) * outer(filtered_data$ch4_surf_prior, prior_bkg_percentile[, 1], "-") / filtered_data$ch4_avg_kernel)
    ch4_numerator <- ch4_enh - ch4_prior_adjustment

    co_enh <- (outer(filtered_data$co, mr_bkg_percentile[, 2], "-") / filtered_data$co_avg_kernel)
    co_prior_adjustment <- ((1 - filtered_data$co_avg_kernel)  * outer(filtered_data$co_surf_prior,  prior_bkg_percentile[, 2], "-") / filtered_data$co_avg_kernel)
    co_denominator <- co_enh - co_prior_adjustment

    prior_mean_impact_co <- mean(co_enh[, 12] / co_denominator[, 12]) # pulls 15th percentile
    prior_mean_impact_ch4 <- mean(ch4_enh[, 12] / ch4_numerator[, 12]) # does not automatically exclude NA

    filtered_data$ch4_numerator <- ch4_numerator
    filtered_data$co_denominator <- co_denominator
    # summed enhancements for ch4 and co: one summed ch4 enhancement and one summed co enhancement per orbit
    summed_ch4_enhancements <- colSums(ch4_numerator) # length of 22
    summed_co_enhancements <- colSums(co_denominator) # length of 22

    mean_ch4_enhancements <- colMeans(ch4_numerator)
    mean_co_enhancements <- colMeans(co_denominator)

    mean_eastward_wind <- mean(filtered_data$eastward_wind)
    mean_northward_wind <- mean(filtered_data$northward_wind)
    # -------------------------------
    # SNR Filtering with CO
    # -------------------------------
    n_points_vec <- dim(filtered_data)[1]

    if (summed_co_enhancements[12] > (n_points_vec * 7)) { # checks 15% background value
      # ratio of summed enhancements
      ch4co_ratio_of_summed_enhancements <- summed_ch4_enhancements / summed_co_enhancements # dim: city_orbits, background options
      # -------------------------------
      # output mean of precisions & averaging kernels, and correlation coeffs
      # -------------------------------
      mean_co_prec <- mean(filtered_data$co_prec, na.rm = TRUE) # length 1
      mean_ch4_prec <- mean(filtered_data$ch4_prec, na.rm = TRUE) # length 1
      prec_cov <- cov(filtered_data$co_prec, filtered_data$co_prec)
      mean_ak_ch4 <- mean(filtered_data$ch4_avg_kernel, na.rm = TRUE) # length 1
      mean_ak_co <- mean(filtered_data$co_avg_kernel, na.rm = TRUE) # length 1
      mean_ch4_albedo <- mean(filtered_data$ch4_albedo) # length 1
      mean_co_albedo <- mean(filtered_data$co_albedo) # length 1
      # correlation coef.
      corr_vec_mr <- cor(filtered_data$ch4, filtered_data$co) # correlation of mixing ratio across the overpass at pixel level, length 1
      corr_vec_enh <- diag(cor(ch4_numerator, co_denominator)) # correlation of the enhancements at pixel level, length 5

      corr_ch4_albedo <- cor(filtered_data$ch4_albedo, filtered_data$ch4, use = "complete.obs")
      corr_co_albedo <- cor(filtered_data$co_albedo, filtered_data$co, use = "complete.obs")

      corr_ch4_aot <- cor(filtered_data$ch4_opt_thickness, filtered_data$ch4, use = "complete.obs")
      corr_co_aot <- cor(filtered_data$co_scat_opt_thickness, filtered_data$co, use = "complete.obs")

      # -------------------------------
      # plot data if it passes p-value or SNR threshold
      # -------------------------------
      if (save_spatial_plots) {
          pdf(paste(path_to_plots, "temp_plots/", city, "_", i, ".pdf", sep = ""))
          # plot mixing ratios
          # only looks at urban domain, not outside/zoomed out at all
          ch4_df2 <- DataFrameforPlot2(filtered_filled_data_ch4$ch4, filtered_filled_data_ch4$ch4_lat, filtered_filled_data_ch4$ch4_lon,
                                       filtered_filled_data_ch4$ch4_lat_b1, filtered_filled_data_ch4$ch4_lat_b2, filtered_filled_data_ch4$ch4_lat_b3, filtered_filled_data_ch4$ch4_lat_b4,
                                       filtered_filled_data_ch4$ch4_lon_b1, filtered_filled_data_ch4$ch4_lon_b2, filtered_filled_data_ch4$ch4_lon_b3, filtered_filled_data_ch4$ch4_lon_b4)

          co_df2 <- DataFrameforPlot2(filtered_filled_data_co$co, filtered_filled_data_co$co_lat, filtered_filled_data_co$co_lon,
                                      filtered_filled_data_co$co_lat_b1, filtered_filled_data_co$co_lat_b2, filtered_filled_data_co$co_lat_b3, filtered_filled_data_co$co_lat_b4,
                                      filtered_filled_data_co$co_lon_b1, filtered_filled_data_co$co_lon_b2, filtered_filled_data_co$co_lon_b3, filtered_filled_data_co$co_lon_b4)

          q <- c(0.05, 0.99)
          p1 <- PlotRegion2(ch4_df2, city_coords, title = paste("CH4", date_ch4), units = "ppb", q = q, printPlot = FALSE, box.coords = city_coords)
          p2 <- PlotRegion2(co_df2, city_coords,  title = paste("CO", date_co), units = "ppb", q = q, printPlot = FALSE, box.coords = city_coords)
          p3 <- cowplot::plot_grid(p1, p2, nrow = 1)
          title_theme <- ggdraw() +
                        draw_label(paste0(city_name, "\nCH4:CO = ",  round(ch4co_ratio_of_summed_enhancements[12], 3),
                                          "\nn_points_vec = ", n_points_vec, "\nch4co mr correlation = ", round(corr_vec_mr, 3),
                                          "\nch4co ehn correlation = ", round(corr_vec_enh[12], 3)),
                                      x = 0.1, hjust = 0)
                                      p4 <- cowplot::plot_grid(title_theme, p3, ncol = 1, rel_heights = c(0.2, 1))
          print(p4)
          dev.off()
        }
      # fix use of date_ch4 and orbit_ch4
      return(data.frame(dates_vec = date_ch4[[1]], orbit = orbit_ch4[1],
                        summed_ch4_enhancements_min = summed_ch4_enhancements[[1]], summed_ch4_enhancements_05 = summed_ch4_enhancements[[2]],  summed_ch4_enhancements_06 = summed_ch4_enhancements[[3]],
                        summed_ch4_enhancements_07 = summed_ch4_enhancements[[4]],  summed_ch4_enhancements_08 = summed_ch4_enhancements[[5]],  summed_ch4_enhancements_09 = summed_ch4_enhancements[[6]],
                        summed_ch4_enhancements_10 = summed_ch4_enhancements[[7]],  summed_ch4_enhancements_11 = summed_ch4_enhancements[[8]],  summed_ch4_enhancements_12 = summed_ch4_enhancements[[9]],
                        summed_ch4_enhancements_13 = summed_ch4_enhancements[[10]], summed_ch4_enhancements_14 = summed_ch4_enhancements[[11]], summed_ch4_enhancements_15 = summed_ch4_enhancements[[12]],
                        summed_ch4_enhancements_16 = summed_ch4_enhancements[[13]], summed_ch4_enhancements_17 = summed_ch4_enhancements[[14]], summed_ch4_enhancements_18 = summed_ch4_enhancements[[15]],
                        summed_ch4_enhancements_19 = summed_ch4_enhancements[[16]], summed_ch4_enhancements_20 = summed_ch4_enhancements[[17]], summed_ch4_enhancements_21 = summed_ch4_enhancements[[18]],
                        summed_ch4_enhancements_22 = summed_ch4_enhancements[[19]], summed_ch4_enhancements_23 = summed_ch4_enhancements[[20]], summed_ch4_enhancements_24 = summed_ch4_enhancements[[21]],
                        summed_ch4_enhancements_25 = summed_ch4_enhancements[[22]],

                        summed_co_enhancements_min = summed_co_enhancements[[1]], summed_co_enhancements_05 = summed_co_enhancements[[2]],  summed_co_enhancements_06 = summed_co_enhancements[[3]],
                        summed_co_enhancements_07 = summed_co_enhancements[[4]],  summed_co_enhancements_08 = summed_co_enhancements[[5]],  summed_co_enhancements_09 = summed_co_enhancements[[6]],
                        summed_co_enhancements_10 = summed_co_enhancements[[7]],  summed_co_enhancements_11 = summed_co_enhancements[[8]],  summed_co_enhancements_12 = summed_co_enhancements[[9]],
                        summed_co_enhancements_13 = summed_co_enhancements[[10]], summed_co_enhancements_14 = summed_co_enhancements[[11]], summed_co_enhancements_15 = summed_co_enhancements[[12]],
                        summed_co_enhancements_16 = summed_co_enhancements[[13]], summed_co_enhancements_17 = summed_co_enhancements[[14]], summed_co_enhancements_18 = summed_co_enhancements[[15]],
                        summed_co_enhancements_19 = summed_co_enhancements[[16]], summed_co_enhancements_20 = summed_co_enhancements[[17]], summed_co_enhancements_21 = summed_co_enhancements[[18]],
                        summed_co_enhancements_22 = summed_co_enhancements[[19]], summed_co_enhancements_23 = summed_co_enhancements[[20]], summed_co_enhancements_24 = summed_co_enhancements[[21]],
                        summed_co_enhancements_25 = summed_co_enhancements[[22]],

                        integ_results_min = ch4co_ratio_of_summed_enhancements[[1]], integ_results_05 = ch4co_ratio_of_summed_enhancements[[2]],  integ_results_06 = ch4co_ratio_of_summed_enhancements[[3]],
                        integ_results_07 = ch4co_ratio_of_summed_enhancements[[4]],  integ_results_08 = ch4co_ratio_of_summed_enhancements[[5]],  integ_results_09 = ch4co_ratio_of_summed_enhancements[[6]],
                        integ_results_10 = ch4co_ratio_of_summed_enhancements[[7]],  integ_results_11 = ch4co_ratio_of_summed_enhancements[[8]],  integ_results_12 = ch4co_ratio_of_summed_enhancements[[9]],
                        integ_results_13 = ch4co_ratio_of_summed_enhancements[[10]], integ_results_14 = ch4co_ratio_of_summed_enhancements[[11]], integ_results_15 = ch4co_ratio_of_summed_enhancements[[12]],
                        integ_results_16 = ch4co_ratio_of_summed_enhancements[[13]], integ_results_17 = ch4co_ratio_of_summed_enhancements[[14]], integ_results_18 = ch4co_ratio_of_summed_enhancements[[15]],
                        integ_results_19 = ch4co_ratio_of_summed_enhancements[[16]], integ_results_20 = ch4co_ratio_of_summed_enhancements[[17]], integ_results_21 = ch4co_ratio_of_summed_enhancements[[18]],
                        integ_results_22 = ch4co_ratio_of_summed_enhancements[[19]], integ_results_23 = ch4co_ratio_of_summed_enhancements[[20]], integ_results_24 = ch4co_ratio_of_summed_enhancements[[21]],
                        integ_results_25 = ch4co_ratio_of_summed_enhancements[[22]],

                        mean_ch4_enhancement_15 = mean_ch4_enhancements[[12]], mean_co_enhancement_15 = mean_co_enhancements[[12]],
                        ch4_mr_bkg_percentile = ch4_mr_bkg_percentile, co_mr_bkg_percentile = co_mr_bkg_percentile, corr_vec_mr = corr_vec_mr, corr_vec_enh = corr_vec_enh[[12]],
                        ch4_prior_bkg_percentile = ch4_prior_bkg_percentile, co_prior_bkg_percentile = co_prior_bkg_percentile,
                        ch4_prior_mean_enh = ch4_prior_mean_enh, co_prior_mean_enh = co_prior_mean_enh, ch4_prior_mean = ch4_prior_mean, co_prior_mean = co_prior_mean, ch4_prior_sd = ch4_prior_sd, co_prior_sd = co_prior_sd,
                        mean_ch4_prec = mean_ch4_prec, mean_co_prec = mean_co_prec, prec_cov = prec_cov,
                        mean_ak_ch4 = mean_ak_ch4, mean_ak_co = mean_ak_co, mean_ch4_albedo = mean_ch4_albedo, mean_co_albedo = mean_co_albedo, corr_ch4_albedo = corr_ch4_albedo, corr_co_albedo = corr_co_albedo,
                        corr_ch4_aot = corr_ch4_aot, corr_co_aot = corr_co_aot,
                        prior_term_ch4 = prior_term_ch4, prior_term_co = prior_term_co, prior_mean_impact_ch4 = prior_mean_impact_ch4, prior_mean_impact_co = prior_mean_impact_co,
                        n_points_vec = n_points_vec, area_fraction = fraction_filtered_data, num_overpasses = 1, filtering_flag = 0, mean_eastward_wind = mean_eastward_wind, mean_northward_wind = mean_northward_wind))
    } else {
      # does not pass CO SNR
      return(data.frame(dates_vec = date_ch4, orbit = orbit_ch4,
                      summed_ch4_enhancements_min = NA,  summed_ch4_enhancements_05 = NA, summed_ch4_enhancements_06 = NA, summed_ch4_enhancements_07 = NA, summed_ch4_enhancements_08 = NA, summed_ch4_enhancements_09 = NA,
                        summed_ch4_enhancements_10 = NA, summed_ch4_enhancements_11 = NA, summed_ch4_enhancements_12 = NA, summed_ch4_enhancements_13 = NA, summed_ch4_enhancements_14 = NA, summed_ch4_enhancements_15 = NA,
                        summed_ch4_enhancements_16 = NA, summed_ch4_enhancements_17 = NA, summed_ch4_enhancements_18 = NA, summed_ch4_enhancements_19 = NA, summed_ch4_enhancements_20 = NA, summed_ch4_enhancements_21 = NA,
                        summed_ch4_enhancements_22 = NA, summed_ch4_enhancements_23 = NA, summed_ch4_enhancements_24 = NA, summed_ch4_enhancements_25 = NA,

                      summed_co_enhancements_min = NA, summed_co_enhancements_05 = NA,  summed_co_enhancements_06 = NA, summed_co_enhancements_07 = NA, summed_co_enhancements_08 = NA, summed_co_enhancements_09 = NA,
                        summed_co_enhancements_10 = NA, summed_co_enhancements_11 = NA, summed_co_enhancements_12 = NA, summed_co_enhancements_13 = NA, summed_co_enhancements_14 = NA, summed_co_enhancements_15 = NA,
                        summed_co_enhancements_16 = NA, summed_co_enhancements_17 = NA, summed_co_enhancements_18 = NA, summed_co_enhancements_19 = NA, summed_co_enhancements_20 = NA, summed_co_enhancements_21 = NA,
                        summed_co_enhancements_22 = NA, summed_co_enhancements_23 = NA, summed_co_enhancements_24 = NA, summed_co_enhancements_25 = NA,

                      integ_results_min = NA,  integ_results_05 = NA, integ_results_06 = NA, integ_results_07 = NA, integ_results_08 = NA, integ_results_09 = NA, integ_results_10 = NA,
                        integ_results_11 = NA, integ_results_12 = NA, integ_results_13 = NA, integ_results_14 = NA, integ_results_15 = NA, integ_results_16 = NA, integ_results_17 = NA,
                        integ_results_18 = NA, integ_results_19 = NA, integ_results_20 = NA, integ_results_21 = NA, integ_results_22 = NA, integ_results_23 = NA,
                        integ_results_24 = NA, integ_results_25 = NA,

                        mean_ch4_enhancement_15 = NA, mean_co_enhancement_15 = NA,
                        ch4_mr_bkg_percentile = NA, co_mr_bkg_percentile = NA, corr_vec_mr = NA, corr_vec_enh = NA,
                        ch4_prior_bkg_percentile = NA, co_prior_bkg_percentile = NA,
                        ch4_prior_mean_enh = NA, co_prior_mean_enh = NA, ch4_prior_mean = NA, co_prior_mean = NA, ch4_prior_sd = NA, co_prior_sd = NA,
                        mean_ch4_prec = NA, mean_co_prec = NA, prec_cov = NA,
                        mean_ak_ch4 = NA, mean_ak_co = NA, mean_ch4_albedo = NA, mean_co_albedo = NA, corr_ch4_albedo = NA, corr_co_albedo = NA,
                        corr_ch4_aot = NA, corr_co_aot = NA,
                        prior_term_ch4 = NA, prior_term_co = NA, prior_mean_impact_ch4 = NA, prior_mean_impact_co = NA,
                        n_points_vec = NA, area_fraction = NA, num_overpasses = 0, filtering_flag = 2, mean_eastward_wind = NA, mean_northward_wind = NA))
    }
  } else {
    # does not pass minimum co-located pixel threshold
    return(data.frame(dates_vec = date_ch4, orbit = orbit_ch4,
                      summed_ch4_enhancements_min = NA,  summed_ch4_enhancements_05 = NA, summed_ch4_enhancements_06 = NA, summed_ch4_enhancements_07 = NA, summed_ch4_enhancements_08 = NA, summed_ch4_enhancements_09 = NA,
                        summed_ch4_enhancements_10 = NA, summed_ch4_enhancements_11 = NA, summed_ch4_enhancements_12 = NA, summed_ch4_enhancements_13 = NA, summed_ch4_enhancements_14 = NA, summed_ch4_enhancements_15 = NA,
                        summed_ch4_enhancements_16 = NA, summed_ch4_enhancements_17 = NA, summed_ch4_enhancements_18 = NA, summed_ch4_enhancements_19 = NA, summed_ch4_enhancements_20 = NA, summed_ch4_enhancements_21 = NA,
                        summed_ch4_enhancements_22 = NA, summed_ch4_enhancements_23 = NA, summed_ch4_enhancements_24 = NA, summed_ch4_enhancements_25 = NA,

                      summed_co_enhancements_min = NA, summed_co_enhancements_05 = NA,  summed_co_enhancements_06 = NA, summed_co_enhancements_07 = NA, summed_co_enhancements_08 = NA, summed_co_enhancements_09 = NA,
                        summed_co_enhancements_10 = NA, summed_co_enhancements_11 = NA, summed_co_enhancements_12 = NA, summed_co_enhancements_13 = NA, summed_co_enhancements_14 = NA, summed_co_enhancements_15 = NA,
                        summed_co_enhancements_16 = NA, summed_co_enhancements_17 = NA, summed_co_enhancements_18 = NA, summed_co_enhancements_19 = NA, summed_co_enhancements_20 = NA, summed_co_enhancements_21 = NA,
                        summed_co_enhancements_22 = NA, summed_co_enhancements_23 = NA, summed_co_enhancements_24 = NA, summed_co_enhancements_25 = NA,

                      integ_results_min = NA,  integ_results_05 = NA, integ_results_06 = NA, integ_results_07 = NA, integ_results_08 = NA, integ_results_09 = NA, integ_results_10 = NA,
                        integ_results_11 = NA, integ_results_12 = NA, integ_results_13 = NA, integ_results_14 = NA, integ_results_15 = NA, integ_results_16 = NA, integ_results_17 = NA,
                        integ_results_18 = NA, integ_results_19 = NA, integ_results_20 = NA, integ_results_21 = NA, integ_results_22 = NA, integ_results_23 = NA,
                        integ_results_24 = NA, integ_results_25 = NA,

                        mean_ch4_enhancement_15 = NA, mean_co_enhancement_15 = NA,
                        ch4_mr_bkg_percentile = NA, co_mr_bkg_percentile = NA, corr_vec_mr = NA, corr_vec_enh = NA,
                        ch4_prior_bkg_percentile = NA, co_prior_bkg_percentile = NA,
                        ch4_prior_mean_enh = NA, co_prior_mean_enh = NA, ch4_prior_mean = NA, co_prior_mean = NA, ch4_prior_sd = NA, co_prior_sd = NA,
                        mean_ch4_prec = NA, mean_co_prec = NA, prec_cov = NA,
                        mean_ak_ch4 = NA, mean_ak_co = NA, mean_ch4_albedo = NA, mean_co_albedo = NA, corr_ch4_albedo = NA, corr_co_albedo = NA,
                        corr_ch4_aot = NA, corr_co_aot = NA,
                        prior_term_ch4 = NA, prior_term_co = NA, prior_mean_impact_ch4 = NA, prior_mean_impact_co = NA,
                        n_points_vec = NA, area_fraction = NA, num_overpasses = 0, filtering_flag = 1, mean_eastward_wind = NA, mean_northward_wind = NA))
  }
}

# add in date information
ch4_co_df$datetimes <- as.Date(split_dates(ch4_co_df$dates_vec)) # from split_dates.R
ch4_co_df$year <- as.numeric(format(as.Date(ch4_co_df$datetimes), "%Y"))
ch4_co_df$month <- as.numeric(format(as.Date(ch4_co_df$datetimes), "%m"))

# drop observations that do not pass filtering requirements
ch4_co_df <- ch4_co_df %>% drop_na()


# -------------------------------
# Adding GFAS FIRE FILTERING
# -------------------------------
temp_ch4co_df  <- ch4_co_df
if (dim(temp_ch4co_df)[1] > 0) {
    # -------------------------------
    #  read in city GFAS data after cropping to each city using automate_GFAS_cropping.sh
    # -------------------------------
    GFAS_CO_nc_path <- paste0(path_to_data, "/GFAS/city_GFAS/", domain_code, "_cropped_GFAS.nc")
    # read in netCDF4 data
    ncin <- nc_open(GFAS_CO_nc_path) # var81, time, lon, lat
    time <- ncvar_get(ncin, "valid_time")
    tunits <- ncatt_get(ncin, "valid_time", "units")
    ntime <- dim(time) #1827
    #var_array <- ncvar_get(ncin, "var81")
    # units: hours since 2019-1-1 00:00:00
    # calendar: proleptic_gregorian
    # convert time units
    cf <- CFtime(tunits$value, calendar = "proleptic_gregorian", time) # convert time to CFtime class
    timestamps <- CFtimestamp(cf)

    # -------------------------------
    # read in city data into rast
    # -------------------------------
    gfas_rast <- terra::rast(GFAS_CO_nc_path)
    area <- terra::cellSize(gfas_rast, mask = FALSE, transform = TRUE, unit = "m")
    gfas_rast_emissions <- gfas_rast * area * 3600 / 1000  # kg/s to Mg/h
    # already cropped to city

    urban_domain <- terra::ext(as.vector(cbind(city_latlon[2] - (delta_lon), city_latlon[2] + (delta_lon),
                                                        city_latlon[1] - (delta_lat), city_latlon[1] + (delta_lat))))
    # -------------------------------
    # sum emissions for each area for each day
    # -------------------------------
    buffer_small_square <- rast(urban_domain, res = 0.1)
    buffer_large_square <- rast(urban_domain * 2, res = 0.1)
    # Prepare a vector to store dates when emissions are above the threshold
    fire_frame <- data.frame(timestamps = timestamps, inner = rep(0, length(timestamps)), outer = rep(0, length(timestamps)))
    # iterate over days
    for (i in 1:nlyr(gfas_rast_emissions)) {
        # Extract the current layer
        layer <- gfas_rast_emissions[[i]]
        # Apply the buffer mask to the current layer
        masked_layer_1 <- crop(layer, buffer_large_square) # used mask with circle
        masked_layer_2 <- crop(layer, buffer_small_square)
        # Sum the values within the buffer, checking against threshold
        sum_emissions_1 <- sum(values(masked_layer_1), na.rm = TRUE)
        sum_emissions_2 <- sum(values(masked_layer_2), na.rm = TRUE)
        fire_frame[i, "outer"] <- sum_emissions_1
        fire_frame[i, "inner"] <- sum_emissions_2
    }
    fire_frame$datetimes <- as.Date(fire_frame$timestamps)
    fire_frame <- subset(fire_frame, select = -c(timestamps))

    # -------------------------------
    # set thresholds
    # -------------------------------
    threshold_1 <- 57 # Mg/h
    threshold_2 <- 23 # Mg/h
    # or set range for sensitivity study


    # -------------------------------
    #merge ch4co and fire_frame on datetimes
    # -------------------------------
    merged_fire_ch4co <- left_join(temp_ch4co_df, fire_frame, by = "datetimes")
    merged_fire_ch4co$outer_fire <- merged_fire_ch4co$outer > threshold_1
    merged_fire_ch4co$inner_fire <- merged_fire_ch4co$inner > threshold_2
    outliers <- boxplot.stats(merged_fire_ch4co$co_mr_bkg_percentile)$out
    outliers_ind <- which(merged_fire_ch4co$co_mr_bkg_percentile %in% c(outliers))
    merged_fire_ch4co$outlier_co_bkg <- merged_fire_ch4co[, "dates_vec"] %in% merged_fire_ch4co[outliers_ind, "dates_vec"]

    # variables to plot: integ_results_15 ,temp_mr_bkg_percentile
    num_overpasses_original <- dim(merged_fire_ch4co)[1]
    ratio_original <- mean(merged_fire_ch4co$integ_results_15)

    # -------------------------------
    # filter and calculate percent of data removed
    # -------------------------------
    merged_fire_ch4co_filtered <- merged_fire_ch4co[(merged_fire_ch4co$outer_fire == "FALSE") & (merged_fire_ch4co$inner_fire == "FALSE"), ]
    num_overpasses_filtered <- dim(merged_fire_ch4co_filtered)[1]
    ratio_filtered <- mean(merged_fire_ch4co_filtered$integ_results_15)
    data_passed_fire_filtering <- num_overpasses_filtered / num_overpasses_original * 100
    ratio_change <- (ratio_filtered - ratio_original) / ratio_original * 100 #ratio is average across all time

    # save filtered TROPOMI data, drop columns that were added
    ch4_co_df <- subset(merged_fire_ch4co_filtered, select = -c(inner, outer, outer_fire, inner_fire))
}
# fire filtering plotting and ratios not adapted for being in this file
paste(city, " -> num overpasses that pass filtering: ", sum(ch4_co_df$num_overpasses, na.rm = TRUE), "out of", dim(city_orbits)[1])



# -------------------------------
# save plots
# -------------------------------

if (save_plots) {
  # find all pdfs in temp_plots file
  temp_list <- list.files(path = paste0(path_to_plots, "temp_plots/"), full.names = TRUE)
  # combine
  if (length(temp_list > 0)) {
    qpdf::pdf_combine(input = temp_list, output = paste(path_to_plots, city, "_overpasses_", Sys.Date(), "_", pixel_coverage, ".pdf", sep = ""))
  #dev.off()
  }
  # remove files in temp folder
  file.remove(temp_list)
}

# -------------------------------
# save data
# -------------------------------
#if (save_data) {
##  save.image(paste0(path_prefix, "Data/", city, type_suffix, Sys.Date(), ".RData")) #added box vs urb to this using plot_type ('box', 'urb')
#}

write.table(ch4_co_df, paste0(path_to_data, "/ch4co/", city, "_ch4co.csv"), row.names = FALSE, sep = ",")