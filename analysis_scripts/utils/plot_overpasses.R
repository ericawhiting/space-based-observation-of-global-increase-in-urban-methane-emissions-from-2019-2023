# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(ggplot2)
library(ncdf4)
library(dplyr)
library(cowplot)
library(stringr)
library(scales)

#path for output plots
path_to_plots <- "../plots/results/" # change your path here

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))
source(paste0(path_to_utils, "calc_dry_column_mixing_ratio.R"))

# -------------------------------
### source: G.Plant Tropomi_ch4ch_city_analysis_single_day.R
DataFrameforPlot2 <- function(inputArray,
                              lat_f, lon_f,
                              latb1, latb2, latb3, latb4,
                              lonb1, lonb2, lonb3, lonb4) {
  ids <- factor(1:length(as.vector(lat_f)))
  lat_bounds_f <- c(rbind(latb1, latb2, latb3, latb4))
  lon_bounds_f <- c(rbind(lonb1, lonb2, lonb3, lonb4))
  df_out <- data.frame(id = rep(ids, each = 4), value = rep(as.vector(inputArray), each = 4),
                       x = as.vector(lon_bounds_f),
                       y = as.vector(lat_bounds_f),
                       lat = rep(as.vector(lat_f), each = 4), lon = rep(as.vector(lon_f), each = 4))
  return(df_out)
}

PlotRegion2 <- function(df, latlon, title = "", units = "", q = c(0, 1), printPlot = FALSE,
                        box.coords = c(NaN, NaN, NaN, NaN), urb_mask_in = NaN, cmap = "RdYlBu") {
  # Plot the given dataset over a geographic region.
  # Best to show box method
  #
  # Args:
  #   df: The dataset, should include the gas data, lat, lon columns, along with bounds
  #   title: The plot title
  #   units: for colorbar label
  #   q: vector of quantiles for colobar, ex. q=c(0.5,0.99)

  subtitle <- paste("Data: MIN =", formatC(min(df$value, na.rm = TRUE), format = "f", digits = 1),
                   "MAX =", formatC(max(df$value, na.rm = TRUE), format = "f", digits = 1))

   p <- ggplot(df, aes(x = x, y = y)) + geom_polygon(aes(fill = value, group = id)) +
        borders("world", xlim = range(df$lon), ylim = range(df$lat),
            colour = "black", linewidth = 0.7, fill = NA) +
    theme_bw() +
    theme(panel.ontop = FALSE, panel.background = element_blank()) +
    scale_fill_distiller(palette = cmap,
                         limits = c(quantile(df$value, q[1], na.rm = TRUE),
                                  quantile(df$value, q[2], na.rm = TRUE)), oob = scales::squish) +
    coord_quickmap(xlim = c(latlon[3], latlon[4]), ylim = c(latlon[1], latlon[2])) +
    geom_rect(aes(xmin = box.coords[3],
                  xmax = box.coords[4],
                  ymin = box.coords[1],
                  ymax = box.coords[2]),
              fill = "transparent", color = "black", linewidth = 1.5) +
    labs(title = title, subtitle = subtitle,
         x = "Longitude", y = "Latitude",
         fill = units)
    if (printPlot) {
        print(p)
    }
  return(p)
}


plot_satellite_overpass <- function(city, orbit_number, title_key_1, title_key_2) {

  # -------------------------------
  # Assign Filters
  # -------------------------------
  pixel_coverage <- 0.2 # 20% area coverage


  land_classification_filter <- TRUE # both CO and CH4
  co_scatter_filter <- TRUE
  qa_filter <- TRUE #both CO and CH4

  #pick ONLY one
  ch4_hu_amt_2016_filter <- FALSE # filters from H. Hu AMT 2016
  ch4_exp_filter <- TRUE # experimental filters from SRON (not official!)

  #CO filter thresholds
  co_qa_filter <- 0 # qa >0
  prec_filt_co <- 5 # precision filter
  t_threshold  <- 0.5  # clear sky optical thickness threshold
  cloud_height_threshold <- 5000 # clear sky cloud height threshold

  #CH4 filters thresholds
  ch4_qa_filter <- 0 # qa >0
  prec_filt_ch4 <- 10 # precision filter
  ch4_aero <- 0.07 #SWIR aerosol optical thickness
  ch4_albedo_limit <- 0.02

  #Define box sizes
  city_latlon <- find_city_center_from_pop_dens(city)
  delta_lat <- read_ideal_boxsize_lat(city)
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

  city_list_path <- "/backup/erwh/TROPOMI_ch4co/city_query_lists/cropped_tropomi/"
  city_list_fname <- paste0(city_list_path, city, "_files.csv")

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
  i <- which(co_orbit_list == orbit_number)
  #CH4
  nc_pathname_ch4 <- city_files[i, 1] # has path in it
  nc_filename_ch4 <- basename(city_files[i, 1])
  nc_ch4 <- nc_open(nc_pathname_ch4)
  date_ch4 <- substr(nc_filename_ch4, 21, 28)

  start_time_ch4 <- substr(nc_filename_ch4, 30, 35)
  orbit_ch4 <- substr(nc_filename_ch4, 53, 57)

  #CO
  nc_pathname_co <- city_files[i, 2] # has path in it
  nc_filename_co <- basename(city_files[i, 1])
  nc_co <- nc_open(nc_pathname_co)
  date_co <- substr(nc_filename_co, 21, 28)
  start_time_co <- substr(nc_filename_co, 30, 35)
  orbit_co <- substr(nc_filename_co, 53, 57)

  dates_vec <- date_ch4[[1]]

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


  # -------------------------------
  # Crop to city
  # -------------------------------
  # cropping to larger than the city to show where urban domain lies, check orbit at true city domain size to see if it will pass filtering!!
  city_coords <- cbind(city_latlon[1] - (delta_lat + 0.15), city_latlon[1] + (delta_lat + 0.15),
                         city_latlon[2] - (delta_lon + 0.15 * lat_adj), city_latlon[2] + (delta_lon + 0.15 * lat_adj))
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
    #print(paste(length(indFilt),' grid cells removed by CO filter'))
  }

  if (ch4_hu_amt_2016_filter) {
    ind_filt2 <- which(ch4_opt_thickness_city * ch4_aero_height_city / ch4_aero_size_city > 120) # threshold from H. Hu 2016 AMT
    # left off here
    ch4_mr_city[ind_filt2] <- NA
    ch4_prec_city[ind_filt2] <- NA
    #print(paste(length(indFilt2),' grid cells removed by CH4 aerosol filter'))

    ind_filt3 <- which(ch4_albedo_city < 0.02) # threshold from H. Hu 2016 AMT
    ch4_mr_city[ind_filt3] <- NA
    ch4_prec_city[ind_filt3] <- NA
    #print(paste(length(indFilt3),' grid cells removed by CH4 albedo filter'))

  }

  if (ch4_exp_filter) {
    #SWIR aerosol optical thickness and precision
    ind_filt6 <- which(ch4_opt_thickness_city >= ch4_aero | ch4_prec_city > prec_filt_ch4)
    ch4_mr_city[ind_filt6] <- NA
    ch4_prec_city[ind_filt6] <- NA
    #print(paste(length(indFilt6),' grid cells removed by exp. CH4 aerosol filter'))
    #Note qa filter is below

    ind_filt3 <- which(ch4_albedo_city < ch4_albedo_limit) # threshold from H. Hu 2016 AMT
    ch4_mr_city[ind_filt3] <- NA
    ch4_prec_city[ind_filt3] <- NA
    #print(paste(length(indFilt3),' grid cells removed by CH4 albedo filter'))
  }

  if (qa_filter) {
    #qa filters
    #CH4
    ind_filt4 <- which(ch4_qa_city <= ch4_qa_filter)
    ch4_mr_city[ind_filt4] <- NA
    ch4_prec_city[ind_filt4] <- NA
    #print(paste(length(ind_filt4), " grid cells CH4 qa filters"))
    #CO
    ind_filt5 <- which(co_qa_city <= co_qa_filter)
    co_mr_city[ind_filt5] <- NA
    co_prec_city[ind_filt5] <- NA
    #print(paste(length(ind_filt5), " grid cells CO qa filters"))
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

  # -------------------------------land_water_classification_ch4 = land_water_classification_ch4,
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
                              ch4_lat_b1 = ch4_lat_b1, ch4_lat_b2 = ch4_lat_b2, ch4_lat_b3 = ch4_lat_b3, ch4_lat_b4 = ch4_lat_b4, ch4_lon_b1 = ch4_lon_b1, ch4_lon_b2 = ch4_lon_b2, ch4_lon_b3 = ch4_lon_b3, ch4_lon_b4 =  ch4_lon_b4,
                              co_lat_b1 = co_lat_b1, co_lat_b2 = co_lat_b2, co_lat_b3 = co_lat_b3, co_lat_b4 = co_lat_b4, co_lon_b1 = co_lon_b1, co_lon_b2 = co_lon_b2, co_lon_b3 = co_lon_b3, co_lon_b4 =  co_lon_b4)
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
                                      ch4_spatial_df, by = c("ch4_lat", "ch4_lon", "ch4_lat_b1", "ch4_lat_b2", "ch4_lat_b3", "ch4_lat_b4", "ch4_lon_b1", "ch4_lon_b2", "ch4_lon_b3", "ch4_lon_b4"), all.y = TRUE)

    filtered_filled_data_co <- merge(filtered_data %>%
                                         select(co, co_lat, co_lon, co_lat_b1, co_lat_b2, co_lat_b3, co_lat_b4, co_lon_b1, co_lon_b2, co_lon_b3, co_lon_b4),
                                     co_spatial_df, by = c("co_lat", "co_lon", "co_lat_b1", "co_lat_b2", "co_lat_b3", "co_lat_b4", "co_lon_b1", "co_lon_b2", "co_lon_b3", "co_lon_b4"), all.y = TRUE)
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


    # save pixel mixing ratio backgrounds calculated with 17 percentile background
    ch4_mr_bkg_percentile <- mr_bkg_percentile[14, 1]
    co_mr_bkg_percentile <- mr_bkg_percentile[14, 2]

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

    filtered_data$ch4_numerator <- ch4_numerator
    filtered_data$co_denominator <- co_denominator
    # summed enhancements for ch4 and co: one summed ch4 enhancement and one summed co enhancement per orbit
    summed_ch4_enhancements <- colSums(ch4_numerator) # length of 21
    summed_co_enhancements <- colSums(co_denominator) # length of 21
    # -------------------------------
    # SNR Filtering with CO
    # -------------------------------
    n_points_vec <- dim(filtered_data)[1]
    # larger domain based on city_coords mask, set extra buffer added in to 0 to look only at urban domain
    if (summed_co_enhancements[12] > (n_points_vec * 7)) { # checks 15% background value
      # ratio of summed enhancements
      ch4co_ratio_of_summed_enhancements <- summed_ch4_enhancements / summed_co_enhancements # dim: city_orbits, background options

      # correlation coef.
      corr_vec_mr <- cor(filtered_data$ch4, filtered_data$co) # correlation of mixing ratio across the overpass at pixel level, length 1
      corr_vec_enh <- diag(cor(ch4_numerator, co_denominator)) # correlation of the enhancements at pixel level, length 5

      # -------------------------------
      # ALBEDO vs MR:
      # -------------------------------
      co_albedo_filtered <- data.frame(albedo = filtered_data$co_albedo)
      co_albedo_filtered$dates <- dates_vec
      co_albedo_filtered$orbit <- orbit_co
      ch4_albedo_filtered <- data.frame(albedo = filtered_data$ch4_albedo)
      ch4_albedo_filtered$dates <- dates_vec
      ch4_albedo_filtered$orbit <- orbit_ch4

      corr_ch4_albedo <- cor(filtered_data$ch4_albedo, filtered_data$ch4, use = "complete.obs")
      corr_co_albedo <- cor(filtered_data$co_albedo, filtered_data$co, use = "complete.obs")

      # plot mixing ratios
      # is zoomed out compared to ideal city boxsize based on adjustment in making city_coords
      ch4_df2 <- DataFrameforPlot2(filtered_filled_data_ch4$ch4, filtered_filled_data_ch4$ch4_lat, filtered_filled_data_ch4$ch4_lon,
                                   filtered_filled_data_ch4$ch4_lat_b1, filtered_filled_data_ch4$ch4_lat_b2, filtered_filled_data_ch4$ch4_lat_b3, filtered_filled_data_ch4$ch4_lat_b4,
                                   filtered_filled_data_ch4$ch4_lon_b1, filtered_filled_data_ch4$ch4_lon_b2, filtered_filled_data_ch4$ch4_lon_b3, filtered_filled_data_ch4$ch4_lon_b4)

      co_df2 <- DataFrameforPlot2(filtered_filled_data_co$co, filtered_filled_data_co$co_lat, filtered_filled_data_co$co_lon,
                                  filtered_filled_data_co$co_lat_b1, filtered_filled_data_co$co_lat_b2, filtered_filled_data_co$co_lat_b3, filtered_filled_data_co$co_lat_b4,
                                  filtered_filled_data_co$co_lon_b1, filtered_filled_data_co$co_lon_b2, filtered_filled_data_co$co_lon_b3, filtered_filled_data_co$co_lon_b4)

      q <- c(0.05, 0.99)
      plot_theme_text <- theme(panel.grid.minor = element_blank(),
                  panel.border = element_blank(),
                  plot.title = element_text(size = 20, face = "bold", family = "Times"),
                  axis.text.y = element_text(size = 16, family = "Times"),
                  axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
                  axis.text.x = element_text(size = 16, family = "Times"),
                  axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
                  legend.title = element_text(size = 16, face = "bold", family = "Times"),
                  legend.text = element_text(size = 16, face = "bold", family = "Times"),
                  text = element_text(family = "Times"),
                  strip.text = element_text(size = 16, family = "Times"))

      domain_city_coords <- cbind(city_latlon[1] - (delta_lat), city_latlon[1] + (delta_lat),
                         city_latlon[2] - (delta_lon), city_latlon[2] + (delta_lon))
      p1 <- PlotRegion2(ch4_df2, city_coords, title = paste("CH4", date_ch4), units = "ppb", q = q, printPlot = FALSE, box.coords = city_coords)
      p2 <- PlotRegion2(co_df2, city_coords,  title = paste("CO", date_co), units = "ppb", q = q, printPlot = FALSE, box.coords = city_coords)
      date_format <- as.Date(as.character(date_ch4), format = "%Y%m%d")

      date_str <- format(date_format, "%Y-%m-%d")
      p1 <- p1 +
            theme_linedraw() +
            plot_theme_text +
            scale_x_continuous(expand = c(0, 0), breaks = seq(floor(city_coords[3]), city_coords[4], 0.5), labels = function(l) sprintf("%.1f°", l)) +
            scale_y_continuous(expand = c(0, 0), breaks = seq(floor(city_coords[1]), city_coords[2], 0.5), labels = function(l) sprintf("%.1f°", l)) +
            labs(fill = "CH4 [ppb]") +
            geom_rect(aes(xmin = domain_city_coords[3],
                          xmax = domain_city_coords[4],
                          ymin = domain_city_coords[1],
                          ymax = domain_city_coords[2]),
                      fill = "transparent", color = "black", linewidth = 1.5) +
              labs(subtitle = NULL, title = bquote(bold("(" * .(title_key_1) * ") TROPOMI CH" [4] * " " * .(date_str)))) +
              theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
              title = element_text(size = 18),
              axis.title.x = element_text(size = 18, hjust = 0.5), axis.title.y = element_text(size = 18, hjust = 0.5),
              legend.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(size = 1)) +
              theme(plot.margin = unit(c(0.05, 0.25, 0.05, 0.1), "inches")) +
              theme(plot.title.position = "plot") +
              guides(fill = guide_colorbar(
                      ticks.linewidth = 1.5,      # Linewidth of ticks
                      ticks.colour = "black",     # Colour of ticks
                      frame.colour = "black",     # Colour of the border frame
                      frame.linewidth = 1.5,      # Linewidth of the border frame
                      draw.ulim = TRUE,          # Draw upper limit tick mark
                      draw.llim = TRUE,          # Draw lower limit tick mark
                      label.theme = element_text(size = 16), # Labels size
                      title.theme = element_text(size = 16),
                      barwidth = unit(0.075, "npc"), # Bar width (as a fraction of the plotting area)
                      barheight = unit(0.5, "npc") # Bar height (as a fraction of the plotting area)
                      ))       # Adjust the size if needed


      p2 <- p2 +
            theme_linedraw() +
            plot_theme_text +
            scale_x_continuous(expand = c(0, 0), breaks = seq(floor(city_coords[3]), city_coords[4], 0.5), labels = function(l) sprintf("%.1f°", l)) +
            scale_y_continuous(expand = c(0, 0), breaks = seq(floor(city_coords[1]), city_coords[2], 0.5), labels = function(l) sprintf("%.1f°", l)) +
            labs(fill = "CO [ppb]") +
            geom_rect(aes(xmin = domain_city_coords[3],
                          xmax = domain_city_coords[4],
                          ymin = domain_city_coords[1],
                          ymax = domain_city_coords[2]),
                      fill = "transparent", color = "black", linewidth = 1.5) +
              labs(subtitle = NULL, title =  bquote(bold("(" * .(title_key_2) * ") TROPOMI CO" * " " * .(date_str)))) +
              theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
              title = element_text(size = 18),
              axis.title.x = element_text(size = 18, hjust = 0.5), axis.title.y = element_text(size = 18, hjust = 0.5),
              legend.margin = margin(10, 10, 10, 10),
              panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
              axis.ticks = element_line(size = 1)) +
              theme(plot.margin = unit(c(0.05, 0.25, 0.05, 0.1), "inches")) +
              theme(plot.title.position = "plot") +
              guides(fill = guide_colorbar(
                      ticks.linewidth = 1.5,      # Linewidth of ticks
                      ticks.colour = "black",     # Colour of ticks
                      frame.colour = "black",     # Colour of the border frame
                      frame.linewidth = 1.5,      # Linewidth of the border frame
                      draw.ulim = TRUE,          # Draw upper limit tick mark
                      draw.llim = TRUE,          # Draw lower limit tick mark
                      label.theme = element_text(size = 16), # Labels size
                      title.theme = element_text(size = 16),
                      barwidth = unit(0.075, "npc"), # Bar width (as a fraction of the plotting area)
                      barheight = unit(0.5, "npc") # Bar height (as a fraction of the plotting area)
                      ))
      legend_p1 <- get_plot_component(p1, "guide-box", return_all = TRUE)
      legend_plot_p1 <- cowplot::ggdraw(legend_p1[[1]])
      legend_p2 <- get_plot_component(p2, "guide-box", return_all = TRUE)
      legend_plot_p2 <- cowplot::ggdraw(legend_p2[[1]])
      ggsave(paste0(city, "_ch4_legend.png"), legend_plot_p1, path = path_to_plots,  height = 3.5, width = 1.6, dpi = 320,  bg = "white")
      ggsave(paste0(city, "_co_legend.png"), legend_plot_p2, path = path_to_plots, height = 3.5, width = 1.6, dpi = 320,  bg = "white")


      p1 <- p1 + guides(fill = "none")
      p2 <- p2 + guides(fill = "none")
      ggsave(paste0(city, "_CH4mr.png"), plot = p1, device = "png", path = path_to_plots, width = 5, height = 5, dpi = 320, bg = "white")
      ggsave(paste0(city, "_COmr.png"), plot = p2, device = "png", path = path_to_plots, width = 5, height = 5, dpi = 320, bg = "white")

        #  dry column mixing ratio
      return()
    }
  }
}
