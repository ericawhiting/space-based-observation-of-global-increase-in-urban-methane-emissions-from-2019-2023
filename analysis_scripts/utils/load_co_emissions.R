# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(terra)
library(stringr)
library(dplyr)
library(abind)

path_to_data <- "../Data/" # change your path here

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))
# ------------------------------------------------------------------------------------------------------------------
# -----------
# USES MONTHLY EDGAR CO FILES, CHECK FOR RECENT CHANGES TO FILES !!
# -----------

assign_season_year <- function(year, month) {
  #' assign a season (DJF, MAM, JJA, SON) and year XXXX to year, month combination
  #' @param year
  #' @param month c(1-12) or c(01-12)
  #' @return season_year in format XXXX_MMM (ex: 2019_DJF)
  # --------
  if (month == 1 || month == 2) {
    season <- "DJF"
    season_year <- paste(year, season, sep = "_")
  } else if (month == 12) { # next years winter
    season <- "DJF"
    season_year <- paste(year + 1, season, sep = "_")
  } else if (month %in% c(3, 4, 5)) {
    season <- "MAM"
    season_year <- paste(year, season, sep = "_")
  } else if (month %in% c(6, 7, 8)) {
    season <- "JJA"
    season_year <- paste(year, season, sep = "_")
  } else if (month %in% c(9, 10, 11)) {
    season <- "SON"
    season_year <- paste(year, season, sep = "_")
  }
  return(season_year)
}

open_edgarv8.1_monthly_file <- function(edgar_y, sector) { # nolint: object_name_linter.
  #' open netcdf4 files for EDGARv8.1 given year, sector, and return opened variable
  #' @param edgar_y is year to pull for EDGAR file names (EDGARv8.1 ends in 2022)
  #' @param sector EDGARv8.1 CO sectors: "AGRICULTURE", "BUILDINGS", "FUEL_EXPLOITATION", "IND_COMBUSTION", "POWER_INDUSTRY", "TRANSPORT", "WASTE"
  #' @return returns an array of emissions (rate not a flux!) of CO emissions in Mg/month from EDGAR for month in year for sector for whole globe
  # --------
  fname <- paste0(path_to_data, "EDGAR/co/monthly_sector_specific_gridmaps/v8.1/", "v8.1_FT2022_AP_CO_", edgar_y, "_bkl_", sector, "_emi.nc")
  nc_file <- ncdf4::nc_open(fname)
  var_name <- names(nc_file$var)
  var <- ncdf4::ncvar_get(nc_file, var_name)
  ncdf4::nc_close(nc_file)
  return(var)
}

return_cell_values <- function(var, urban_domain_rast, m, sf_flag) {
  #' given a nc var from EDGARv8.1 CO monthly, sectoral file and a given domain, return emissions values from cropped domain
  #' @param var variable name from ncdf4 file
  #' @param urban_domain_rast the terra extent of the spatial domain
  #' @param m month
  #' @param sf_flag if the urban domain outline is not a rast, it will be a sf (csf outline) and a mask must be applied in cropping
  # --------
  var_1 <- var[, ncol(var):1, m] # ncvar_get reads in as lat lon, want lon lat # nolint: seq_linter.
  var_2 <- t(var_1)
  r_grid <- terra::rast(var_2, crs = "EPSG:4326")
  terra::ext(r_grid) <- c(-180, 180, -90, 90)
  cropped_gas <- terra::crop(r_grid, urban_domain_rast, mask = sf_flag)
  emissions_rast <- values(cropped_gas) # Mg/month or Tonne/month
  return(emissions_rast)
}

generate_random_samples <- function(mean_value, n_samples = 1000, sigma) {
  #' randomly sample from normal distribution 1000 times where mean is set by cell value and sigma is given
  #' @param mean_value the center of the normal distribution from which to sample
  #' @param n_samples the number of times to randomly select, default value is 1000
  #' @param sigma the standard deviation of the normal distribution from which to sample, a value not a fraction or percent
  #' @return random values of length n_samples, randomly sampled from normal distribution
  # --------
  rnorm(n_samples, mean = mean_value, sd = sigma)
}

load_all_month_sectors_co <- function(y, urban_domain_rast, sf_flag) {
  #' load and sum all edgar co values across 8 sectors for each month for each grid cell, combine into dataframe for year
  #' @param y year to pull edgar file (v8.1 goes until 2022)
  #' @param urban_domain_rast the terra extent of the spatial domain
  #' @param sf_flag if the urban domain outline is not a rast, it will be a sf (csf outline) and a mask must be applied in cropping within return_cell_values
  #' @return values_of_urban_rast monthly emissions in Tonnes or Mg per month of dimension # gridcells by # months in year of file
  # --------
  var_1 <- open_edgarv8.1_monthly_file(y, "AGRICULTURE")
  var_2 <- open_edgarv8.1_monthly_file(y, "BUILDINGS")
  var_3 <- open_edgarv8.1_monthly_file(y, "FUEL_EXPLOITATION")
  var_4 <- open_edgarv8.1_monthly_file(y, "IND_COMBUSTION")
  var_5 <- open_edgarv8.1_monthly_file(y, "IND_PROCESSES")
  var_6 <- open_edgarv8.1_monthly_file(y, "POWER_INDUSTRY")
  var_7 <- open_edgarv8.1_monthly_file(y, "TRANSPORT")
  var_8 <- open_edgarv8.1_monthly_file(y, "WASTE")

  monthly_sum_across_sectors <- var_1 + var_2 + var_3 + var_4 + var_5 + var_6 + var_7 + var_8
  month_list <- seq(1, 12)
  values_of_urban_rast <- c()
  for (m in month_list) { # gives array where columns are each month and rows are values of grid cells
      values_of_urban_rast <-  cbind(values_of_urban_rast, return_cell_values(monthly_sum_across_sectors, urban_domain_rast, m, sf_flag))
  }
  return(values_of_urban_rast) # monthly emissions are Mg / month, Tonnes / month
}

load_all_months_sectors_individually_co <- function(city, y, urban_domain_rast, sf_flag = FALSE) {
  var_1 <- open_edgarv8.1_monthly_file(y, "AGRICULTURE")
  var_2 <- open_edgarv8.1_monthly_file(y, "BUILDINGS")
  var_3 <- open_edgarv8.1_monthly_file(y, "FUEL_EXPLOITATION")
  var_4 <- open_edgarv8.1_monthly_file(y, "IND_COMBUSTION")
  var_5 <- open_edgarv8.1_monthly_file(y, "IND_PROCESSES")
  var_6 <- open_edgarv8.1_monthly_file(y, "POWER_INDUSTRY")
  var_7 <- open_edgarv8.1_monthly_file(y, "TRANSPORT")
  var_8 <- open_edgarv8.1_monthly_file(y, "WASTE")
  version <- "v8.1"
  month_list <- seq(1, 12)
  values_of_urban_rast <- data.frame()
  for (m in month_list) { # gives array where columns are each month and rows are values of grid cells
      values_of_urban_rast <-  rbind(values_of_urban_rast,
                                     data.frame(city = city, year = y, month = m, version = version,
                                        AGRICULTURE = colSums(return_cell_values(var_1, urban_domain_rast, m, sf_flag)),
                                        BUILDINGS = colSums(return_cell_values(var_2, urban_domain_rast, m, sf_flag)),
                                        FUEL_EXPLOITATION = colSums(return_cell_values(var_3, urban_domain_rast, m, sf_flag)),
                                        IND_COMBUSTION = colSums(return_cell_values(var_4, urban_domain_rast, m, sf_flag)),
                                        IND_PROCESSES = colSums(return_cell_values(var_5, urban_domain_rast, m, sf_flag)),
                                        POWER_INDUSTRY = colSums(return_cell_values(var_6, urban_domain_rast, m, sf_flag)),
                                        TRANSPORT = colSums(return_cell_values(var_7, urban_domain_rast, , sf_flag)),
                                        WASTE = colSums(return_cell_values(var_8, urban_domain_rast, m, sf_flag))))
  }
  row.names(values_of_urban_rast) <- NULL
  return(values_of_urban_rast) # monthly emissions are Mg / month, Tonnes / month
}

load_htap_co_emissions <- function(inv_version, year) {
  # not a FLUX, in tons/y
  HTAP_data_path <- paste0(path_to_data, "HTAPv3/annual/")
  fname <- paste0(HTAP_data_path, "edgar_HTAPv3_", year, "_CO.nc") # annual file name

  nc_file <- ncdf4::nc_open(fname) # can use str(nc_file$dim) to look at properties
  var_name <- names(nc_file$var)
  var <- 0
  for (i in seq_len(14)) {
    var <- var + ncdf4::ncvar_get(nc_file, var_name[i]) # metric ton/yr
  }
  var <- var[, ncol(var):1] # ncvar_get reads in as lat lon, want lon lat
  var <- t(var)
  ncdf4::nc_close(nc_file)
  r_grid <- terra::rast(var, crs = "EPSG:4326")
  terra::ext(r_grid) <- c(-180, 180, -90, 90)
  return(r_grid)
}

scale_edgar_by_htap <- function(urban_domain_rast, sf_flag = FALSE) {
  #' use city total HTAPv3 to scale city total EDGAR v8.1 CO sum for cities that have unrealistic CO values (shg)
  #' @param urban_domain_rast terra rast to crop edgar and htap grids down to in order to find city sum
  #' @param sf_flag if the urban domain outline is not a rast, it will be a sf (csf outline) and a mask must be applied in cropping within return_cell_values
  #' @return co_adjustment_scale fraction to multiply by edgar value in find_co_emissions_from_monthly
  # --------
    # Load EDGAR (Tonnes/month, aggregate to year is Tonnes/y)
    edgar_2018_city_monthly <- load_all_month_sectors_co(2018, urban_domain_rast, sf_flag) # dim: # grid cells by 12 months; unit: Tonnes/month
    edgar_2018_city <- rowSums(edgar_2018_city_monthly) # dim: # grid cells; unit: Tonnes/yr
    # find co emissions for whole city for 2019 from EDGAR CO
    edgar_2018 <- sum(edgar_2018_city, na.rm = TRUE)
    # Load HTAP (Tonnes/y)
    htap_2018_globe <- load_htap_co_emissions("HTAPv3_CO", 2018)
    htap_2018_city <- terra::crop(htap_2018_globe, urban_domain_rast, mask = sf_flag)
    htap_2018 <- terra::global(htap_2018_city, "sum", na.rm = TRUE)$sum

    co_adjustment_scale <- htap_2018 / edgar_2018 # used in Shanghai only
    return(co_adjustment_scale)
}

combine_monthly_emissions <- function(co_monthly_emissions_array, year_month_df, obs_year_month_list, aggregate, all_months) { # nolint: cyclocomp_linter.
  #' find average EDGAR CO emissions from monthly files for given time step and with given obs constraint
  #' @param co_monthly_emissions_array df of co emissions from monthly edgar files with dim # gridcells for domain vs # all months across years of study
  #' @param year_month_df df of all months from start to end of study with columns year, month
  #' @param obs_year_month_list df of unique months of observations with columns year, month, season_year
  #' @param aggregate annual, seasonal, monthly: the time step over which we aggregate and find mean emissions
  #' @param all_months TRUE or FALSE, if TRUE then looks at all months of given time step aggregation, if false, uses months in which we had obs
  #' @return co_emissions_domain df of mean co emissions for given timesteps in kg/s
  # --------
  co_emissions_domain <- data.frame()
  if (aggregate == "annual") {
    unique_years <- unique(year_month_df$year)
    for (year in unique_years) {
        # find indicies for unique year
        year_indices <- which(year_month_df$year == year)
        # subset aray for year
        subset_array <- co_monthly_emissions_array[, year_indices, drop = FALSE]
        if (all_months != TRUE) { # subset for months in which we have observations
            subset_array <- subset_array[, obs_year_month_list[obs_year_month_list$year == year, "month"]$month, drop = FALSE]
        }
        # mean across desired months for each grid cell
        cell_mean <- apply(subset_array, c(1), mean) * 12 * 1000 / (365.25 * 24 * 3600)
        domain_mean <- sum(cell_mean, na.rm = TRUE)
        co_emissions_domain <- rbind(co_emissions_domain, data.frame(year, domain_mean))
    }
  } else if (aggregate == "seasonal") {
    year_month_df$season_year <- mapply(assign_season_year, year_month_df$year, year_month_df$month)
    # find unique seasons
    unique_seasons <- unique(year_month_df$season_year)
    for (season_year in unique_seasons) {
        # find indices for unique season and year
        season_indices <- which(year_month_df$season_year == season_year)
        subset_array <- co_monthly_emissions_array[, season_indices, drop = FALSE]
        if (all_months != TRUE) { # subset for months in which we have observations
            months_season <- year_month_df[year_month_df$season_year == season_year, c("year", "month", "season_year")]
            obs_season <- obs_year_month_list[(obs_year_month_list$season_year == season_year), ]
            matching_indices <- which((months_season$year %in% obs_season$year) &
                                      (months_season$month %in% obs_season$month))
            subset_array <- subset_array[, matching_indices, drop = FALSE]
        }
        cell_mean <- apply(subset_array, c(1), mean) * 12 * 1000 / (365.25 * 24 * 3600)
        domain_mean <- sum(cell_mean, na.rm = TRUE)
        co_emissions_domain <- rbind(co_emissions_domain, data.frame(season_year, domain_mean))
    }
  } else if (aggregate == "monthly") {
    unique_years <- unique(year_month_df$year)
    for (year in unique_years) {
        for (month in seq(1, 12)) {
            # find indicies for unique year
            year_month_indices <- which((year_month_df$year == year) & (year_month_df$month == month))
            # subset aray for year and month
            subset_array <- co_monthly_emissions_array[, year_month_indices, drop = FALSE]
            if ((all_months != TRUE) && (dim(obs_year_month_list[(obs_year_month_list$year == year) & (obs_year_month_list$month == month), ])[1] == 0)) {
                subset_array <- matrix(NA, c(dim(subset_array)))  # replace with NAN
            }
            cell_mean <- subset_array * 12 * 1000 / (365.25 * 24 * 3600)
            domain_mean <- sum(cell_mean, na.rm = TRUE)
            co_emissions_domain <- rbind(co_emissions_domain, data.frame(year, month, domain_mean))
        }
    }
  } else {
      print("Aggregation method not developed")
  }
  # return realizations of co emissions with random noise added at the grid level
  return(co_emissions_domain)
}

find_co_emissions_from_monthly <- function(city, urban_domain_rast, year_list, obs_year_month_list, aggregate, all_months, uncertainty_method, sf_flag = FALSE) {
  #' uses edgarv8.1 monthly co files to list out numerical realizations of mean co emissions for given time aggregation level for given domain with errror
  #' @param city city code from which to search for EDGAR uncertainty value
  #' @param urban_domain_rast the terra extent of the spatial domain
  #' @param year_list a list of years to evaluate edgar over
  #' @param obs_year_month_list df of unique months of observations with columns year, month, season_year
  #' @param aggregate annual, seasonal, monthly: the time step over which we aggregate and find mean emissions
  #' @param all_months TRUE or FALSE, if TRUE then looks at all months of given time step aggregation, if false, uses months in which we had obs
  #' @param uncertainty_method full_uncertainty or SRON_uncertainty, for adding noise to CO emissions, must be a fraction
  #' @param sf_flag if the urban domain outline is not a rast, it will be a sf (csf outline) and a mask must be applied in cropping within return_cell_values
  #' default sf_flag is FALSE
  #' @return co_emissions_realizations df of relevant time steps and numerical realizations of co emissions for that time period + noise
  # --------
  arrays_list <- list()
  year_month_df <- data.frame(year = numeric(), month = numeric())
  for (y in year_list) {
    if (y == 2023 || y == 2024) {
      edgar_y <- 2022 # run 2022 for edgar
    } else {
      edgar_y <- y
    }
    # 2d array with dim: #grid cell, #month in year
    array_for_year <- load_all_month_sectors_co(edgar_y, urban_domain_rast, sf_flag)
    # append array for each year to list
    arrays_list[[y]] <- array_for_year
    year_month_df <- rbind(year_month_df, data.frame(year = rep(y), month = seq(1, 12)))
  }
  all_years_co <- do.call(abind, (c(arrays_list, along = 2))) # # grid cells vs total num months
  if (city == "shg") {
    scale <- scale_edgar_by_htap(urban_domain_rast)
    print("scaled by HTAP")
    all_years_co_scaled <- all_years_co * scale
    all_years_co <- all_years_co_scaled
  }
  # find mean co emissions for time period of interest (aggregate)
  co_emissions_combined <- combine_monthly_emissions(all_years_co, year_month_df, obs_year_month_list, aggregate, all_months)
  if (uncertainty_method == "full_uncertainty") {
    # add uncertainty based on EDGAR v4.3.2 SI
    sigma_value <- find_EDGAR_uncertainty(city) # from urban_boundaries.R
    # sigma from find_EDGAR_uncertainty.R is a fraction, must be multiplied by mean value
  } else if (uncertainty_method == "SRON_uncertainty") {
    sigma_value <- 0.1210457 # SRON CO CSF study, average CO change
  } else if (uncertainty_method == "no_uncertainty") {
    sigma_value <- 0
  }
  co_emissions_realizations <- data.frame()
  for (time_step in seq_len(dim(co_emissions_combined)[1])) {
      samples_list <- generate_random_samples(co_emissions_combined[time_step, "domain_mean"], n_samples = 1000, sigma = sigma_value * co_emissions_combined[time_step, "domain_mean"])
      co_emissions_realizations <- rbind(co_emissions_realizations,
                                            data.frame(co_emissions_combined[time_step, ] %>% select(-c("domain_mean")),
                                                       t(samples_list)))
  }

  return(co_emissions_realizations) # first column is year or season year, rest are 1000 realizations for the sum (annual,seasonal, etc.)
}