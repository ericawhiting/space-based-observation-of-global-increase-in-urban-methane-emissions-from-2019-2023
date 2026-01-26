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
source(paste0(path_to_utils, "load_co_emissions.R")) # return_cell_values and combine_monthly_emissions functions
# ------------------------------------------------------------------------------------------------------------------

open_edgar_ch4_monthly_file <- function(edgar_y, sector, version) { # nolint: object_name_linter.
  #' open netcdf4 files for EDGARv8.0 CH4 given year, sector, and return opened variable
  #' @param edgar_y is year to pull for EDGAR file names (EDGARv8.0 ends in 2022)
  #' @param sector EDGARv8.0 CH4 sectors: "AGRICULTURE", "BUILDINGS", "FUEL_EXPLOITATION", "IND_COMBUSTION", "POWER_INDUSTRY", "TRANSPORT", "WASTE"
  #' @return returns an array of emissions (rate not a flux!) of CO emissions in Mg/month from EDGAR for month in year for sector for whole globe
  # --------
  if (version == "v8") {
    fname <- paste0(path_to_data, "EDGAR/ch4/monthly_sector_specific_gridmaps/v8/", "v8.0_FT2022_GHG_CH4_", edgar_y, "_", sector, "_emi.nc")
  } else if (version == "2024") {
    fname <- paste0(path_to_data, "EDGAR/ch4/monthly_sector_specific_gridmaps/2024/", "EDGAR_2024_GHG_CH4_", edgar_y, "_bkl_", sector, "_emi.nc")
  } else {
    print("edgar ch4 monthly file not found")
  }
  nc_file <- ncdf4::nc_open(fname)
  var_name <- names(nc_file$var)
  var <- ncdf4::ncvar_get(nc_file, var_name)
  ncdf4::nc_close(nc_file)
  return(var)
}


load_all_month_sectors_ch4 <- function(y, urban_domain_rast, version) {
  #' load and sum all edgar ch4 values across 8 sectors for each month for each grid cell, combine into dataframe for year
  #' @param y year to pull edgar file (v8.0 goes until 2022)
  #' @param urban_domain_rast the terra extent of the spatial domain
  #' @return values_of_urban_rast monthly emissions in Tonnes or Mg per month of dimension # gridcells by # months in year of file
  # --------
  if (version == "v8") {
    # accounting for changes made to edgar v8 ch4 files
    var_1 <- open_edgar_ch4_monthly_file(y, "bkl_AGRICULTURE", version)
  } else {
    var_1 <- open_edgar_ch4_monthly_file(y, "AGRICULTURE", version)
  }
  #var_1 <- open_edgar_ch4_monthly_file(y, "bkl_AGRICULTURE", version)
  var_2 <- open_edgar_ch4_monthly_file(y, "BUILDINGS", version)
  var_3 <- open_edgar_ch4_monthly_file(y, "FUEL_EXPLOITATION", version)
  var_4 <- open_edgar_ch4_monthly_file(y, "IND_COMBUSTION", version)
  var_5 <- open_edgar_ch4_monthly_file(y, "IND_PROCESSES", version)
  var_6 <- open_edgar_ch4_monthly_file(y, "POWER_INDUSTRY", version)
  var_7 <- open_edgar_ch4_monthly_file(y, "TRANSPORT", version)
  var_8 <- open_edgar_ch4_monthly_file(y, "WASTE", version)

  monthly_sum_across_sectors <- var_1 + var_2 + var_3 + var_4 + var_5 + var_6 + var_7 + var_8
  month_list <- seq(1, 12)
  values_of_urban_rast <- c()
  for (m in month_list) { # gives array where columns are each month and rows are values of grid cells
      values_of_urban_rast <-  cbind(values_of_urban_rast, return_cell_values(monthly_sum_across_sectors, urban_domain_rast, m))
  }
  return(values_of_urban_rast) # monthly emissions are Mg / month, Tonnes / month
}

load_all_months_sectors_individually_ch4 <- function(city, y, urban_domain_rast, version) {
  if (version == "v8") {
      # accounting for changes made to edgar v8 ch4 files
      var_1 <- open_edgar_ch4_monthly_file(y, "bkl_AGRICULTURE", version)
  } else {
      var_1 <- open_edgar_ch4_monthly_file(y, "AGRICULTURE", version)
  }
  var_2 <- open_edgar_ch4_monthly_file(y, "BUILDINGS", version)
  var_3 <- open_edgar_ch4_monthly_file(y, "FUEL_EXPLOITATION", version)
  var_4 <- open_edgar_ch4_monthly_file(y, "IND_COMBUSTION", version)
  var_5 <- open_edgar_ch4_monthly_file(y, "IND_PROCESSES", version)
  var_6 <- open_edgar_ch4_monthly_file(y, "POWER_INDUSTRY", version)
  var_7 <- open_edgar_ch4_monthly_file(y, "TRANSPORT", version)
  var_8 <- open_edgar_ch4_monthly_file(y, "WASTE", version)

  month_list <- seq(1, 12)
  values_of_urban_rast <- data.frame()
  for (m in month_list) { # gives array where columns are each month and rows are values of grid cells
      values_of_urban_rast <-  rbind(values_of_urban_rast,
                                     data.frame(city = city, year = y, month = m, version = version,
                                        AGRICULTURE = colSums(return_cell_values(var_1, urban_domain_rast, m)),
                                        BUILDINGS = colSums(return_cell_values(var_2, urban_domain_rast, m)),
                                        FUEL_EXPLOITATION = colSums(return_cell_values(var_3, urban_domain_rast, m)),
                                        IND_COMBUSTION = colSums(return_cell_values(var_4, urban_domain_rast, m)),
                                        IND_PROCESSES = colSums(return_cell_values(var_5, urban_domain_rast, m)),
                                        POWER_INDUSTRY = colSums(return_cell_values(var_6, urban_domain_rast, m)),
                                        TRANSPORT = colSums(return_cell_values(var_7, urban_domain_rast, m)),
                                        WASTE = colSums(return_cell_values(var_8, urban_domain_rast, m))))
  }
  row.names(values_of_urban_rast) <- NULL
  return(values_of_urban_rast) # monthly emissions are Mg / month, Tonnes / month
}

find_ch4_emissions_from_monthly <- function(city, urban_domain_rast, year_list, obs_year_month_list, version, aggregate, all_months) {
  #' uses edgarv8.0 monthly ch4 files to list out numerical realizations of mean ch4 emissions for given time aggregation level for given domain
  #' @param city city code from which to search for EDGAR uncertainty value
  #' @param urban_domain_rast the terra extent of the spatial domain
  #' @param year_list a list of years to evaluate edgar over
  #' @param obs_year_month_list df of unique months of observations with columns year, month, season_year
  #' @param version v8 or 2024 version of EDGAR ch4
  #' @param aggregate annual, seasonal, monthly: the time step over which we aggregate and find mean emissions
  #' @param all_months TRUE or FALSE, if TRUE then looks at all months of given time step aggregation, if false, uses months in which we had obs
  #' @return ch4_emissions_realizations df of relevant time steps and numerical realizations of ch4 emissions for that time period
  # --------
  arrays_list <- list()
  year_month_df <- data.frame(year = numeric(), month = numeric())
  for (y in year_list) {
    if ((version == "v8") && (y == 2023 || y == 2024)) {
      edgar_y <- 2022 # run 2022 edgar for 2023 and 2024
    } else if ((version == "2024") && (y == 2024)) {
      edgar_y <- 2023 # run 2023 edgar for 2024
    } else {
      edgar_y <- y
    }
    # 2d array with dim: #grid cell, #month in year
    array_for_year <- load_all_month_sectors_ch4(edgar_y, urban_domain_rast, version)
    # append array for each year to list
    arrays_list[[y]] <- array_for_year
    year_month_df <- rbind(year_month_df, data.frame(year = rep(y), month = seq(1, 12)))
  }
  all_years_ch4 <- do.call(abind, (c(arrays_list, along = 2))) # # grid cells vs total num months
  # find mean ch4 emissions for time period of interest (aggregate)
  ch4_emissions_combined <- combine_monthly_emissions(all_years_ch4, year_month_df, obs_year_month_list, aggregate, all_months)
  return(ch4_emissions_combined) # first column is year or season year, next is aggregated mean
}
