# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(dplyr)
library(foreach)

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "load_co_emissions.R"))

path_to_data <- "../Data/" # change your path here

co_emissions_path_full_uncertainty <- paste0(path_to_data, "EDGAR_co_emissions_realizations/full_uncertainty/")
co_emissions_path_SRON_uncertainty <- paste0(path_to_data, "EDGAR_co_emissions_realizations/SRON_uncertainty/")
co_emissions_path_no_uncertainty  <-  paste0(path_to_data, "EDGAR_co_emissions_realizations/no_uncertainty/")

# -------------------------------
# Load City List
# -------------------------------
city_list <- read.csv(paste0(path_to_data, "city_information/city_codes.csv"), header = FALSE)$V1

# -------------------------------
# Iterate Over Cities
# -------------------------------
foreach(i = seq_along(city_list)) %dopar% {
    # Load City
    city <- city_list[i]
    ch4co_df <- read.csv(paste0(path_to_data, "ch4co/", city, "_ch4co.csv"))

    if (dim(ch4co_df)[1] > 0) { # only runs for cities there are obs for, could change this
        sf_flag <- FALSE
        # Load Urban Domain
        city_latlon <- find_city_center_from_pop_dens(city)
        delta_lat <- read_ideal_boxsize_lat(city)
        lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0)) # used to adjust delta longitudes for boxes based on where city is latitude wise
        delta_lon <- delta_lat * lat_adj
        coastal_flag <- find_city_coastal_flag(city) # return 1 if near coast or considered coastal city by eye
        delta_lat <- delta_lat + 0.25 * delta_lat * coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land
        delta_lon <- delta_lon + 0.25 * delta_lon * coastal_flag
        urban_domain_rast <- rast(terra::ext(as.vector(cbind(city_latlon[2] - delta_lon, city_latlon[2] + delta_lon, city_latlon[1] - delta_lat, city_latlon[1] + delta_lat))))

        # Create List of Month, Year of Observations
        obs_year_month_list <- ch4co_df %>% select(c("year", "month")) %>% group_by(year, month) %>% unique()
        obs_year_month_list$season_year <- mapply(assign_season_year, obs_year_month_list$year, obs_year_month_list$month)
        year_list <- seq(2019, 2024)


        # -------------------------------
        # ANNUAL
        # uncomment which one you want to run
        # -------------------------------
        # No Uncertainty, All Months
        # annual_co_emissions_no_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "annual", all_months = TRUE, uncertainty_method = "no_uncertainty", sf_flag = sf_flag)
        # write.table(annual_co_emissions_no_uncertainty, paste0(co_emissions_path_no_uncertainty, "annual/", city, "_realizations_annual_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",")

        # # Full Uncertainty, All Months
        # annual_co_emissions_full_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "annual", all_months = TRUE, uncertainty_method = "full_uncertainty", sf_flag = sf_flag)
        # write.table(annual_co_emissions_full_uncertainty, paste0(co_emissions_path_full_uncertainty, "annual/", city, "_realizations_annual_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",")
        # SRON Uncertainty, All Months
        annual_co_emissions_SRON_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "annual", all_months = TRUE, uncertainty_method = "SRON_uncertainty", sf_flag = sf_flag)
        write.table(annual_co_emissions_SRON_uncertainty, paste0(co_emissions_path_SRON_uncertainty, "annual/", city, "_realizations_annual_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",")
        # Full Uncertainty, Months with Obs Only
        # annual_co_emissions_full_uncertainty_matchmos <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "annual", all_months = FALSE, uncertainty_method = "full_uncertainty", sf_flag = sf_flag)
        # write.table(annual_co_emissions_full_uncertainty_matchmos, paste0(co_emissions_path_full_uncertainty, "annual/", city, "_realizations_annual_co_emissions_matchmos_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",")
        # # SRON Uncertainty, Months with Obs Only
        # annual_co_emissions_SRON_uncertainty_matchmos <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "annual", all_months = FALSE, uncertainty_method = "SRON_uncertainty", sf_flag = sf_flag)
        # write.table(annual_co_emissions_SRON_uncertainty_matchmos, paste0(co_emissions_path_SRON_uncertainty, "annual/", city, "_realizations_annual_co_emissions_matchmos_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",")

        # -------------------------------
        # SEASONAL
        # -------------------------------
        # seasonal_co_emissions_full_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "seasonal", all_months = TRUE, uncertainty_method = "full_uncertainty") # nolint: line_length_linter, object_length_linter, object_name_linter.
        # write.table(seasonal_co_emissions_full_uncertainty, paste0(co_emissions_path_full_uncertainty, "seasonal/", city, "_realizations_seasonal_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",") # nolint: line_length_linter, object_length_linter, object_name_linter.

        # seasonal_co_emissions_full_uncertainty_matchmos <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "seasonal", all_months = FALSE, uncertainty_method = "full_uncertainty") # nolint: line_length_linter, object_length_linter, object_name_linter.
        # write.table(seasonal_co_emissions_full_uncertainty_matchmos, paste0(co_emissions_path_full_uncertainty, "seasonal/", city, "_realizations_seasonal_co_emissions_matchmos_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",") # nolint: line_length_linter, object_length_linter, object_name_linter.

        # -------------------------------
        # MONTHLY
        # -------------------------------
        # monthly_co_emissions_full_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "monthly", all_months = TRUE, uncertainty_method = "full_uncertainty") # nolint: line_length_linter, object_length_linter, object_name_linter.
        # write.table(monthly_co_emissions_full_uncertainty, paste0(co_emissions_path_full_uncertainty, "monthly/", city, "_realizations_monthly_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",") # nolint: line_length_linter, object_length_linter, object_name_linter.
        # monthly_co_emissions_no_uncertainty <- find_co_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, aggregate = "monthly", all_months = TRUE, uncertainty_method = "no_uncertainty") # nolint: line_length_linter, object_length_linter, object_name_linter.
        # write.table(monthly_co_emissions_no_uncertainty, paste0(co_emissions_path_no_uncertainty, "monthly/", city, "_realizations_monthly_co_emissions_from_monthly_edgarv8.1.csv"), row.names = FALSE, sep = ",") # nolint: line_length_linter, object_length_linter, object_name_linter.
    }
     print(paste0(city, " done with aggregating co emissions"))
}