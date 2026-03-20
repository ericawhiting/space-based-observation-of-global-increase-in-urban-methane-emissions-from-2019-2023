# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(dplyr)
library(foreach)

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "load_ch4_emissions.R"))

path_to_data <- "../Data/" # change your path here

ch4_emissions_path <- paste0(path_to_data, "EDGAR_ch4_emissions/")
# may be averaged by month, year, etc

# set EDGAR CH4 version:
version <- "v8" # or "2024"
if (version == "v8") {
  version_name <- "v8.0"
} else if (version == "2024") {
  version_name <- "v2024"
} else {
  print("unknown EDGAR ch4 version")
}

# -------------------------------
# Load City List
# -------------------------------
city_list <- read.csv(paste0(path_to_data, "city_information/city_codes.csv"), header = FALSE)$V1
# -------------------------------
# Iterate Over Cities
# -------------------------------
foreach(i = seq_along(city_list)) %do% {
    # Load City
    city <- city_list[i]
    # Check if there are observations, if yes, then save out CH4 emissions
    ch4co_df <- read.csv(paste0(path_to_data, "ch4co/", city, "_ch4co.csv"))
    if (dim(ch4co_df)[1] > 0) {
        # Load Urban Domain
        city_latlon <- find_city_center_from_pop_dens(city)
        delta_lat <- read_ideal_boxsize_lat(city)
        lat_adj <- 1 / cos(city_latlon[1] * (pi / 180.0)) # used to adjust delta longitudes for boxes based on where city is latitude wise
        delta_lon <- delta_lat * lat_adj
        coastal_flag <- find_city_coastal_flag(city) # return 1 if near coast or considered coastal city by eye
        delta_lat <- delta_lat + 0.25 * delta_lat * coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land
        delta_lon <- delta_lon + 0.25 * delta_lon * coastal_flag
        urban_domain_rast <- rast(terra::ext(as.vector(cbind(city_latlon[2] - delta_lon,
                                                             city_latlon[2] + delta_lon,
                                                             city_latlon[1] - delta_lat,
                                                             city_latlon[1] + delta_lat))))

        # Create List of Month, Year of Observations
        obs_year_month_list <- ch4co_df %>% select(c("year", "month")) %>% group_by(year, month) %>% unique()
        obs_year_month_list$season_year <- mapply(assign_season_year, obs_year_month_list$year, obs_year_month_list$month)
        year_list <- seq(2019, 2024)

        # -------------------------------
        # ANNUAL
        # -------------------------------
        annual_ch4_emissions <- find_ch4_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, version, aggregate = "annual", all_months = TRUE)
        write.table(annual_ch4_emissions, paste0(ch4_emissions_path, "annual/", city, "_annual_ch4_emissions_from_monthly_edgar", version_name, ".csv"), row.names = FALSE, sep = ",")

        # Months that we have Obs in Only
        # annual_ch4_emissions_matchmos <- find_ch4_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, version, aggregate = "annual", all_months = FALSE)
        # write.table(annual_ch4_emissions_matchmos, paste0(ch4_emissions_path, "annual/", city, "_annual_ch4_emissions_matchmos_from_monthly_edgar", version_name, ".csv"), row.names = FALSE, sep = ",")

        # -------------------------------
        # MONTHLY
        # -------------------------------
        #monthly_ch4_emissions <- find_ch4_emissions_from_monthly(city, urban_domain_rast, year_list, obs_year_month_list, version, aggregate = "monthly", all_months = TRUE)
        #write.table(monthly_ch4_emissions, paste0(ch4_emissions_path, "monthly/", city, "_monthly_ch4_emissions_from_monthly_edgar", version_name, ".csv"), row.names = FALSE, sep = ",")
    }
    print(paste0(city, " done with aggregating ch4 emissions"))
}
