# Load Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(foreach)
  library(doParallel)
  library(viridis)
})

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))
source(paste0(path_to_utils, "split_dates.R")) # assign_season_year
source(paste0(path_to_utils, "CSF_analysis.R"))

path_to_data <- "../Data/" # change your path here

# ---------------------------------------------------------------------------------------------
tg_y_conversion <- 3600 * 24 * 365.25 * 1e-9 # from kg/s

# ---------------------------------------------------------------------------------------------
# Load CH4CO Data
# ---------------------------------------------------------------------------------------------
ch4co_data_path <- paste0(path_to_data, "ch4co/") # filename format: paste0(city, "_ch4co.csv")

# Read in CH4 CO Data for C40 cities
no_pass_list <- c() # for cities that had 0 overpasses
ch4co_df <- data.frame() # initialize df for combining all C40 city information
# iterate through all csv data from city analysis, append to df, add city code, name, and region
for (file in list.files(ch4co_data_path, pattern = "*.csv")) {
  city <- str_split(file, "_")[[1]][1]
  temp_ch4co_df <- read.csv(paste0(ch4co_data_path, file), header = TRUE)
  if (dim(temp_ch4co_df)[1] > 0) {
    temp_ch4co_df$city <- city
  } else {
    no_pass_list <- append(no_pass_list, city)
  }
  ch4co_df <- rbind(ch4co_df, temp_ch4co_df)
}

# -------------------------------
# Filter CH4CO Data
# -------------------------------
# Remove 2024
ch4co_df <- ch4co_df %>% filter(year < 2024)
# Metropolitan Areas / dealing with double counting
ch4co_df <- ch4co_df[(ch4co_df$city != "tok") & (ch4co_df$city != "eku") &
                     (ch4co_df$city != "tsh") & (ch4co_df$city != "hk") &
                     (ch4co_df$city != "rot") & (ch4co_df$city != "zhe"), ]

# FILTERING (cities must have >=5 overpasses per year)
removed_cities_5 <- ch4co_df %>% group_by(city, year) %>% summarise(num_overpasses = sum(num_overpasses)) %>% filter(num_overpasses < 5)
ch4co_df <- anti_join(ch4co_df, removed_cities_5, by = c("city", "year"))
# Filter out negative ch4 enh
ch4co_df <- ch4co_df %>% filter(summed_ch4_enhancements_15 > 0) # ensure all ch4 enhancements are positive

# check if has data for 2023, necessary for all city total
removed_cities_no2023 <- ch4co_df %>%
                          group_by(city) %>%
                          summarize(has_2023_data = any(year == 2023)) %>%
                          filter(!has_2023_data) %>%
                          select(city)
ch4co_df <- anti_join(ch4co_df, removed_cities_no2023, by = c("city"))

# remove C40 cities that do not have data across all years because that is where nonC40 cities are used here
unique_years <- seq(2019, 2023) #unique(ch4co_df$year)
ch4co_df$C40 <- apply(ch4co_df[c("city")], 1, function(x) is_C40(x))
removed_cities_nonC40 <- ch4co_df %>%
                          filter(C40 == 0) %>%
                          group_by(city, year) %>%
                          summarize(data_count = n()) %>%
                          group_by(city) %>%
                          summarize(years_with_data = n_distinct(year)) %>%
                          filter(years_with_data < length(unique_years)) %>%
                          select(city)
ch4co_df <- anti_join(ch4co_df, removed_cities_nonC40, by = c("city"))


# Add Datetime information
ch4co_df$yrmo <- format(as.Date(paste(ch4co_df$year, ch4co_df$month, "01", sep = "-")))
ch4co_df$season_year <- mapply(assign_season_year, ch4co_df$year, ch4co_df$month)
# Create list of cities to iterate through
cities <- unique(ch4co_df$city)

# ---------------------------------------------------------------------------------------------
# Produce Realizations of CH4CO ratio
# ---------------------------------------------------------------------------------------------
find_ch4co_realizations <- function(ch4_data, co_data, prec_data, n_data) {
    length_to_sample <- dim(ch4_data)[1]
    bkg_option <- seq(10, 20, 1)
    df_resample_ch4co <- data.frame(NULL)
    # create realizations of ratio using uncertainty on ratio
    for (i in seq_len(1000)) {
        ch4co_data_temp <- data.frame(NULL)
        # sample dates and backgrounds
        date_sample <- sample(seq_len(length_to_sample), length_to_sample, replace = TRUE)
        bkg_sample <- sample(bkg_option, length_to_sample, replace = TRUE)

        for (j in seq_len(length_to_sample)) {
            bkg <- bkg_sample[j]
            k <- date_sample[j]

            bkg_name_ch4 <- paste0("summed_ch4_enhancements_", bkg)
            bkg_name_co <- paste0("summed_co_enhancements_", bkg)

            co <- as.numeric(co_data[k, paste0(bkg_name_co)])
            ch4 <- as.numeric(ch4_data[k, paste0(bkg_name_ch4)])
            sd_co <- as.numeric(prec_data[k, "mean_co_prec"])
            sd_ch4 <- as.numeric(prec_data[k, "mean_ch4_prec"]) * 2 # based on comparison to TCCON ground sites (SRON 2020)
            cov_ch4co <- 0
            if ((co > (n_data[k, "n_points_vec"] * 7)) && (ch4 > 0)) {
                uncert <- ch4 / co * sqrt((sd_ch4 / ch4)^2 + (sd_co / co)^2 - cov_ch4co / (ch4 * co))
                noisy_ch4co <- rnorm(1, mean = ch4 / co, sd = uncert)
                ch4co_data_temp <- rbind(ch4co_data_temp, noisy_ch4co)
            } else {
                ch4co_data_temp <- rbind(ch4co_data_temp, NA)
            }
        }
        df_resample_ch4co <- rbind(df_resample_ch4co, t(ch4co_data_temp))
    }
    # return annual mean for each annual resampling
    return(data.frame(ch4co = as.numeric(rowMeans(df_resample_ch4co, na.rm = TRUE))))
    # alternatively could resample by randomly select one ratio from each column 1000 times
}

# ---------------
# Prep CH4CO Dataframes
ch4_enh_df <- ch4co_df %>% select(c(datetimes, year, month, city, season_year,
                          summed_ch4_enhancements_10, summed_ch4_enhancements_11, summed_ch4_enhancements_12, summed_ch4_enhancements_13,
                          summed_ch4_enhancements_14, summed_ch4_enhancements_15, summed_ch4_enhancements_16, summed_ch4_enhancements_17,
                          summed_ch4_enhancements_18, summed_ch4_enhancements_19, summed_ch4_enhancements_20))
co_enh_df <- ch4co_df %>% select(c(datetimes, year, month,  city, season_year,
                          summed_co_enhancements_10, summed_co_enhancements_11, summed_co_enhancements_12, summed_co_enhancements_13,
                          summed_co_enhancements_14, summed_co_enhancements_15, summed_co_enhancements_16, summed_co_enhancements_17,
                          summed_co_enhancements_18, summed_co_enhancements_19, summed_co_enhancements_20))
ch4co_enh_df <- ch4co_df %>% select(c(datetimes, year, month, city, season_year, integ_results_10, integ_results_11, integ_results_12, integ_results_13,
                            integ_results_14, integ_results_15, integ_results_16, integ_results_17, integ_results_18, integ_results_19,
                            integ_results_20))
prec_df <- ch4co_df %>% select(c(datetimes, year, month, city, season_year,  mean_ch4_prec, mean_co_prec))
n_point_df <- ch4co_df %>% select(c(datetimes, year, month, city, season_year, n_points_vec))

# ---------------
# Iterate over cities and find CH4CO realizations

# remove parallelization call and change %dopar% to %do% if needed
num_cores <-  30
registerDoParallel(30)

# remove parallelization call and change %dopar% to %do% if needed
annual_ch4co_realizations <- foreach(c = seq_along(cities), .combine = rbind) %dopar% {
    # Assign City
    city <- cities[c]
    # Assign Empty dataframe
    ch4co_realizations <- data.frame(NULL)
    # Iterate over years
    for (y in unique(ch4_enh_df[ch4_enh_df$city == city, ]$year)) {
        ch4_data <- ch4_enh_df[(ch4_enh_df$city == city) & (ch4_enh_df$year == y), ]
        co_data <- co_enh_df[(co_enh_df$city == city) & (co_enh_df$year == y), ]
        prec_data <- prec_df[(prec_df$city == city) & (prec_df$year == y), ]
        n_data <- n_point_df[(n_point_df$city == city) & (n_point_df$year == y), ]

        length_to_sample <- dim(ch4_data)[1] # how many overpasses are in one year
        bstrap <- find_ch4co_realizations(ch4_data, co_data, prec_data, n_data)
        ch4co_realizations <- rbind(ch4co_realizations,
                                    data.frame(city = rep(city, 3), y = rep(y, 3), length_to_sample = rep(length_to_sample, 3),
                                               variable = c("ch4co"), t(bstrap)))
        print(paste0("done with: ", city, " ", y))
    }
    return(ch4co_realizations)
}

# ---------------------------------------------------------------------------------------------
# Produce Realizations of Annual CO Emissions from EDGARv8.1
# ---------------------------------------------------------------------------------------------
co_emissions_path_full_uncertainty <- paste0(path_to_data, "EDGAR_co_emissions_realizations/revised/full_uncertainty/")
co_emissions_path_SRON_uncertainty <- paste0(path_to_data, "EDGAR_co_emissions_realizations/revised/SRON_uncertainty/")

# ANNUAL FULL UNCERT
co_full_realizations <- foreach(c = seq_along(cities), .combine = rbind) %do% {
    city <- cities[c]
    annual_co_emissions <- cbind(data.frame(city = city),
                                 read.csv(paste0(co_emissions_path_full_uncertainty, "annual/", city, "_realizations_annual_co_emissions_from_monthly_edgarv8.1.csv")))
    return(annual_co_emissions)
}

# ANNUAL SRON UNCERT
co_SRON_realizations <- foreach(c = seq_along(cities), .combine = rbind) %do% {
    city <- cities[c]
    annual_co_emissions <- cbind(data.frame(city = city),
                                 read.csv(paste0(co_emissions_path_SRON_uncertainty, "annual/", city, "_realizations_annual_co_emissions_from_monthly_edgarv8.1.csv")))
    return(annual_co_emissions)
}

# ANNUAL FULL UNCERT MATCHMOS
co_full_realizations_matchmos <- foreach(c = seq_along(cities), .combine = rbind) %do% {
    city <- cities[c]
    annual_co_emissions <- cbind(data.frame(city = city),
                                 read.csv(paste0(co_emissions_path_full_uncertainty, "annual/", city, "_realizations_annual_co_emissions_matchmos_from_monthly_edgarv8.1.csv")))
    return(annual_co_emissions)
}

# ANNUAL SRON UNCERT MATCHMOS
co_SRON_realizations_matchmos <- foreach(c = seq_along(cities), .combine = rbind) %do% {
    city <- cities[c]
    annual_co_emissions <- cbind(data.frame(city = city),
                                 read.csv(paste0(co_emissions_path_SRON_uncertainty, "annual/", city, "_realizations_annual_co_emissions_matchmos_from_monthly_edgarv8.1.csv")))
    return(annual_co_emissions)
}

# ---------------------------------------------------------------------------------------------
# CH4 Emissions from Realizations
# check length here
# ---------------------------------------------------------------------------------------------
combine_ch4co_and_co_annual_realizations <- function(ch4co_annual_realizations, co_annual_realizations) {
  mol_frac <- (16.04 / 28.01) # molCH4 / molCO * gCH4 / molCH4
  annual_ch4 <- data.frame(NULL)
  for (city in unique(ch4co_annual_realizations$city)) {
      region <- find_city_region(city)
      for (year in unique(ch4co_annual_realizations[ch4co_annual_realizations$city == city, "y"])) {
          ch4co_array <- ch4co_annual_realizations[(ch4co_annual_realizations$city == city) & (ch4co_annual_realizations$y == year), c(4:1003)]
          num_obs <- ch4co_annual_realizations[(ch4co_annual_realizations$city == city) & (ch4co_annual_realizations$y == year), "length_to_sample"]
          co_array <- co_annual_realizations[(co_annual_realizations$city == city) & (co_annual_realizations$year == year), c(3:1002)] # first column is city, second is year
          if (dim(ch4co_array)[1] > 0) {
              bstrap_data <- as.numeric(ch4co_array[]) * mol_frac * as.numeric(co_array[]) # kg/s # both ch4co and co are already randomly noisy, do not need to re-randomize
              annual_ch4 <- rbind(annual_ch4, data.frame(city, year, region, num_obs, t(bstrap_data)))
          } else {
              annual_ch4 <- rbind(annual_ch4, data.frame(city, year, region, num_obs = NA, t(array(NA, 1000))))
          }
      }
  }
  return(annual_ch4)
}

# find ch4 realizations
annual_ch4_full <- combine_ch4co_and_co_annual_realizations(annual_ch4co_realizations, co_full_realizations)
annual_ch4_full_matchmos <- combine_ch4co_and_co_annual_realizations(annual_ch4co_realizations, co_full_realizations_matchmos)
annual_ch4_SRON <- combine_ch4co_and_co_annual_realizations(annual_ch4co_realizations, co_SRON_realizations)
annual_ch4_SRON_matchmos <- combine_ch4co_and_co_annual_realizations(annual_ch4co_realizations, co_SRON_realizations_matchmos)

# pare down CO files to match obs lengths:
annual_co_SRON <- co_SRON_realizations %>%
                            semi_join(annual_ch4_full, by = c("city", "year"))

annual_co_SRON_matchmos <- co_SRON_realizations_matchmos %>%
                            semi_join(annual_ch4_full, by = c("city", "year"))

# -------------------
# City Ratio and Emissions (Detection of Change) (emissions as % change relative to 2019)
# for each city, pull mean ch4co, CI, mean emissions detection over change and CI
# -------------------
city_annual_df <- data.frame(NULL)
for (ii in seq_len(dim(annual_ch4_full)[1])) {
  city_annual_df <- rbind(city_annual_df,
                          data.frame(annual_ch4_full[ii, "city"], # city data
                                     annual_ch4_full[ii, "year"],
                                     # annual_ch4[ii, "m"], #num months
                                     annual_ch4_full[ii, "region"],
                                     annual_ch4_full[ii, "num_obs"],
                                     # ch4co data
                                     mean(as.numeric(annual_ch4co_realizations[ii, c(4:1003)]), na.rm = TRUE),
                                     quantile(as.numeric(annual_ch4co_realizations[ii, c(4:1003)]), 0.025, na.rm = TRUE),
                                     quantile(as.numeric(annual_ch4co_realizations[ii, c(4:1003)]), 0.975, na.rm = TRUE),
                                     # ch4 data
                                     mean(as.numeric(annual_ch4_full[ii, c(4:1003)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_full[ii, c(4:1003)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_full[ii, c(4:1003)]), 0.975, na.rm = TRUE) * tg_y_conversion,
                                     # ch4 matchmos
                                     mean(as.numeric(annual_ch4_full_matchmos[ii, c(4:1003)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_full_matchmos[ii, c(4:1003)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_full_matchmos[ii, c(4:1003)]), 0.975, na.rm = TRUE) * tg_y_conversion,
                                     # ch4 data sron uncert
                                     mean(as.numeric(annual_ch4_SRON[ii, c(4:1003)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_SRON[ii, c(4:1003)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_SRON[ii, c(4:1003)]), 0.975, na.rm = TRUE) * tg_y_conversion,
                                     # ch4 matchmos sron uncert
                                     mean(as.numeric(annual_ch4_SRON_matchmos[ii, c(4:1003)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_SRON_matchmos[ii, c(4:1003)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_ch4_SRON_matchmos[ii, c(4:1003)]), 0.975, na.rm = TRUE) * tg_y_conversion,
                                     # co sron
                                     mean(as.numeric(annual_co_SRON[ii, c(3:1002)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_co_SRON[ii, c(3:1002)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_co_SRON[ii, c(3:1002)]), 0.975, na.rm = TRUE) * tg_y_conversion,
                                     # co sron matchmos
                                     mean(as.numeric(annual_co_SRON_matchmos[ii, c(3:1002)]), na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_co_SRON_matchmos[ii, c(3:1002)]), 0.025, na.rm = TRUE) * tg_y_conversion,
                                     quantile(as.numeric(annual_co_SRON_matchmos[ii, c(3:1002)]), 0.975, na.rm = TRUE) * tg_y_conversion
                                     ))
}
rownames(city_annual_df) <- NULL
colnames(city_annual_df) <- c("city", "year", "region", "num_obs", "ch4co", "ch4co_lower", "ch4co_upper",
                              "ch4_full", "ch4_full_lower", "ch4_full_upper", "ch4_full_matchmos", "ch4_full_lower_matchmos", "ch4_full_upper_matchmos",
                              "ch4_sron", "ch4_sron_lower", "ch4_sron_upper", "ch4_sron_matchmos", "ch4_sron_lower_matchmos", "ch4_sron_upper_matchmos",
                              "co_sron", "co_sron_lower", "co_sron_upper", "co_sron_matchmos", "co_sron_lower_matchmos", "co_sron_upper_matchmos")

# Emissions in Tg/y
write.table(city_annual_df, paste0(path_to_data, "results/annual_ch4_emissions.csv"), row.names = FALSE, sep = ",")




# -------------------
# Regional Emissions Sums
# -------------------
annual_ch4_SRON$C40 <- apply(annual_ch4_SRON[c("city")], 1, function(x) is_C40(x))

annual_ch4_C40_SRON_allyears <- annual_ch4_SRON %>%
                                     filter(year != 2024) %>%
                                     group_by(region, city) %>%
                                     filter(length(unique(year)) >= length(seq(2019, 2023))) %>%
                                     ungroup() %>%
                                     filter(C40 == 1)
regional_realizations <- data.frame(NULL)
annual_regional_sum_df <- data.frame(NULL)
for (region in unique(annual_ch4_C40_SRON_allyears$region)) {
  for (year in c(2019, 2020, 2021, 2022, 2023)) {
    bstrap <- c()
    ch4_year <- as.data.frame(annual_ch4_C40_SRON_allyears[(annual_ch4_C40_SRON_allyears$year == year) & (annual_ch4_C40_SRON_allyears$region == region), c(5:1004)])
    bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE)) # sum in each column bc already randomized
    regional_realizations <- rbind(regional_realizations,
                                  data.frame(year, region, n = dim(ch4_year)[1], t(bstrap * tg_y_conversion)))
    annual_regional_sum_df <- rbind(annual_regional_sum_df,
                                data.frame(year, region, dim(ch4_year)[1],
                                           mean(bstrap, na.rm = TRUE), sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
  }
}
rownames(annual_regional_sum_df) <- NULL
colnames(annual_regional_sum_df) <- c("year", "region", "n", "mean", "sd", "lower", "upper")
annual_regional_sum_df <- annual_regional_sum_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_regional_sum_df
# Emissions in Tg/y


# -------------------
# Regional Emissions Sums MATCHMONTHS
# -------------------
annual_ch4_SRON_matchmos$C40 <- apply(annual_ch4_SRON_matchmos[c("city")], 1, function(x) is_C40(x))

annual_ch4_C40_SRON_matchmos_allyears <- annual_ch4_SRON_matchmos %>%
                                     filter(year != 2024) %>%
                                     group_by(region, city) %>%
                                     filter(length(unique(year)) >= length(seq(2019, 2023))) %>%
                                     # filter(length(unique(year)) == length(unique(annual_ch4_SRON$year))) %>%
                                     ungroup() %>%
                                     filter(C40 == 1)

annual_regional_sum_matchmos_df <- data.frame(NULL)
for (region in unique(annual_ch4_C40_SRON_matchmos_allyears$region)) {
  for (year in c(2019, 2020, 2021, 2022, 2023)) {
    bstrap <- c()
    ch4_year <- as.data.frame(annual_ch4_C40_SRON_matchmos_allyears[(annual_ch4_C40_SRON_matchmos_allyears$year == year) & (annual_ch4_C40_SRON_matchmos_allyears$region == region), c(5:1004)])
    bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE)) # sum in each column bc already randomized

    annual_regional_sum_matchmos_df <- rbind(annual_regional_sum_matchmos_df,
                                data.frame(year, region, dim(ch4_year)[1],
                                           mean(bstrap, na.rm = TRUE), sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
  }
}
rownames(annual_regional_sum_matchmos_df) <- NULL
colnames(annual_regional_sum_matchmos_df) <- c("year", "region", "n", "mean", "sd", "lower", "upper")
annual_regional_sum_matchmos_df <- annual_regional_sum_matchmos_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_regional_sum_matchmos_df
# Emissions in Tg/y


# -------------------
# C40 Network Sum
# -------------------
annual_ch4_SRON$C40 <- apply(annual_ch4_SRON[c("city")], 1, function(x) is_C40(x))

annual_ch4_SRON_allyears <- annual_ch4_SRON %>%
                                filter(year != 2024) %>%
                                group_by(C40, region, city) %>%
                                filter(length(unique(year)) >= length(seq(2019, 2023))) %>%
                                ungroup()

annual_network_sum_df <- data.frame(NULL)
network_realizations <- data.frame(NULL)
for (C40 in unique(annual_ch4_SRON_allyears$C40)) {
  for (year in c(2019, 2020, 2021, 2022, 2023)) {
    bstrap <- c()
    ch4_year <- as.data.frame(annual_ch4_SRON_allyears[(annual_ch4_SRON_allyears$year == year) & (annual_ch4_SRON_allyears$C40 == C40), c(5:1004)])
    bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE)) # sum in each column bc already randomized
    network_realizations <- rbind(network_realizations,
                                  data.frame(year, C40, n = dim(ch4_year)[1], t(bstrap * tg_y_conversion)))
    annual_network_sum_df <- rbind(annual_network_sum_df,
                                   data.frame(year, C40, dim(ch4_year)[1],
                                              mean(bstrap, na.rm = TRUE), sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
  }
}
rownames(annual_network_sum_df) <- NULL
colnames(annual_network_sum_df) <- c("year", "C40", "n", "mean", "sd", "lower", "upper")
annual_network_sum_df <- annual_network_sum_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_network_sum_df


# -------------------
# Emissions of all cities that we have obs for over all years
# -------------------
annual_ch4_SRON_allyears_allcities <- annual_ch4_SRON %>%
                                select(-c(C40)) %>%
                                filter(year != 2024) %>%
                                group_by(region, city) %>% # same as above but without C40 differentiation
                                filter(length(unique(year)) >= length(seq(2019, 2023))) %>%
                                ungroup()

annual_sum_over_time_df <- data.frame(NULL) # all together
all_city_realizations <- data.frame(NULL)
for (year in c(2019, 2020, 2021, 2022, 2023)) {
  bstrap <- c()
  ch4_year <- as.data.frame(annual_ch4_SRON_allyears_allcities[(annual_ch4_SRON_allyears_allcities$year == year), c(5:1004)])
  bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE)) # sum in each column bc already randomized
  all_city_realizations <- rbind(all_city_realizations,
                                  data.frame(year, n = dim(ch4_year)[1], t(bstrap * tg_y_conversion)))
  annual_sum_over_time_df <- rbind(annual_sum_over_time_df,
                                  data.frame(year, dim(ch4_year)[1],
                                            mean(bstrap, na.rm = TRUE), sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
}

rownames(annual_sum_over_time_df) <- NULL
colnames(annual_sum_over_time_df) <- c("year", "n", "mean", "sd", "lower", "upper")
annual_sum_over_time_df <- annual_sum_over_time_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_sum_over_time_df
# Emissions in Tg/y


# -------------------
# C40 Network Sum MATCHMONTHS
# -------------------

annual_ch4_SRON_matchmos_allyears <- annual_ch4_SRON_matchmos %>%
                                filter(year != 2024) %>%
                                group_by(C40, region, city) %>%
                                filter(length(unique(year)) >= length(seq(2019, 2023))) %>%
                                ungroup()

annual_network_sum_matchmos_df <- data.frame(NULL)
for (C40 in unique(annual_ch4_SRON_matchmos_allyears$C40)) {
  for (year in c(2019, 2020, 2021, 2022, 2023)) {
    bstrap <- c()
    ch4_year <- as.data.frame(annual_ch4_SRON_matchmos_allyears[(annual_ch4_SRON_matchmos_allyears$year == year) & (annual_ch4_SRON_matchmos_allyears$C40 == C40), c(5:1004)])
    bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE)) # sum in each column bc already randomized

    annual_network_sum_matchmos_df <- rbind(annual_network_sum_matchmos_df,
                                   data.frame(year, C40, dim(ch4_year)[1],
                                              mean(bstrap, na.rm = TRUE), sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
  }
}
rownames(annual_network_sum_matchmos_df) <- NULL
colnames(annual_network_sum_matchmos_df) <- c("year", "C40", "n", "mean", "sd", "lower", "upper")
annual_network_sum_matchmos_df <- annual_network_sum_matchmos_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_network_sum_matchmos_df
# Emissions in Tg/y


# -------------------
# CH4 Emissions of All Cities by Year, full uncert
# -------------------
# not to be used with years together!
annual_total_sum_df <- data.frame(NULL)
for (year in c(2019, 2020, 2021, 2022, 2023)) {
  bstrap <- c()
  ch4_year <- annual_ch4_full[(annual_ch4_full$year == year), c(5:1004)]
  bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE))

  annual_total_sum_df <- rbind(annual_total_sum_df,
                          data.frame(year, dim(ch4_year)[1], mean(bstrap, na.rm = TRUE),
                                     sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
}
rownames(annual_total_sum_df) <- NULL
colnames(annual_total_sum_df) <- c("year", "n", "mean", "sd", "lower", "upper")
annual_total_sum_df <- annual_total_sum_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_total_sum_df
# Emissions in Tg/y


# -------------------
# CH4 Emissions of All Cities by Year, full uncert MATCHMONTHS
# -------------------
# not to be used with years together!
annual_total_sum_matchmos_df <- data.frame(NULL)
for (year in c(2019, 2020, 2021, 2022, 2023)) {
  bstrap <- c()
  ch4_year <- annual_ch4_full_matchmos[(annual_ch4_full_matchmos$year == year), c(5:1004)]
  bstrap <- c(bstrap, colSums(ch4_year, na.rm = TRUE))

  annual_total_sum_matchmos_df <- rbind(annual_total_sum_matchmos_df,
                          data.frame(year, dim(ch4_year)[1], mean(bstrap, na.rm = TRUE),
                                     sd(bstrap, na.rm = TRUE), quantile(bstrap, 0.025, na.rm = TRUE), quantile(bstrap, 0.975, na.rm = TRUE)))
}
rownames(annual_total_sum_matchmos_df) <- NULL
colnames(annual_total_sum_matchmos_df) <- c("year", "n", "mean", "sd", "lower", "upper")
annual_total_sum_matchmos_df <- annual_total_sum_matchmos_df %>%
                        mutate(mean = mean * tg_y_conversion, sd = sd * tg_y_conversion, lower = lower * tg_y_conversion, upper = upper * tg_y_conversion)
annual_total_sum_matchmos_df
# Emissions in Tg/y


# -------------------------------
# SAVE DATA
# -------------------------------
# city ch4co, ch4 emissions with change over time
paste0(path_to_data, "results/annual_ch4_emissions.csv")
write.table(city_annual_df, paste0(path_to_data, "results/city_ch4co_and_emissions.csv"), row.names = FALSE, sep = ",")

# regional sum emissions (change over time)
write.table(annual_regional_sum_df, paste0(path_to_data, "results/regional_ch4_emissions.csv"), row.names = FALSE, sep = ",")


# network sum emissions (change over time)
write.table(annual_network_sum_df, paste0(path_to_data, "results/network_ch4_emissions.csv"), row.names = FALSE, sep = ",")
write.table(annual_network_sum_matchmos_df, paste0(path_to_data, "results/network_ch4_emissions_matchmos.csv"), row.names = FALSE, sep = ",")

# all cities change over time (same set of cities each year)
write.table(annual_sum_over_time_df, paste0(path_to_data, "results/allcities_over_time_ch4_emissions.csv"), row.names = FALSE, sep = ",")



# all cities per year urban contribution (absolute uncertainty on EDGAR CO) (city count varies by year)
write.table(annual_total_sum_df, paste0(path_to_data, "results/annual_sum_ch4_emissions_fulluncert.csv"), row.names = FALSE, sep = ",")
write.table(annual_total_sum_matchmos_df, paste0(path_to_data, "results/annual_sum_ch4_emissions_matchmos_fulluncert.csv"), row.names = FALSE, sep = ",")


# save realizations for analysis:
# -------------------------------
# ch4co realizations
write.table(annual_ch4co_realizations, paste0(path_to_data, "results/realizations/reannual_ch4co_realizations.csv"), row.names = FALSE, sep = ",")


# ch4 emission realizations
write.table(annual_ch4_full, paste0(path_to_data, "results/realizations/annual_ch4_realizations_fulluncert.csv"), row.names = FALSE, sep = ",")
write.table(annual_ch4_SRON, paste0(path_to_data, "results/realizations/annual_ch4_realizations_SRONuncert.csv"), row.names = FALSE, sep = ",")

# ch4 emissions realizations grouping for difference distributions
write.table(network_realizations, paste0(path_to_data, "results/realizations/network_realizations.csv"), row.names = FALSE, sep = ",")
write.table(regional_realizations, paste0(path_to_data, "results/realizations/regional_realizations.csv"), row.names = FALSE, sep = ",")
write.table(all_city_realizations, paste0(path_to_data, "results/realizations/all_city_realizations.csv"), row.names = FALSE, sep = ",")