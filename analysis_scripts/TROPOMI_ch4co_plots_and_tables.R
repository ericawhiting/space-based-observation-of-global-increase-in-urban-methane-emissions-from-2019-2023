# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(kableExtra)
library(reshape2)
library(ggh4x)
library(tidyverse)

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))
source(paste0(path_to_utils, "plot_overpasses.R"))

path_to_data <- "../Data/" # change your path here

#path for output plots
path_to_plots <- "../plots/results/" # change your path here
# -------------------------------
# Plotting
# -------------------------------
regional_palette <- c("#9A9999", "#E69F00", "#56B4E9", "#009E73", "#d7cb24", "#0072B2", "#AA4499")
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

# -------------------------------
# SI TABLE
# -------------------------------
city_df <- read.csv(paste0(path_to_data, "results/city_ch4co_and_emissions.csv"))
city_df <- city_df %>% filter(city != "gyq") # not including this city in analysis

# load edgar v8 and edgar v2024 data
#####
edgar_v8_ch4_df <- data.frame()
edgar_v2024_ch4_df <- data.frame()
city_list <- unique(city_df$city)
for (city in city_list) {
    # edgar v8 CH4
    ch4_v8_file_name_format <- list.files(paste0(path_to_data, "EDGAR_ch4_emissions/annual"),
                                           pattern = paste0(city, "_annual_ch4_emissions_from_monthly_edgarv8.0.csv"), full.names = TRUE)
    if (length(ch4_v8_file_name_format) == 1) { # essentially if file exists
        edgar_v8_ch4_df_temp <- read.csv(ch4_v8_file_name_format)
        edgar_v8_ch4_df_temp$city <- city
        edgar_v8_ch4_df_temp <- edgar_v8_ch4_df_temp %>% rename(EDGARv8_CH4 = domain_mean)
        edgar_v8_ch4_df <- rbind(edgar_v8_ch4_df, edgar_v8_ch4_df_temp)
    }
    # edgar 2024 ch4
    ch4_v2024_file_name_format <- list.files(paste0(path_to_data, "EDGAR_ch4_emissions/annual"),
                                           pattern = paste0(city, "_annual_ch4_emissions_from_monthly_edgarv2024.csv"), full.names = TRUE)
    if (length(ch4_v2024_file_name_format) == 1) { # essentially if file exists
        edgar_v2024_ch4_df_temp <- read.csv(ch4_v2024_file_name_format)
        edgar_v2024_ch4_df_temp$city <- city
        edgar_v2024_ch4_df_temp <- edgar_v2024_ch4_df_temp %>% rename(EDGARv2024_CH4 = domain_mean)
        edgar_v2024_ch4_df <- rbind(edgar_v2024_ch4_df, edgar_v2024_ch4_df_temp)
    }
}

tg_y_conversion <- 3600 * 24 * 365.25 * 1e-9 # from kg/s
edgar_ch4_df <- merge(edgar_v8_ch4_df, edgar_v2024_ch4_df, by = c("city", "year"), all = TRUE) %>%
                    mutate(EDGARv8_CH4 = EDGARv8_CH4 * tg_y_conversion,
                           EDGARv2024_CH4 = EDGARv2024_CH4 * tg_y_conversion) %>%
                    filter(year %in% c(2019, 2020, 2021, 2022, 2023)) %>%
                    filter(city %in% unique(city_df$city))
#####
# Load City Information
city_df <- merge(city_df, edgar_ch4_df, by = c("year", "city"), all.x = TRUE, all.y = TRUE)
city_df$region <- apply(city_df[c("city")], 1, function(x)  find_city_region(x))
city_df$C40 <- apply(city_df[c("city")], 1, function(x) is_C40(x))
city_df$name <- apply(city_df[c("city")], 1, function(x) find_city_name(x))
city_df[c("lat", "lon")] <- t(sapply(city_df$city, function(city) unlist(find_city_center_from_pop_dens(city))))
city_df$delta_lat <- apply(city_df[c("city")], 1, function(x) read_ideal_boxsize_lat(x))
city_df$coastal_flag <- apply(city_df[c("city")], 1, function(x) find_city_coastal_flag(x))
city_df$delta_lat <- city_df$delta_lat + 0.25 * city_df$delta_lat * city_df$coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land

# ch4 edgar totals
# this does not match initial paper submission due to edgar updates
city_df %>% filter(num_obs > 0) %>% group_by(year) %>% summarize(n = n(), total_v8 = sum(EDGARv8_CH4), total_v2024 = sum(EDGARv2024_CH4))


# Helper functions
construct_uncertainty <- function(value, lower, upper) {
  if (any(is.na(c(value, lower, upper)))) {
    return(NA_character_)
  }
  return(paste0(value, " (", lower, ", ", upper, ")"))
}

add_degree <- function(x) {
  if_else(is.na(x), NA_character_, paste0(x, "°"))
}

si_table_data <- city_df  %>% select(c(name, region, C40, lat, lon, delta_lat, year, num_obs, ch4co, ch4co_lower, ch4co_upper, ch4_sron, ch4_sron_lower, ch4_sron_upper, EDGARv8_CH4, EDGARv2024_CH4))
si_table_data[si_table_data$region == "South Africa", "region"] <- "Southern Africa"

si_table_data <- si_table_data %>% mutate(across(where(is.numeric), ~ round(.x, 2)))
si_table_data$C40 <- replace(si_table_data$C40, si_table_data$C40 == 0, "No")
si_table_data$C40 <- replace(si_table_data$C40, si_table_data$C40 == 1, "Yes")
si_table_data <- si_table_data %>%
  mutate(delta_lat = add_degree(delta_lat))
edgar_v8 <- si_table_data %>%
             select(name, region, lat, lon, delta_lat, C40, year, EDGARv8_CH4) %>%
             mutate(EDGARv8_CH4 = if_else(year >  2022, NA, EDGARv8_CH4)) %>%
             pivot_wider(names_from = year, values_from = EDGARv8_CH4) %>%
             mutate(Type = "EDGARv8 CH4 [Tg]") %>%
             mutate(across(`2019`:`2023`, as.character))
edgar_v2024 <- si_table_data %>%
             select(name, region, lat, lon, delta_lat, C40, year, EDGARv2024_CH4) %>%
             pivot_wider(names_from = year, values_from = EDGARv2024_CH4) %>%
             mutate(Type = "EDGARv2024 CH4 [Tg]") %>%
             mutate(across(`2019`:`2023`, as.character))

n_data <- si_table_data %>%
            select(name, region, lat, lon, delta_lat, C40, year, num_obs) %>%
            pivot_wider(names_from = year, values_from = num_obs) %>%
            mutate(Type = "# Obs")
n_data[is.na(n_data)] <- 0
n_data <- n_data %>%
  mutate(across(c("2019", "2020", "2021", "2022", "2023"), as.character))

ch4co_data <- si_table_data %>%
  select(name, region, lat, lon, delta_lat, C40, year, ch4co, ch4co_lower, ch4co_upper) %>%
  pivot_wider(names_from = year, values_from = c(ch4co, ch4co_lower, ch4co_upper)) %>%
  rowwise() %>%
  mutate(ch4co_2019 = construct_uncertainty(ch4co_2019, ch4co_lower_2019, ch4co_upper_2019),
         ch4co_2020 = construct_uncertainty(ch4co_2020, ch4co_lower_2020, ch4co_upper_2020),
         ch4co_2021 = construct_uncertainty(ch4co_2021, ch4co_lower_2021, ch4co_upper_2021),
         ch4co_2022 = construct_uncertainty(ch4co_2022, ch4co_lower_2022, ch4co_upper_2022),
         ch4co_2023 = construct_uncertainty(ch4co_2023, ch4co_lower_2023, ch4co_upper_2023),
         data_all_years = !is.na(ch4co_2019) & !is.na(ch4co_2020) & !is.na(ch4co_2021) & !is.na(ch4co_2022) & !is.na(ch4co_2023)) %>%
  ungroup() %>%
  select(-c(data_all_years, ch4co_lower_2019, ch4co_upper_2019, ch4co_lower_2020, ch4co_upper_2020,
            ch4co_lower_2021, ch4co_upper_2021, ch4co_lower_2022, ch4co_upper_2022, ch4co_lower_2023, ch4co_upper_2023)) %>%
  rename("2019" = ch4co_2019, "2020" = ch4co_2020, "2021" = ch4co_2021, "2022" = ch4co_2022, "2023" = ch4co_2023) %>%
  mutate(Type = "CH4:CO")

ch4_emissions_data <- si_table_data %>%
  rename(ch4 = ch4_sron, ch4_lower = ch4_sron_lower, ch4_upper = ch4_sron_upper) %>%
  select(name, region, lat, lon, delta_lat, C40, year, ch4, ch4_lower, ch4_upper) %>%
  pivot_wider(names_from = year, values_from = c(ch4, ch4_lower, ch4_upper)) %>%
  rowwise() %>%
  mutate(ch4_2019 = construct_uncertainty(ch4_2019, ch4_lower_2019, ch4_upper_2019),
         ch4_2020 = construct_uncertainty(ch4_2020, ch4_lower_2020, ch4_upper_2020),
         ch4_2021 = construct_uncertainty(ch4_2021, ch4_lower_2021, ch4_upper_2021),
         ch4_2022 = construct_uncertainty(ch4_2022, ch4_lower_2022, ch4_upper_2022),
         ch4_2023 = construct_uncertainty(ch4_2023, ch4_lower_2023, ch4_upper_2023),
         data_all_years = !is.na(ch4_2019) & !is.na(ch4_2020) & !is.na(ch4_2021) & !is.na(ch4_2022) & !is.na(ch4_2023)) %>%
  ungroup() %>%

  select(-c(data_all_years, ch4_lower_2019, ch4_upper_2019, ch4_lower_2020, ch4_upper_2020, ch4_lower_2021, ch4_upper_2021, ch4_lower_2022, ch4_upper_2022, ch4_lower_2023, ch4_upper_2023)) %>%
  rename("2019" = ch4_2019, "2020" = ch4_2020, "2021" = ch4_2021, "2022" = ch4_2022, "2023" = ch4_2023) %>%
  mutate(Type = "CH4 [Tg]")

si_table <- bind_rows(n_data, ch4co_data, ch4_emissions_data, edgar_v8, edgar_v2024) %>% arrange(name)

df_of_annual_inv_data <- si_table

si_table <- si_table %>%
  select(name, region, C40, lat, lon, delta_lat, Type, "2019", "2020", "2021", "2022", "2023")

# move region and lat lon into row 2 and 3
si_table <- si_table %>%
  group_by(name) %>%
  mutate(name = case_when(row_number() == 1 ~ as.character(name),
                          row_number() == 2 ~ paste0(as.character(region)),
                          row_number() == 3 ~ paste0("(", as.character(lat), "°, ", as.character(lon), "°)"),
                          row_number() == 4 ~ "",
                          row_number() == 5 ~ ""),
         delta_lat = case_when(row_number() == 1 ~ as.character(delta_lat),
                               row_number() != 1 ~ ""),
         C40 = case_when(row_number() == 1 ~ as.character(C40),
                         row_number() != 1 ~ "")
                            ) %>%
  ungroup() %>%
  select(-c("region", "lat", "lon"))

# Alternating color logic
si_table <- si_table %>%
  mutate(row_pair = rep(rep(1:2, each = 5), length.out = nrow(si_table)),
         background = ifelse(row_pair %% 2 == 0, "gray", "white"))

si_table <- si_table %>%
   rename("Name, Region, Location" = name, "Box" = delta_lat)
# Generate the LaTeX table without the auxiliary columns
latex_table <- kable(si_table %>% select(-row_pair, -background),
                     format = "latex", booktabs = TRUE, longtable = TRUE,
                     caption = "CH4CO and CH4 Emissions Data", linesep = "") %>%
  kable_styling(latex_options = c("hold_position", "repeat_header"), full_width = FALSE, font_size = 7) %>%
  row_spec(which(si_table$background == "gray"), background = "#dad6d6") %>%
  row_spec(which(si_table$background == "white"), background = "white") %>%
  row_spec(0, bold = TRUE)

# Save the table
save_kable(latex_table, file = paste0(path_to_data, "results/SI_table.tex"), self_contained = TRUE)

# -------------------------------
# FIGURE 1: Map of Cities
# -------------------------------
city_df <- read.csv(paste0(path_to_data, "results/city_ch4co_and_emissions.csv")) %>% filter(city != "gyq") # not including this city in analysis
city_df$C40 <- apply(city_df[c("city")], 1, function(x) is_C40(x))
city_df$name <- apply(city_df[c("city")], 1, function(x) find_city_name(x))
city_df[c("lat", "lon")] <- t(sapply(city_df$city, function(city) unlist(find_city_center_from_pop_dens(city))))
city_df$delta_lat <- apply(city_df[c("city")], 1, function(x) read_ideal_boxsize_lat(x))
city_df$coastal_flag <- apply(city_df[c("city")], 1, function(x) find_city_coastal_flag(x))
city_df$delta_lat <- city_df$delta_lat + 0.25 * city_df$delta_lat * city_df$coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land

ggplot(data = city_df %>%
group_by(city) %>%
  filter(length(unique(year)) == length(unique(city_df$year))) %>%
ungroup() %>%
group_by(region, year) %>%
summarize(mean_ch4co = mean(ch4co, na.rm = TRUE))) +
geom_line(aes(x = year, y = mean_ch4co, color = region), linewidth = 1.5) +
geom_point(aes(x = year, y = mean_ch4co, fill = region, color = region), size = 3) +
scale_color_manual(values = regional_palette) +
scale_fill_manual(values = regional_palette) + theme_linedraw() + plot_theme_text


# sort into c40 all years, c40 not all years, nonc40 all years
ch4_completedata <- city_df %>%
  group_by(city) %>%
  filter(length(unique(year)) == length(unique(city_df$year))) %>%
  filter(C40 == 1) %>%
  ungroup()
ch4_completedata$complete <- 1
ch4_incompletedata <- city_df %>%
  group_by(city) %>%
  filter(length(unique(year)) != length(unique(city_df$year))) %>%
  filter(C40 == 1) %>%
  ungroup()
ch4_incompletedata$complete <- 0
non_c40_cities <- city_df %>%
  group_by(city) %>%
  filter(length(unique(year)) == length(unique(city_df$year))) %>%
  filter(C40 == 0) %>%
  ungroup()

# plot world map of c40 cities, note if cities have data across all years (>= 5 data points per year)
ch4_data_combined <- bind_rows(
  dplyr::select(ch4_completedata, lon, lat, region, complete),
  dplyr::select(ch4_incompletedata, lon, lat, region, complete)
)
merged_c40 <- merge(ch4_completedata %>% group_by(region) %>% summarize(n = length(unique(city))),
                    ch4_incompletedata %>% group_by(region) %>% summarize(n = length(unique(city))),
                  by = "region", all = TRUE)

merged_c40[is.na(merged_c40$n.y), "n.y"] <- 0
merged_c40$n <-  merged_c40$n.x + merged_c40$n.y
merged_c40

ch4_data_combined$region <- factor(ch4_data_combined$region,
                                      levels = c("Central East Asia", "Africa", "East, Southeast Asia & Oceania", "North America", "Latin America", "South & West Asia",  "Europe"),
                                      labels = c("Central East Asia (9, 7)", "Africa (6, 2)", "East, Southeast Asia & Oceania (6, 6)", "North America (16, 11)", "Latin America (8, 5)",
                                                  "South & West Asia (10, 8)",  "Europe (16, 12)"))

world_map <- ggplot() + borders("world", color = "black", fill = "white")
city_map <- world_map +
  theme_linedraw() +
  geom_point(data = non_c40_cities, aes(x = lon, y = lat),
             fill = "black", color = "black", shape = 21, size = 3) +
  geom_point(data = ch4_data_combined, aes(x = lon, y = lat, fill = region, shape = as.factor(complete)),
             color = "black", stroke = 1, size = 4) +

  scale_fill_manual(values = c("#ff9924", "#9A9999", "#56B4E9", "#0072B2",  "#d7cb24", "#AA4499", "#009E73"),
                    guide = guide_legend(override.aes = list(shape = 22, size = 8))) +  # Force the legend to display as colors
  scale_shape_manual(values = c(24, 21), labels = c("One or More Years of Data", "Data for All Years (2019-2023)"),
                     name = "Data Availability") +
  scale_x_continuous(breaks = c(-180, -90, 0, 90, 180), labels = function(l) sprintf("%1.f°", l)) +
  scale_y_continuous(breaks = c(90, -60, -30, 0, 30, 60, 90), labels = function(l) sprintf("%1.f°", l)) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 70), expand = FALSE) +
  plot_theme_text +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  labs(x = "Longitude", y = "Latitude", fill = "C40 Region", shape = "Years of Data", title = "(a)") +
  theme(plot.title.position = "plot")
city_map
ggsave("worldmap_cities.png", plot = city_map, device = "png", path = path_to_plots, width = 19, height = 9, dpi = 320, bg = "white")

plot_satellite_overpass("nyc", "21073", "b", "c")
plot_satellite_overpass("del", "21720", "d", "e")

# Read in CH4 CO Data for C40 cities
#ch4co_data_path <- paste0(path_prefix, "Data/C40_ch4co_fire/C40_ch4co_fire/")
ch4co_data_path <- paste0(path_to_data, "ch4co/") # changed bc fire filtering is now within analysis rather than separate in initial analysis period

no_pass_list <- c() # for cities that had 0 overpasses :(
ch4co_df <- data.frame() # initialize df for combining all C40 city information
# iterate through all csv data from city analysis, append to df, add city code, name, and region
for (file in list.files(ch4co_data_path, pattern = "*.csv")) {
  city <- str_split(file, "_")[[1]][1]
  region <- find_city_region(city)
  name <- find_city_name(city)
  latlon <- find_city_center_from_pop_dens(city)
  temp_ch4co_df <- read.csv(paste0(ch4co_data_path, file), header = TRUE)
  temp_ch4co_df  <- temp_ch4co_df %>% drop_na() # #, -slope_results_2)))
  temp_ch4co_df$datetimes <- as.Date(temp_ch4co_df$datetimes)
  if (dim(temp_ch4co_df)[1] > 0) {
    temp_ch4co_df$city <- city
    temp_ch4co_df$region <- region
    temp_ch4co_df$name <- name
    temp_ch4co_df$lat <- latlon[1]
    temp_ch4co_df$lon <- latlon[2]

  } else {
    no_pass_list <- append(no_pass_list, city)
  }
  ch4co_df <- rbind(ch4co_df, temp_ch4co_df)
}

# Metropolitan Areas / dealing with double counting
ch4co_df <- ch4co_df[(ch4co_df$city != "tok") & (ch4co_df$city != "eku") & (ch4co_df$city != "tsh") & (ch4co_df$city != "hk") & (ch4co_df$city != "rot") & (ch4co_df$city != "zhe"), ]

# Add Datetime information: year, month, season
ch4co_df$year <- as.numeric(format(as.Date(ch4co_df$datetimes), "%Y"))
ch4co_df$month <- as.numeric(format(as.Date(ch4co_df$datetimes), "%m"))
ch4co_df$yrmo <- format(as.Date(paste(ch4co_df$year, ch4co_df$month, "01", sep = "-")))

# remove cities with < 5 overpasses per year
removed_cities_5 <- ch4co_df %>% group_by(region, city, name, year) %>% summarise(sum_overpasses = sum(num_overpasses)) %>% filter(sum_overpasses < 5)
ch4co_df <- anti_join(ch4co_df, removed_cities_5, by = c("city", "year"))

#NYC
nyc_ch4co_df <- ch4co_df[ch4co_df$city == "nyc", ]
start_date <- as.Date("2019-01-01")
end_date <- as.Date("2024-01-01")
ch4co_nyc <- ggplot(data = nyc_ch4co_df) +
  geom_point(aes(x = datetimes, y = integ_results_15), color = "#0072B2", fill = "#0072B2", size = 3) +
  scale_x_date(breaks = seq(from = start_date, to = end_date, by = "1 year"),
               limits = c(start_date, end_date),
               date_labels = "%Y-%m") +
  theme_linedraw() +
  plot_theme_text +
  labs(title = "(f) Enhancement Ratio", x = "Date", y = bquote(bold(frac(Sigma ~ (CH[4])[enh], Sigma ~ (CO)[enh])))) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16),  axis.title.y = element_text(size = 18),
        plot.title.position = "plot", plot.margin = unit(c(0.05, 0.21, 0.05, 0.05), "inches"))
ch4co_nyc
ggsave("ch4co_nyc.png", plot = ch4co_nyc, device = "png", path = path_to_plots, width = 12, height = 3, dpi = 320, bg = "white")

# DEL

del_ch4co_df <- ch4co_df[ch4co_df$city == "del", ]
start_date <- as.Date("2019-01-01")   # Start on March 1, 2023
end_date <- as.Date("2024-01-01")
ch4co_del <- ggplot(data = del_ch4co_df) +
  geom_point(aes(x = datetimes, y = integ_results_15), color = "#AA4499", fill = "#AA4499", size = 3) +
  scale_x_date(breaks = seq(from = start_date, to = end_date, by = "1 year"),
               limits = c(start_date, end_date),
               date_labels = "%Y-%m") +
  theme_linedraw() +
  plot_theme_text +
  labs(title = "(g) Enhancement Ratio", x = "Date", y = bquote(bold(frac(Sigma ~ (CH[4])[enh], Sigma ~ (CO)[enh])))) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 18),
        plot.title.position = "plot",
        plot.margin = unit(c(0.05, 0.21, 0.05, 0.05), "inches"))
ch4co_del
ggsave("ch4co_del.png", plot = ch4co_del, device = "png", path = path_to_plots, width = 12, height = 3, dpi = 320, bg = "white")

# -------------------------------
# FIGURE 2: City CH4CO and CH4 Emissions
# -------------------------------
city_df <- read.csv(paste0(path_to_data, "results/city_ch4co_and_emissions.csv"))
city_df$C40 <- apply(city_df[c("city")], 1, function(x) is_C40(x))
city_df$name <- apply(city_df[c("city")], 1, function(x) find_city_name(x))
city_df[c("lat", "lon")] <- t(sapply(city_df$city, function(city) unlist(find_city_center_from_pop_dens(city))))
city_df$delta_lat <- apply(city_df[c("city")], 1, function(x) read_ideal_boxsize_lat(x))
city_df$coastal_flag <- apply(city_df[c("city")], 1, function(x) find_city_coastal_flag(x))
city_df$delta_lat <- city_df$delta_lat + 0.25 * city_df$delta_lat * city_df$coastal_flag # make box slightly larger for coastal cities so less sensitive to how many pixels are lost to land

test_cases <- city_df[city_df$city %in% c("dur", "joh", "bei", "shg", "seo", "yok", "ams", "par", "rio", "sao", "chi", "nyc", "dha", "del"), ] %>%
              select(city, year, region, ch4co, ch4co_lower, ch4co_upper, ch4_sron, ch4_sron_lower, ch4_sron_upper)
test_cases_relative <- test_cases %>%
                       group_by(city) %>%
                       mutate(ch4_relative = ch4_sron / first(ch4_sron[year == 2019]),
                              lower_lim_rel = ch4_sron_lower / ch4_sron * ch4_relative,
                              upper_lim_rel = ch4_sron_upper / ch4_sron * ch4_relative,
                              ch4_relative = (ch4_relative - 1) * 100,
                              lower_lim_rel = (lower_lim_rel - 1) * 100,
                              upper_lim_rel = (upper_lim_rel - 1) * 100)

test_cases_ch4co <- melt(test_cases_relative %>%
                        select(city, year, region, ch4co, ch4_relative) %>%
                        rename(CH4_from_EDGARv8.1_CO_relative =  ch4_relative), id = c("city", "year", "region")) %>%
                    rename(ch4co_variable = variable, ch4co_value = value)
test_cases_ch4co$city <- factor(test_cases_ch4co$city, levels = c("bei", "shg", "dur", "joh", "seo", "yok", "chi", "nyc",  "rio", "sao",  "dha", "del", "ams", "par"))
test_cases_lower <- melt(test_cases_relative %>%
                          select(city, year, region, ch4co_lower, lower_lim_rel) %>%
                          rename(ch4co = ch4co_lower, CH4_from_EDGARv8.1_CO_relative = lower_lim_rel), id = c("city", "year", "region")) %>%
                    rename(ch4co_variable = variable, lower_value = value)
test_cases_upper <- melt(test_cases_relative %>%
                          select(city, year, region, ch4co_upper, upper_lim_rel) %>%
                          rename(ch4co = ch4co_upper, CH4_from_EDGARv8.1_CO_relative = upper_lim_rel), id = c("city", "year", "region")) %>%
                    rename(ch4co_variable = variable, upper_value = value)

test_cases_reformat_1 <- merge(test_cases_ch4co, test_cases_lower, on =  c("city", "year", "region", "ch4co_variable"))
test_cases_reformat <- merge(test_cases_reformat_1, test_cases_upper, on =  c("city", "year", "region", "ch4co_variable"))

levels(test_cases_reformat$ch4co_variable)[levels(test_cases_reformat$ch4co_variable) == "ch4co"] <- "CH4:CO"
levels(test_cases_reformat$ch4co_variable)[levels(test_cases_reformat$ch4co_variable) == "CH4_from_EDGARv8.1_CO_relative"] <- "CH4 Emissions\n(relative to 2019 in %)"
# Relevel the factor
test_cases_reformat$city <- factor(test_cases_reformat$city,
                              levels = c("bei", "shg", "dur", "joh", "seo", "yok", "chi", "nyc",  "rio", "sao",  "dha", "del", "ams", "par"),
                              labels = c("Beijing          ",
                                        "Shanghai          ",
                                        "Durban         ",
                                        "Johannesburg   ",
                                        "Seoul            ",
                                        "Yokohama         ",
                                        "Chicago         ",
                                        "New York City   ",
                                        "Rio de Janeiro   ",
                                        "São Paulo        ",
                                        "Dhaka               ",
                                        "Delhi               ",
                                        "Amsterdam      ",
                                        "Paris          "))
test_cases_reformat$region <- factor(test_cases_reformat$region,
                                      levels = c("Central East Asia", "Africa", "East, Southeast Asia & Oceania", "North America", "Latin America", "South & West Asia",  "Europe"),
                                      labels = c("Central East Asia", "Southern Africa", "East, Southeast Asia & Oceania", "North America", "Latin America", "South & West Asia",  "Europe"))


test_case_plot <- ggplot(data =  test_cases_reformat) +
  facet_grid(vars(ch4co_variable), vars(region), scales = "free_y", switch = "y", labeller = labeller(region = label_wrap_gen(width = 15.1), .rows = NULL)) +
  geom_ribbon(aes(x = year, y = ch4co_value, ymin = (lower_value), ymax = (upper_value), group = city, color = city, fill = city), alpha = 0.65) +
  geom_line(aes(x = year, y = ch4co_value, group = city, color = city), linewidth = 1.0, alpha = 1) +
  geom_point(aes(x = year, y = ch4co_value, group = city, color = city, fill = city, shape = city), size = 3, alpha = 1) +
  scale_color_manual(values = c("#6e3302", "#ffcd94", "#1311118c", "#fffdfd",
                                "#000000", "#a0dcff", "#1b264d",   "#3596cf", "#28251c", "#fff784",
                                "#fcacef", "#461f40", "#01110d", "#5fe8c4")) +
  scale_fill_manual(values = c("#ab4d00", "#ff9924", "#2d2d2d8c", "#C6C5C5",
                                "#08283b", "#6fc3f3", "#0d1841",  "#035c8f", "#756426", "#f9eb27",
                                "#f45cdb", "#AA4499", "#023d2d", "#009E73")) +
  scale_shape_manual(values = c(21, 24, 21, 24, 21, 24, 21, 24, 21, 24, 21, 24, 21, 24)) +
  theme_linedraw() + plot_theme_text +
  facetted_pos_scales(
    y = list(scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3)),
            scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75), limits = c(-80, 80), labels = function(x) paste0(x, "%")))) +
  guides(y = NULL, color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(x = "Year", y = bquote(bold(atop("CH"[4] ~ " Emissions                    ", "(relative to 2019 in %)"~~~~~~~~"CH"[4]:"CO"))),
      color = "", fill = "", shape = "") +
  theme(axis.text.x = element_text(angle = 90, size = 14),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 16, face = "bold", margin = margin(5, 5, 50, 5)),
        text = element_text(family = "Times"),
        legend.direction = "horizontal",
        legend.position = c(0.5, 1.06), legend.background = element_blank(),
        axis.title.y = element_text(vjust = 0.5, hjust = -2.3, margin = margin(t = 30), size = 16),
            legend.text = element_text(color = "white", size = 14), strip.text.y = element_blank())#, #margin = margin(r = 16, unit = "pt")))
test_case_plot
ggsave("figure_2.png", plot = test_case_plot, device = "png", path = path_to_plots, width = 16, height = 7, dpi = 320, bg = "white")

# print ch4 emissions to add onto figure in powerpoint
test_cases %>% filter(year == 2019) %>% select(city, region, ch4_sron) %>% arrange(region)

# -------------------------------
# FIGURE 3: Regional CH4 Emissions
# -------------------------------
regional_sum_df <- read.csv(paste0(path_to_data, "results/regional_ch4_emissions.csv"))
# cities with 5 years of data, c40 only

regional_sum_df_relative <- regional_sum_df %>%
                              group_by(region) %>%
                              mutate(ch4_relative = mean / first(mean[year == 2019]),
                                     lower_lim_rel = lower / mean * ch4_relative,
                                     upper_lim_rel = upper / mean * ch4_relative,
                                     ch4_relative = (ch4_relative - 1) * 100,
                                     lower_lim_rel = (lower_lim_rel - 1) * 100,
                                     upper_lim_rel = (upper_lim_rel - 1) * 100)

# make Africa label --> South Africa
regional_sum_df_relative$region <- factor(regional_sum_df_relative$region,
                                      levels = c("Central East Asia", "Africa", "East, Southeast Asia & Oceania", "North America", "Latin America", "South & West Asia",  "Europe"),
                                      labels = c("Central East Asia (7)", "Southern Africa (2)", "East, Southeast Asia & Oceania (6)", "North America (11)", "Latin America (5)", "South & West Asia (8)",  "Europe (12)"))
regional_sum_plot <- ggplot(data = regional_sum_df_relative, aes(x = year, y = ch4_relative, color = region, fill = region, group = region)) +
  geom_ribbon(aes(ymin = (lower_lim_rel), ymax = (upper_lim_rel)), alpha = 0.5) +
  scale_color_manual(values = c("#ff9924", "#9A9999", "#56B4E9", "#0072B2",  "#d7cb24", "#AA4499", "#009E73")) +
  scale_fill_manual(values = c("#ff9924", "#9A9999", "#56B4E9", "#0072B2",  "#d7cb24", "#AA4499", "#009E73")) +
  facet_wrap(~region, nrow = 1,  labeller = labeller(region = label_wrap_gen(width = 15))) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(-40, 44), breaks = c(-30, -15, 0, 15, 30)) +
  geom_line(linewidth = 2.5) +
  geom_point(size = 3.5) +
  theme_linedraw() +
  plot_theme_text + theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
  axis.title.x = element_text(vjust = 0.15)) +
  guides(fill = "none", color = "none") +
  labs(x = "Year", y = bquote(bold(atop("CH"[4] ~ " Emissions", "(relative to 2019 in %)"))), fill = "C40 Region", color = "C40 Region") +
  theme(strip.text = element_text(face = "bold", size = 16),
        text = element_text(family = "Times")) +
  theme(panel.spacing = unit(1, "lines")) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  theme(plot.margin = unit(c(0.05, 0.21, 0.05, 0.05),
                                "inches")) +
  theme(plot.title = element_text(hjust = 0, vjust = -4), plot.title.position = "plot")

regional_sum_plot

ggsave("figure_3.png", plot = regional_sum_plot, device = "png", path = path_to_plots, width = 20, height = 4.7, dpi = 320, bg = "white")

# print ch4 emissions to add onto figure in powerpoint
regional_sum_df %>% filter(year == 2019) %>% select(region, mean) %>% arrange(region)

# -------------------------------
# FIGURE 4:  Network CH4 Emissions
# -------------------------------
# network sum emissions (change over time)
network_sum_df <- read.csv(paste0(path_to_data, "results/network_ch4_emissions.csv"))

network_sum_df
change_df <- network_sum_df %>%
  group_by(C40) %>%
    filter(year != 2019) %>%
    mutate(change = 1 + ((mean - lag(mean)) / lag(mean))) %>%
  ungroup() %>%
  filter(!is.na(change))
geometric_mean(as.vector(change_df[change_df$C40 == 1, "change"]$change))

network_sum_df_relative <- network_sum_df %>%
                              group_by(C40) %>%
                              mutate(ch4_relative = mean / first(mean[year == 2019]),
                                     lower_lim_rel = lower / mean * ch4_relative,
                                     upper_lim_rel = upper / mean * ch4_relative,
                                     ch4_relative = (ch4_relative - 1) * 100,
                                     lower_lim_rel = (lower_lim_rel - 1) * 100,
                                     upper_lim_rel = (upper_lim_rel - 1) * 100)

network_sum_df_relative$C40 <- factor(network_sum_df_relative$C40, levels = c(0, 1), labels = c("Non-C40 Cities (21)", "C40 Cities (51)"))

network_sum_plot <- ggplot(data = network_sum_df_relative, aes(x = year, y = ch4_relative, color = as.factor(C40), shape = as.factor(C40), fill = as.factor(C40), group = as.factor(C40))) +
  geom_ribbon(aes(ymin = (lower_lim_rel), ymax = (upper_lim_rel)), alpha = 0.4) +
  scale_color_manual(values = c("red", "blue")) +
  scale_fill_manual(values = c("red", "blue")) +
  scale_shape_manual(values = c(21, 24)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(-16, 16), breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_line(linewidth = 2) +
  geom_point(size = 4) +
  theme_linedraw() +
  plot_theme_text + theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
  axis.title.x = element_text(vjust = 0.15)) +
  labs(x = "Year", y = bquote(bold(atop("CH"[4] ~ " Emissions", "(relative to 2019 in %)"))), fill = "C40", color = "C40", shape = "C40") +
  theme(strip.text = element_text(face = "bold", size = 18),
        text = element_text(family = "Times"),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.margin = unit(c(0.05, 0.21, 0.05, 0.05), "inches"),
        plot.title = element_text(hjust = 0, vjust = -4), plot.title.position = "plot",
        legend.position = c(0.76, 0.15), legend.background = element_rect(fill = "white", linewidth = 0.5, linetype = "solid", color = "black"))

network_sum_plot
ggsave("figure_4.png", plot = network_sum_plot, device = "png", path = path_to_plots, width = 9, height = 7, dpi = 320, bg = "white")

network_sum_df %>% filter(year == 2019) %>% select(year, C40, mean)
network_sum_df %>% filter(year == 2023) %>% select(year, C40, mean)

# 2020-2023
# nonc40: 11.99% --> 12%
# c40: 0.095% --> 10%

# -------------------------------
# urban contribution (total error)
# -------------------------------
urban_contribution <- read.csv(paste0(path_to_data, "results/annual_sum_ch4_emissions_fulluncert.csv"))
urban_contribution

# -------------------------------
# allcities sum df:
# -------------------------------
# allcities sum emissions (change over time across all cities there is consistent data for)
allcities_sum_total_df <- read.csv(paste0(path_to_data, "results/allcities_over_time_ch4_emissions.csv"))
allcities_sum_total_df


geometric.mean <- function(x, na.rm = TRUE) {
  exp(mean(log(x), na.rm = na.rm))
}

mean(c(1 + (23.5-23.2)/23.2, 1 + (23.8 - 23.5)/ 23.5, 1 + (23.9 - 23.8)/ 23.8))
#<1%
mean(c(1 + (23.5-24.3)/ 24.3, 1 + (24.4-23.5)/23.5, 1 + (25.2 - 24.4)/24.4, 1 + (25.9-25.2)/25.2))
#1.6%, though changing by often 3% each year
# EDGAR
geometric.mean(c(1 + (23.5-23.2)/23.2, 1 + (23.8 - 23.5)/ 23.5, 1 + (23.9 - 23.8)/ 23.8)) -1
# CH4 Obs
geometric.mean(c(1 + (23.5-24.3)/ 24.3, 1 + (24.4-23.5)/23.5, 1 + (25.2 - 24.4)/24.4, 1 + (25.9-25.2)/25.2)) -1
# 3.3% if absolute value of change

# change
change_df <- si_table_data %>%
    filter(C40 == "No")%>%
    filter(!is.na(num_obs)) %>%
    group_by(name) %>%
    filter(length(unique(year)) == length(unique(si_table_data$year))) %>%
    ungroup() %>%
   # filter(year != 2023) %>%
    group_by(name) %>%
    arrange(year, .by_group = TRUE) %>%
    mutate(pct_change_8 = (EDGARv8_CH4 - lag(EDGARv8_CH4)) / lag(EDGARv8_CH4) * 100,
           pct_change = (ch4 - lag(ch4)) / lag(ch4) * 100)

total_pct_change_df <- change_df %>%
  group_by() %>%
  summarize(avg_pct_change_inv = mean(pct_change_8, na.rm = TRUE), sd_pct_change_inv = sd(pct_change_8, na.rm = TRUE),
            avg_pct_change_obs = mean(pct_change, na.rm = TRUE), sd_pct_change_obs = sd(pct_change, na.rm = TRUE))
total_pct_change_df

# 2.35% (19.0%) for non C40 cities
# 2.45% (16.2%) for C40 cities

# 1.35% (7.12%) for inventory of c40 cities

# ------------
# DIFFERENCE DISTRIBUTIONS
# ch4 emissions realizations grouping for difference distributions
network_realizations <- read.csv(paste0(path_to_data, "results/realizations/network_realizations.csv"))
regional_realizations <- read.csv(paste0(path_to_data, "results/realizations/regional_realizations.csv"))
annual_ch4_SRON <- read.csv(paste0(path_to_data, "results/realizations/annual_ch4_realizations_SRONuncert.csv"))
all_city_realizations <- read.csv(paste0(path_to_data, "results/realizations/all_city_realizations.csv"))
# histograms for each year
df_long <- network_realizations %>%
  pivot_longer(
    cols = 4:1003, # Adjust if your actual columns are named differently
    names_to = "measurement",
    values_to = "value"
  )

df_long$C40 <- factor(df_long$C40, levels = c(0, 1), labels = c("Non-C40 Cities (21)", "C40 Cities (51)"))
histogram1 <- ggplot(df_long, aes(x = value, fill = factor(year))) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  facet_wrap(~C40, nrow = 1, scales = "free") +
  labs(
    fill = "Year",
    x = bquote(bold("CH"[4] ~ " Emissions Realizations [Tg/y]")),
    y = "Count"
    #title = bquote(bold("CH"[4]~" Emissions Realizations"))
  ) +
  scale_fill_manual(values = regional_palette) +
  theme_linedraw() + plot_theme_text +
  theme(panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave("SI_figure_histograms.png", plot = histogram1, device = "png", path = path_to_plots, width = 12, height = 4, dpi = 320, bg = "white")



# Difference Distribution
df_wide <- network_realizations %>%
  pivot_longer(
    cols = 4:1003,
    names_to = "measurement",
    values_to = "value"
  ) %>%
  select(year, C40, measurement, value) %>%
  pivot_wider(
    names_from = year,
    values_from = value,
    names_prefix = "value_"
  )

year_pairs <- list(
  c(2023, 2019),
  c(2023, 2020),
  c(2023, 2022),
  c(2022, 2021),
  c(2021, 2020),
  c(2020, 2019)
)
n_iter <- 5000
mc_sample_group <- function(df_group, year_pairs, n_iter) {
  result <- list()
  for (pair in year_pairs) {
    y <- paste0("value_", unlist(pair))

    v1 <- df_group[[y[1]]]
    v2 <- df_group[[y[2]]]

    mc_v1 <- sample(v1, n_iter, replace = TRUE)
    mc_v2 <- sample(v2, n_iter, replace = TRUE)
    result[[paste0("diff_", pair[1], "_", pair[2])]] <- mc_v1 - mc_v2
  }
  result
}

df_mc_C40 <- df_wide %>%
  group_by(C40) %>%
  group_modify(~{
    mc_results <- mc_sample_group(.x, year_pairs, n_iter = 5000)
    tibble(
      diff_pair = names(mc_results),
      mc_values = mc_results
    )
  }) %>%
  unnest(mc_values)

# tidy up for plotting
df_mc_C40$C40_name <- factor(df_mc_C40$C40, levels = c(0, 1), labels = c("Non-C40 Cities (21)", "C40 Cities (51)"))
df_mc_C40$year_group_name <- factor(df_mc_C40$diff_pair,
                                      levels = c("diff_2023_2019", "diff_2023_2020", "diff_2023_2022", "diff_2022_2021", "diff_2021_2020", "diff_2020_2019"),
                                      labels = c("2023 - 2019", "2023 - 2020", "2023 - 2022", "2022 - 2021", "2021 - 2020", "2020 - 2019"))
df_mc_C40 <- df_mc_C40 %>%
                mutate(section = ifelse(year_group_name %in% c("2023 - 2022", "2022 - 2021", "2021 - 2020", "2020 - 2019"), "Individual Year", "Multiple Years"))
df_mc_C40$section <- factor(df_mc_C40$section, levels = c("Individual Year", "Multiple Years"))

ci_95_C40 <- df_mc_C40 %>%
                group_by(section, C40, C40_name, diff_pair, year_group_name) %>%
                summarize(mean = mean(mc_values),
                          lower_90 = quantile(mc_values, 0.05),
                          upper_90 = quantile(mc_values, 0.95),
                          lower_95 = quantile(mc_values, 0.025),
                          upper_95 = quantile(mc_values, 0.975),
                          ci_90 = ifelse((lower_90 < 0) & (upper_90 > 0), "contains_0", "no"),
                          ci_95 = ifelse((lower_95 < 0) & (upper_95 > 0), "contains_0", "no"))
ci_95_C40 %>% filter(section == "Multiple Years")

# look at change manually

ci_95_C40[(ci_95_C40$C40 == 0) & (ci_95_C40$diff_pair == "diff_2023_2020"), "mean"]

histogram_2 <- ggplot(data = df_mc_C40, aes(x = mc_values, fill = year_group_name)) +
geom_histogram(position = "identity", alpha = 0.75, bins = 50) +
  facet_grid(C40_name ~ section, scales = "free_y") +
  labs(
    fill = "Subtraction of Years",
    y = "Count",
    x = bquote(bold("Difference in CH"[4] ~ " Emissions Realizations"))
  ) +
  scale_fill_manual(values = regional_palette) +
  theme_linedraw() + plot_theme_text + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave("SI_difference_histograms.png", plot = histogram_2, device = "png", path = path_to_plots, width = 16, height = 6, dpi = 320, bg = "white")

# could also look at this by region, or for 72 cities in total