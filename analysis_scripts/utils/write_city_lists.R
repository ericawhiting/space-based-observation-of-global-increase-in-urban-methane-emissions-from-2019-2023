# you will need to adjust the paths and matching. Our files were cropped closer to the city level to reduce storage needs with a city domain code added to the file name and stored by year
# -------------------------------
# Load Libraries & Set Paths
# -------------------------------
library(stringr)
library(foreach)
library(doParallel)

path_to_utils <- "./utils/" # change your path here
source(paste0(path_to_utils, "urban_boundaries.R"))

path_to_data <- "../Data/" # change your path here

# -------------------------------
# Define the root directory where your data is stored
root_directory <- "TROPOMI_city_subsets/"

# Create a list of cities
cities <- read.csv(paste0(path_to_data, "city_information/city_codes.csv"), header = FALSE)


# remove parallelization call and change %dopar% to %do% if needed
registerDoParallel(16)

foreach(city = cities[, 1]) %dopar% {
    # Create a list to store the data for the current city
    domain_code <- return_city_boundary_flags(city)
    temp_ch4_files <- c()
    temp_ch4_orbits <- c()
    temp_co_files <- c()
    temp_co_orbits <- c()

    # Iterate through the directories [Your file structure may look different!!]
    for (year_directory in c("2019", "2020", "2021", "2022", "2023")) {
        year_path <- paste0(path_to_data, root_directory, year_directory, "/")
        for (file in list.files(year_path, pattern = sprintf(".*__CH4__.*-%s.nc", domain_code))) {
            orbit <- str_match(file, "_(\\d{5})_")[2]
            temp_ch4_files <- append(temp_ch4_files, paste0(root_directory, year_directory, "/", file))
            temp_ch4_orbits <- append(temp_ch4_orbits, orbit)
        }
        for (file in list.files(year_path, pattern = sprintf(".*__CO__.*-%s.nc", domain_code))) {
            orbit <- str_match(file, "_(\\d{5})_")[2]
            temp_co_files <- append(temp_co_files, paste0(root_directory, year_directory, "/", file))
            temp_co_orbits <- append(temp_co_orbits, orbit)
        }
    }
    # merge co and ch4
    ch4_files_df <- data.frame(ch4 = temp_ch4_files, orbit = temp_ch4_orbits)
    co_files_df <- data.frame(co = temp_co_files, orbit = temp_co_orbits)
    merged_df <- merge(ch4_files_df, co_files_df, how = "inner", on = orbit)
    merged_df <- subset(merged_df, select = -c(orbit))
    write.table(merged_df, paste0(path_to_data, "/city_information/city_list_of_tropomi_files/", city, "_files.csv"), row.names = FALSE, sep = ",")
}
