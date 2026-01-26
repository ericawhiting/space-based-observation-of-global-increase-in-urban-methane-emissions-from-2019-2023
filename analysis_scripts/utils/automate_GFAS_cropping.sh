#!/bin/bash

echo "Cropping GFAS"

#already converted GFAS grib to nc using cdo
#cdo -f nc copy GFAS_CO_2019_2023.grib GFAS_CO_2019_2023.nc

# The path to the CSV file
csv_file="../Data/city_information/city_bounds_larger.csv" # larger than domain
# PI value to be used in the bc calculations for converting degrees to radians
PI=$(echo "4*a(1)" | bc -l)
# Skip the header line using tail -n +2 to process from the second line
tail -n +2 "$csv_file" | while IFS=',' read -r city lat lon degrees_latitude; do
    city=$(echo "$city" | tr -d '"') # Remove quotes from each variable
    echo "Cropping $city"

    # Calculating the latitude adjustment factor lat_adj
    lat_rad=$(echo "$lat * $PI / 180" | bc -l) # Convert latitude from degrees to radians
    lat_adj=$(echo "1 / c($lat_rad)" | bc -l)  # Calculate adjustment factor

    # Adjust the degrees for longitude using lat_adj
    degrees_longitude=$(echo "$degrees_latitude * $lat_adj" | bc -l)
    # Now calculate the bounding box with the adjusted longitude
    min_lon=$(echo "$lon - $degrees_longitude" | bc -l) # Adjusting the longitude bounds
    max_lon=$(echo "$lon + $degrees_longitude" | bc -l)
    min_lat=$(echo "$lat - $degrees_latitude" | bc -l)  # Original latitude bounds
    max_lat=$(echo "$lat + $degrees_latitude" | bc -l)

    # Using cdo to crop the global file to the city's bounding box
    cdo sellonlatbox,$min_lon,$max_lon,$min_lat,$max_lat ../Data/GFAS/GFAS_CO_2018_2024.nc "../Data/GFAS/city_GFAS/${city}_cropped_GFAS.nc"

    # check exit status of cdo to ensure the operation was successful
    if [ $? -eq 0 ]; then
        echo "Cropped global data for $city saved to ${city}_cropped_GFAS.nc"
    else
        echo "Failed to crop global data for $city."
    fi
done