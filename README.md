# Code supporting Space-based observation of global increase in urban methane emissions from 2019-2023
### authors of code: Erica Whiting, Genevieve Plant
### contact: erwh@umich.edu
### publication doi: 
### code doi (please see Zenodo for updates to version): 10.5281/zenodo.17582608

### overview
These scripts and files were used to read TROPOMI CH4 and CO measurement data, calculate enhancement ratios for selected urban domains, and calculate annual CH4 emissions as is described in the publication "Space-based observation of global increase in urban methane emissions from 2019-2023" 

### code citation
Whiting, E. and Plant, G. 2026, Space-based observation of global increase in urban methane emissions from 2019-2023 Code, v1.0, Zenodo, doi:10.5281/zenodo.17582608

### summary of files provided
- analysis_scripts/
    - utils/
        - <u>automate_GFAS_cropping.sh</u> - crop GFAS wildfire co emissions netcdf4 data to area surrounding city to filter fire events
        - <u>calc_dry_column_mixing_ratio.R</u> - calculate dry column mixing ratio from TROPOMI CO total column 
        - <u>crop_pop_density.R</u> - crop population density down to use when finding city center in find_city_center_from_pop_dens in urban_boundaries.R (can skip this by using city_center_from_pop_density.csv instead and adjusting utils/urban_boundaries/find_city_center_from_pop_dens)
        - <u>find_box.R</u> - find ideal boxsize for each city (delta_lat and delta_lon are from center of box to edge) (can skip by using Data/city_information/ideal_boxsize.csv)
        - <u>load_ch4_emissions.R</u> - functions to load EDGARv8.0 and EDGARv2024 CH4 emissions for each city at monthly, seasonal, or annual level, uses some functions from load_co_emissions.R
        - <u>load_co_emissions.R</u> - functions to load EDGARv8.1 CO emissions for each city at monthly, seasonal, or annual level
        - <u>orbits_to_filter.R</u> - TROPOMI orbits to filter out when VIIRS was known to be down
        - <u>plot_overpasses.R</u> - plot TROPOMI CH4 and CO overpasses for city, used in TROPOMI_ch4co_city_analysis.R
        - <u>split_dates.R</u> - date and time related utility functions
        - <u>urban_boundaries.R</u> - all utilities related city location, size, characteristics, uses files in Data/city_information
        - <u>write_city_lists.R</u> - makes two column list of TROPOMI CH4 and CO file names per city aligned by orbit number
    - <u>aggregate_EDGAR_ch4_emissions.R</u> - create realizations for city CH4 emissions
    - <u>aggregate_EDGAR_co_emissons.R</u> - create realizations for city CO emissions based on uncertainty of choice (full magnitude, change in time, no uncertainty)
    - <u>TROPOMI_ch4co_annual_mean_and_uncertainties.R</u>- combines ch4co realizations and co emissions realizations based on uncertainty selected to calculate annual urban CH4 emissions for each city and emissions totals at the regional and total scale from 2019 to 2023
    - <u>TROPOMI_ch4co_city_analysis.R</u> - open daily TROPOMI CH4 and CO measurements to calculate daily enhancement ratios for each urban domain after filtering
    - <u>TROPOMI_ch4co_plots_and_tables.R</u> - used to create SI table, figures in main text, and calculations
- Data/
    - city_information/
        - <u>city_bounds_larger.csv</u> - used to crop TROPOMI files when downloading to a box larger than the urban domain, expanded for areas with multiple cities and for cropping GFAS data around city
        - <u>city_center_from_pop_density.csv</u> - has lat lon location of city center for each city based on the peak in population density
        - <u>city_codes.csv</u> - the city codes to run for the analysis scripts
        - <u>city_information.csv</u> - compiled city information (name, country, C40 region, lat, lon, latlon, coastal flag, etc)
        - <u>ideal_boxsize.csv</u> - the box size (before coastal adjustment) used in this work
    - results/
        - <u>allcities_over_time_ch4_emissions.csv</u> - summed annual CH4 emissions of 72 cities for which we have obs from 2019-2023 and 95% confidence intervals
        - <u>annual_sum_ch4_emissions_fulluncert.csv</u> - summed annual CH4 emissions of all cities observed in each year (# cities varies and which cities have obs vary by year) and 95% confidence intervals
        - <u>city_ch4co_and_emisisons.csv</u> - annual enhancement ratio (ch4co) 95% confidence intervals, CH4 emissions with differenct uncertainties based on CO emissions (full from EDGAR CO uncertainty magnitude and sron which describes confidence in potential co emissions change over time that isn't already accounted for in CO inventory) and 95% confidence intervals for each city for each year with region and number of observations
        - <u>network_ch4_emissions.csv</u> - summed annual CH4 emissions for C40 network vs non-C40 cities that we have observations over and 95% confidence intervals
        - <u>regional_ch4_emissions.csv</u> - summed annual CH4 emissions for C40 region and 95% confidence intervals



### supporting data (see data availability statement of publication for access)
- TROPOMI L2 CH4 and CO data (2019-2023)
- EDGAR monthly gridded CO emissions data (v8.1 2019-2022)
- EDGAR monthly gridded CH4 emissions data (v8.0 2019-2022, v2024 2019-2023)
- CAMS GFAS Wildfire flux of carbon monoxide (CO) daily averages(2019-2023)
- GHSL data for pop density (can skip, see notes above)

in the file structure below, an asterisk (*) indicates where the scripts will expect these data to be before running

### expected structure
```bash
space-based-observation-of-global-increase-in-urban-methane-emissions-from-2019-2023/
├── analysis_scripts/
│   ├── utils/
│   │   ├── automate_GFAS_cropping.sh
│   │   ├── calc_dry_column_mixing_ratio.R
│   │   ├── crop_pop_density.R
│   │   ├── find_box.R
│   │   ├── load_ch4_emissions.R
│   │   ├── load_co_emissions.R
│   │   ├── orbits_to_filter.R
│   │   ├── plot_overpasses.R
│   │   ├── split_dates.R
│   │   ├── urban_boundaries.R
│   │   └── write_city_lists.R
│   ├── aggregate_EDGAR_ch4_emissions.R
│   ├── aggregate_EDGAR_co_emissons.R
│   ├── TROPOMI_ch4co_annual_mean_and_uncertainties.R
│   ├── TROPOMI_ch4co_city_analysis.R
│   └── TROPOMI_ch4co_plots_and_tables.R
├── Data/
│   ├── ch4co/
│   ├── city_information/
│   │   ├── pop_density_cropped/
│   │   ├── city_bounds_larger.csv
│   │   ├── city_center_from_pop_density.csv
│   │   ├── city_codes.csv
│   │   ├── city_information.csv
│   │   └── ideal_boxsize.csv
│   ├── EDGAR/
│   │   ├── ch4/
│   │   │   └── monthly_sector_specific_gridmaps/
│   │   │   │   ├── 2024/*
│   │   │   │   └── v8/*
│   │   └── co/
│   │   │   └── monthly_sector_specific_gridmaps/
│   │   │   │   └── v8.1/*
│   ├── EDGAR_ch4_emissions/
│   │   └── annual/
│   ├── EDGAR_co_emissions_realizations/
│   │   ├── full_uncertainty/
│   │   │   │   └── annual/
│   │   ├── no_uncertainty/
│   │   │   │   └── annual/
│   │   ├── SRON_uncertainty/
│   │   │   │   └── annual/
│   ├── GFAS/
│   │   └── global_GFAS_data/*
│   │   └── city_GFAS/
│   ├── GHSL/* or can skip (see note above)
│   ├── HTAPv3/
│   │   └── annual/* or can skip (see note above)
│   ├── results/
│   │   ├── realizations/
│   │   │   ├── annual_ch4_realizations_SRONuncert.csv
│   │   │   ├── all_city_realizations.csv
│   │   │   ├── network_realizations.csv
│   │   │   └── regional_realizations.csv
│   │   ├── allcities_over_time_ch4_emissions.csv
│   │   ├── annual_sum_ch4_emissions_fulluncert.csv
│   │   ├── city_ch4co_and_emisisons.csv
│   │   ├── network_ch4_emissions.csv
│   │   └── regional_ch4_emissions.csv
│   └── TROPOMI_city_subsets/
│   │   ├── 2019/*
│   │   ├── 2020/*
│   │   ├── 2021/*
│   │   ├── 2022/*
│   │   └── 2023/*
├── plots/
│   ├── results/
│   ├── tropomi_overpasses/
│   │   └──temp_plots/
└── README.md
```