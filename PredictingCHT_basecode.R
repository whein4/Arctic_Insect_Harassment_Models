##################################################################################
##################################################################################
##################################################################################


# This R code is designed to predict insect harassment
# using the CHT insect indices introduced in Hein et al. 2025.
# Please refer to the documentation provided on Github
# https://github.com/whein4/Arctic_Insect_Harassment_Models/blob/main/README.md

# Contact: william.hein@mail.mcgill.ca


##################################################################################
##################################################################################
##################################################################################


# In this code we start from scratch and end with predicted mosquito and oestrid fly harassment multidimensional rasters
# Step 1: Downloading ERA5 rasters 
# Step 2: Manipulate ERA5 rasters
# Step 3: Download DSM and Calculate TPI
# Step 4: Standardize rasters
# Step 5: Create a time of day raster
# Step 6: Predict Using CHT model
# Step 7: Visualize the results


# Install necessary libraries
install.packages("ecmwfr")
devtools::install_github("https://github.com/ErikKusch/KrigR")
install.packages("lubridate")
install.packages("terra")
install.packages("devtools")
install.packages("zoo")
install.packages("mgcv")
install.packages("httr")
install.packages("jsonlite")
install.packages("raster")
install.packages("progressr")
install.packages("bgeva")
install.packages("ggplot2")
install.packages("shiny")
install.packages("leaflet")


# Load necessary libraries
library(ecmwfr)
library(KrigR)
library(lubridate)
library(terra)
library(devtools)
library(zoo)
library(mgcv)
library(bgeva)
library(httr)
library(jsonlite)
library(raster)
library(progressr) 
library(ggplot2)
library(shiny)
library(leaflet)


# Set your working and output directory
setwd("C:/your/working/directory")
output_dir <- "C:/your/output/directory"


##################################################################
              # Step 1 Downloading ERA5 rasters 
##################################################################


################################
      # Define parameters
###############################
# Set your ERA5 API key (https://cds.climate.copernicus.eu/profile)
# The one provided below will not work, you need to create your own
wf_set_key(key = "a966427f-5ff3-4s72-q737-65d4c92fee31")

# Set your study year 
# ERA5 allows for analysis back as far as 1950, and up to about a month before the present day
year <- "2018"

# Set your study area
# Note that the CHT insect indices were developed on the Yukon North Slope
# They may not predict as accurately in significantly different ecosystems
area <- c(70, -141, 68.5, -138)

# Define your ERA5 search period and data manipulation parameters
# Note you should only be editing the year of these parameters, unless you have good reason to do so
months <- c("05", "06", "07", "08")
snowmonth <- c("05", "06", "07")
precipmonth <- c("04", "05", "06", "07", "08")
days <- sprintf("%02d", 1:31)
times <- sprintf("%02d:00", 0:23)
start_date <- as.POSIXct("2018-05-01 00:00", tz = "UTC")
index_start_date <- as.POSIXct("2018-06-01 12:00", tz = "UTC")
index_end_date <- as.POSIXct("2018-08-31 12:00", tz = "UTC")
index_start_date_hourly <- as.POSIXct("2018-06-01 00:00", tz = "UTC")
index_end_date_hourly <- as.POSIXct("2018-08-31 23:00", tz = "UTC")
end_date <- as.POSIXct("2018-08-31 23:00", tz = "UTC")
snow_melt_end <- as.POSIXct("2018-07-31 23:00", tz = "UTC")
snow_start_date <- as.POSIXct(paste0(year, "-05-01 00:00"), tz = "UTC")
precipitation_start <- as.POSIXct("2018-04-01 00:00", tz = "UTC")



################
# Download temperature raster
###############
# Define request parameters
request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "2m_temperature",
    year = year,
    month = months,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_2m_temperature_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,  
    path = output_dir         
  )
  
  #unzip the file
  nc_file <- paste0("era5_2m_temperature_land_", year, ".nc")
  zip_file <- paste0("era5_2m_temperature_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  temp_raw <- rast(nc_file)
  
  # Define the time sequence
  time_points <- seq(from = start_date, to = end_date, by = "hour")
  terra::time(temp_raw) <- time_points


# Repeat to download code for other covariates:
################
# Download U component of wind raster
###############
request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "10m_u_component_of_wind",
    year = year,
    month = months,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_10m_u_wind_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,   
    path = "."         
  )
  
  #unzip the file
  nc_file <- paste0("era5_10m_u_wind_land_", year, ".nc")
  zip_file <- paste0("era5_10m_u_wind_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  wind_u_raw <- rast(nc_file)
  
  
  # Define the time sequence
  time_points <- seq(from = start_date, to = end_date, by = "hour")
  terra::time(wind_u_raw) <- time_points

################
# Download V component of wind raster
###############
  request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "10m_v_component_of_wind",
    year = year,
    month = months,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_10m_v_wind_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,  
    path = "."         
  )
  
  #unzip the file
  nc_file <- paste0("era5_10m_v_wind_land_", year, ".nc")
  zip_file <- paste0("era5_10m_v_wind_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  wind_v_raw <- rast(nc_file)
  
  
  # Define the time sequence
  time_points <- seq(from = start_date, to = end_date, by = "hour")
  terra::time(wind_v_raw) <- time_points

################
# Download soil moisture raster
###############
request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "volumetric_soil_water_layer_1",
    year = year,
    month = months,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_soil_moisture_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,   
    path = "."         
  )
  
  #unzip the file
  nc_file <- paste0("era5_soil_moisture_land_", year, ".nc")
  zip_file <- paste0("era5_soil_moisture_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  soil_raw <- rast(nc_file)
  
  
  # Define the time sequence
  time_points <- seq(from = start_date, to = end_date, by = "hour")
  terra::time(soil_raw) <- time_points

################
# Download snow depth raster
###############
request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "snow_depth",
    year = year,
    month = snowmonth,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_snow_depth_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,   
    path = "."        
  )
  
  #unzip the file
  nc_file <- paste0("era5_snow_depth_land_", year, ".nc")
  zip_file <- paste0("era5_snow_depth_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  snow_raw <- rast(nc_file)
  
  # Define the time sequence
  time_points <- seq(from = snow_start_date, to = snow_melt_end, by = "hour")
  terra::time(snow_raw) <- time_points


################
# Download precipitation raster
###############
request <- list(
    dataset_short_name = "reanalysis-era5-land",
    product_type = "reanalysis",
    variable = "total_precipitation",
    year = year,
    month = precipmonth,
    day = days,  
    time = times, 
    area = area, 
    format = "netcdf",
    target = paste0("era5_total_precipitation_land_", year, ".nc")
  )
  
  file_path <- wf_request(
    request = request,
    transfer = TRUE,   
    path = "."         
  )
  
  #unzip the file
  nc_file <- paste0("era5_total_precipitation_land_", year, ".nc")
  zip_file <- paste0("era5_total_precipitation_land_", year,".zip")
  
  if (!file.exists(nc_file)) {
    # Unzip and get the list of extracted files
    files <- unzip(zip_file, exdir = ".")
    # Rename the first extracted file to your desired name if needed
    if (basename(files[1]) != nc_file) {
      file.rename(files[1], nc_file)
    }
  }
  
  # Load the NetCDF file
  precip_raw <- rast(nc_file)
  
  # Define the time sequence
  time_points <- seq(from = precipitation_start, to = end_date, by = "hour")
  terra::time(precip_raw) <- time_points
  
  
  


# Check the rasters loaded correctly
print(temp_raw)
print(snow_raw)
print(wind_v_raw)
print(wind_u_raw)
print(precip_raw)
print(soil_raw)










#######################################################################################
                      # Step 2 Manipulate ERA5 rasters
######################################################################################
# We need to manipulate the raw ERA5 rasters to become the covariates used in the CHT insect indices

######################################
# Celsius to kelvin
######################################
# Convert temperature from Kelvin to Celsius
temp_celsius <- temp_raw - 273.15

######################################
# Temp to cumulative growing degree day
######################################
# We need to calculate cumulative growing degree day using temperature
# Aggregate the hourly raster to daily minimum
temp_min_daily <- tapp(temp_celsius, "days", fun = min)

# Aggregate the hourly raster to daily maximum
temp_max_daily <- tapp(temp_celsius, "days", fun = max)

# Compute the daily average temperature
gdd_raster <- (temp_min_daily + temp_max_daily) / 2

# Base temperature = 0 C*
# Set all negative values to zero
gdd_raster[gdd_raster < 0] <- 0

# Calculate cumulative GDD
gdd_cumulative <- app(gdd_raster, function(x) {
  cumsum(x)  # Compute cumulative sum across layers for each pixel
})

# Inspect the resulting cumulative GDD raster
print(gdd_cumulative)

######################################
# Calculate wind speed
######################################
# To get wind speed we need to calculate magnitude of the U and V components of wind
# Calculate the square root of the sum of squares
wind_speed <- sqrt(wind_u_raw^2 + wind_v_raw^2)

# Inspect the resulting wind speed raster
print(wind_speed)

######################################
# Calculate average summer soil moisture
######################################
# Calculate the average across all layers
soil_avg <- mean(soil_raw, na.rm = TRUE)

# Inspect the resulting raster
print(soil_avg)

######################################
# Calculate days since snow melt from snow depth 
######################################
# We used ERA5 snow depth >0.2m to define "snow melt"
# Aggregate the snow_raw raster to daily minima
snow_daily_min <- tapp(snow_raw, "days", fun = min)

# Find the layer index where snow depth < 0.2 for the first time
first_melt_layer <- app(snow_daily_min, function(x) {
  melt_index <- which(x < 0.2)[1]  # Find the first occurrence of snow depth < 0.2
  if (!is.na(melt_index)) {
    melt_index  # Return the layer index
  } else {
    NA  
  }
})

# For each pixel we should have layer number corresponding to the first day where snow depth < 0.2
print(first_melt_layer)

# Now we need to convert that into a "days since snow melt" raster
# Define the start and end dates (the year does not matter here, we will fix that down the road)
start_date <- as.Date("2022-05-01")
end_date <- as.Date("2022-08-31")
all_dates <- seq(start_date, end_date, by = "day")
num_days <- length(all_dates)

# Create an empty raster stack with the same dimensions and resolution as snow_daily_min
index_stack <- rast(snow_daily_min, nlyr = num_days)

# Assign sequential indices (1, 2, ...) to the raster stack
values(index_stack) <- rep(1:num_days, each = ncell(index_stack))

# Assign date-based names to the layers
index_stack <- setNames(index_stack, paste0("d_", all_dates))

# Expand first_melt_layer to match the number of layers in the index_stack
expanded_first_melt <- rast(lapply(1:num_days, function(i) first_melt_layer))

# Ensure the expanded raster has the same layer names as the index stack
expanded_first_melt <- setNames(expanded_first_melt, paste0("d_", all_dates))

# Subtract the expanded first melt raster from the index stack
days_since_snow_melt <- index_stack - expanded_first_melt

# Inspect the resulting "days since snow melt" raster
print(days_since_snow_melt)

######################################
#Calculate 2 week sum of precipitation
######################################
# Aggregate to daily totals
precip_daily <- tapp(precip_raw, "days", fun = sum)

# Extract raster values as a matrix
precip_values <- as.matrix(precip_daily)

# Apply rolling sum for each row (pixel-wise operation)
rolling_sums <- t(apply(precip_values, 1, function(x) {
  rollapply(x, width = 15, FUN = sum, align = "right", fill = NA)
}))

# Create a new raster with the rolling sum values
precip_2wk <- setValues(precip_daily, rolling_sums)

# Inspect the result
print(precip_2wk)










#######################################################################################
                          # Step 3 Download DSM and Calculate TPI
######################################################################################

###########################
# If you have already calculated and saved TPI - skip to "Load TPI"
###########################

# If this is your first calculation, you will need to download DSM tiles
# Download ALOS DSM tiles within your study area from https://www.eorc.jaxa.jp/ALOS/en/dataset/aw3d30/aw3d30_e.htm
# Save all data files to a single file

####################
# Merge DSM tiles
###################
# Set the path to the folder containing the DSM tiles
folder_path <- "C:/your/DSM/file/location"

# List all DSM files that end in 'DSM'
dsm_files <- list.files(folder_path, pattern = "DSM", full.names = TRUE)

# Check the list of files
print(dsm_files)

# Read all DSM tiles into a list of SpatRasters
dsm_rasters <- lapply(dsm_files, rast)

# Use the first tile as a reference
reference_raster <- dsm_rasters[[1]]

# Align all rasters while preserving extent
dsm_rasters_aligned <- lapply(dsm_rasters, function(r) {
  project(r, reference_raster, method = "bilinear", align = TRUE)
})

# Merge the rasters into one, aligning them based on their extents
merged_dsm <- do.call(mosaic, c(dsm_rasters_aligned, fun = "mean"))

#visual check
plot(merged_dsm)  


################
# Resample merged DSM
################
# As described in Hein et al. (2025) we used a coarser version (30m to 1 km) of the DSM to calculate TPI
# Reproject DSM raster to UTM (assuming UTM zone 7N is appropriate for your area)
utm_crs <- "+proj=utm +zone=7 +datum=WGS84"
dsm_utm <- project(merged_dsm, utm_crs)

# Resample DSM to 1 km resolution
dsm_1km <- resample(dsm_utm, rast( ext = ext(dsm_utm), resolution = 1000), method = "bilinear")

################
# Calculate TPI
################
# Calculate the TPI over a 5 km neighborhood (5000 meters)
# Define the radius in pixels for a 5km neighborhood (1km resolution, so 5 pixels)
radius_pixels <- 5

# Create a focal window matrix with a size corresponding to a 5 km radius
focal_window <- matrix(1, nrow = radius_pixels * 2 + 1, ncol = radius_pixels * 2 + 1)

# Calculate the mean elevation within a 5 km radius for each pixel
mean_elevation_raster <- focal(dsm_1km, w = focal_window, fun = mean, na.rm = TRUE, pad = TRUE)

# Calculate the TPI by subtracting the DSM from the mean elevation raster
tpi_raster <- mean_elevation_raster - dsm_1km

# Plot the TPI raster
plot(tpi_raster, main = "Topographic Position Index (TPI)")

# Save TPI Raster for future use
writeRaster(tpi_raster, "TPI_5km_raster.tif", overwrite = TRUE)






################
# Load TPI
################
# Skip to this once TPI has been calculate once to save time
# Load TPI raster
tpi_raster <- rast("TPI_5km_raster.tif")
# Plot 5km TPI raster
plot(tpi_raster)








#########################################################################################################
                            # Step 4 Standardize rasters
#########################################################################################################
# In order to predict using the CHT insect indice's GAM models, all rasters must be formatted identically

######################################
# subset 2 week sum of precip
######################################
#We do not need the month of May from the 2 week sum of precip raster
precip_2wk <- subset(precip_2wk, 31:153)

######################################
# standardize the coordinate reference system
######################################
# Define the common CRS (WGS84)
common_crs <- "+proj=longlat +datum=WGS84 +no_defs"

# Reproject all rasters to the common CRS
temp_celsius_proj <-  terra::project(temp_celsius, common_crs)
wind_speed_proj <- project(wind_speed, common_crs)
precip_proj <- project(precip_2wk, common_crs)
gdd_proj <- project(gdd_cumulative, common_crs)
dssm_proj <- project(days_since_snow_melt, common_crs)
soil_proj <- project(soil_avg, common_crs)
tpi_proj <- project(tpi_raster, common_crs)




######################################
# standardize the number of layers
######################################
# We need all rasters to have the same number of layers
# This is where we can subset the rasters to save processing time and effort
# Select either option 1 or option 2, depending on your needs and computing capabilities

# Option 1: Daily metric
# This code will subset the rasters to one hour per day, essentially becoming daily rasters, rather then hourly
# This is recommended if your computer is limited in processing power (</= 16 Gb of RAM)

#######################
# Calculating daily metric of insect harassment
############################
# We only need 12:00 for hourly rasters (roughly peak daily harassment)
# Note you need to consider the time zone you are predicting on, the standard here is Yukon time (PST: UTC-7)
# If we want PST = 12:00, we need to start at UTC = 19:00
# If Layer 1 = 05/01 00:00 UTC, then layer 745 = 06/01 00:00 UTC
# That means 06/01 19:00 UTC = 764th layer
# If we select the 764th layer and every 24th layer after that, we will get a daily harassment raster at 12:00 PST
layer_indices <- seq(764, 2948, by = 24)

# Subset each raster individually
temp_celsius_proj_daily <- subset(temp_celsius_proj, layer_indices)
wind_speed_proj_daily <- subset(wind_speed_proj, layer_indices)

# since precipitation, GDD, and DSSM are daily rasters, we need to start at the 32nd layer (06/01)
precip_proj_daily <- subset(precip_proj, 32:123)
gdd_proj_daily <- subset(gdd_proj, 32:123)
dssm_proj_daily <- subset(dssm_proj, 32:123)

# Function to expand single layer rasters
expand_single_layer <- function(single_layer, n_layers) {
  # Duplicate the single layer n_layers times
  expanded_raster <- terra::rast(replicate(n_layers, single_layer))
  
  return(expanded_raster)
}

# Expand soil and tpi rasters to 92 layers (number of layers in the other rasters)
soil_proj_daily <- expand_single_layer(soil_proj, 92)
tpi_proj_daily <- expand_single_layer(tpi_proj, 92)

# Check to ensure layer numbers are consistent 
print(temp_celsius_proj_daily)
print(wind_speed_proj_daily)
print(precip_proj_daily)
print(gdd_proj_daily)
print(dssm_proj_daily)
print(soil_proj_daily)
print(tpi_proj_daily)
# End option 1




# Option 2: Hourly metric
# The following code will retain the temporal extent of the original ERA5 rasters (hourly)
# This is recommended if your computer has significant processing power (</= 32 Gb of RAM)

#############
# Calculating hourly metric of insect harassment
#############
# We only need the 06/01 - 08/31 period. 
# If Layer 1 = 05/01 00:00 UTC, then layer 745 = 06/01 00:00 UTC
# Define the sequence of layers to subset
layer_indices <- seq(745, 2952)

# Subset each raster individually
temp_celsius_proj_daily <- subset(temp_celsius_proj, layer_indices)
wind_speed_proj_daily <- subset(wind_speed_proj, layer_indices)

# Since precipitation, GDD, and DSSM are daily rasters, we need to start at the 32nd layer (06/01)
precip_proj_daily <- subset(precip_proj, 32:123)
gdd_proj_daily <- subset(gdd_proj, 32:123)
dssm_proj_daily <- subset(dssm_proj, 32:123)


# We need to expand precipitation, GDD, and DSSM to match the number of layers from above by expanding each daily layer by 24
expand_to_hourly <- function(raster, repetitions) {
  n_layers <- terra::nlyr(raster)
  expanded_layers <- lapply(seq_len(n_layers), function(i) {
    replicate(repetitions, raster[[i]])
  })
  expanded_raster <- terra::rast(do.call(c, expanded_layers))
  return(expanded_raster)
}

# Number of repetitions per layer
repetitions_per_layer <- 24

# Expand each raster to 2208 layers
precip_proj_daily <- expand_to_hourly(precip_proj_daily, repetitions_per_layer)
gdd_proj_daily <- expand_to_hourly(gdd_proj_daily, repetitions_per_layer)
dssm_proj_daily <- expand_to_hourly(dssm_proj_daily, repetitions_per_layer)


# Function to expand single layer rasters to match the number of layers to the other rasters
expand_single_layer <- function(single_layer, n_layers) {
  # Duplicate the single layer n_layers times
  expanded_raster <- terra::rast(replicate(n_layers, single_layer))
  
  return(expanded_raster)
}

# Expand soil and tpi rasters to 2208 layers (number of layers in the other rasters)
soil_proj_daily <- expand_single_layer(soil_proj, 2208)
tpi_proj_daily <- expand_single_layer(tpi_proj, 2208)

# Check to ensure layer numbers are consistent 
print(temp_celsius_proj_daily)
print(wind_speed_proj_daily)
print(precip_proj_daily)
print(gdd_proj_daily)
print(dssm_proj_daily)
print(soil_proj_daily)
print(tpi_proj_daily)
#End option 2



#############################
# Standardize Extent
#############################
# We need to ensure all rasters have the same extent, they should be close already but close is not good enough
# List rasters
raster_list <- list(temp_celsius_proj_daily, wind_speed_proj_daily, soil_proj_daily, precip_proj_daily, 
                    dssm_proj_daily, tpi_proj_daily, gdd_proj_daily)

# Extract extents
extents <- lapply(raster_list, terra::ext)

# Find the minimum bounding box (overlapping extent)
min_extent <- Reduce(terra::intersect, extents)

# Verify the minimum extent
print(min_extent)

# Crop all rasters to the smallest extent
temp_celsius_proj_crop <- terra::crop(temp_celsius_proj_daily, min_extent)
wind_speed_proj_crop <- terra::crop(wind_speed_proj_daily, min_extent)
precip_proj_crop <- terra::crop(precip_proj_daily, min_extent)
gdd_proj_crop <- terra::crop(gdd_proj_daily, min_extent)
dssm_proj_crop <- terra::crop(dssm_proj_daily, min_extent)
soil_proj_crop <- terra::crop(soil_proj_daily, min_extent)
tpi_proj_crop <- terra::crop(tpi_proj_daily, min_extent)

# Check extents after cropping
print(terra::ext(temp_celsius_proj_crop))
print(terra::ext(wind_speed_proj_crop))
print(terra::ext(precip_proj_crop))
print(terra::ext(gdd_proj_crop))
print(terra::ext(dssm_proj_crop))
print(terra::ext(soil_proj_crop))
print(terra::ext(tpi_proj_crop))

###################################
# Standardize resolutions
###################################
# Next we need to make sure the resolutions match
# Use the finest resolution raster (TPI) as the reference raster
reference_raster <- tpi_proj_crop

# Resample all rasters to align with the reference raster
temp_celsius_proj_res_crop <- terra::resample(temp_celsius_proj_crop, reference_raster)
wind_speed_proj_res_crop <- terra::resample(wind_speed_proj_crop, reference_raster)
soil_proj_res_crop <- terra::resample(soil_proj_crop, reference_raster)
precip_proj_res_crop <- terra::resample(precip_proj_crop, reference_raster)
dssm_proj_res_crop <- terra::resample(dssm_proj_crop, reference_raster)
gdd_proj_res_crop <- terra::resample(gdd_proj_crop, reference_raster)
tpi_proj_res_crop <- tpi_proj_crop 

# Check extents after resampling
print(terra::ext(temp_celsius_proj_res_crop))
print(terra::ext(wind_speed_proj_res_crop))
print(terra::ext(soil_proj_res_crop))
print(terra::ext(precip_proj_res_crop))
print(terra::ext(dssm_proj_res_crop))
print(terra::ext(gdd_proj_res_crop))
print(terra::ext(tpi_proj_res_crop))

# Check resolution
print(terra::res(temp_celsius_proj_res_crop))
print(terra::res(wind_speed_proj_res_crop))
print(terra::res(soil_proj_res_crop))
print(terra::res(precip_proj_res_crop))
print(terra::res(dssm_proj_res_crop))
print(terra::res(gdd_proj_res_crop))
print(terra::res(tpi_proj_res_crop))












#########################################################################################################
                                # Step 5 Create a time of day raster
#########################################################################################################
# We need to create a time of day raster to predict oestrid fly harassment
# We will need to create the raster based on how you subsetted the ERA5 covariates (daily vs hourly)
# Select either option 1 or option 2, depending on your selection when setting layer numbers

# Option 1: Daily time of day raster 
# If you selected option 1 for subsetting the ERA5 covariates (step 4)

#######################
# Creating a daily time of day raster 
############################
# We only need 12:00 for our time of day raster 
# Extract raster parameters from the reference raster
nrows <- nrow(temp_celsius_proj_res_crop)
ncols <- ncol(temp_celsius_proj_res_crop)
nlayers <- nlyr(temp_celsius_proj_res_crop)
extent_obj <- ext(temp_celsius_proj_res_crop)
crs_obj <- crs(temp_celsius_proj_res_crop)

# Create a new raster with the same dimensions
time_proj_res_crop <- rast(
  nrows = nrows, 
  ncols = ncols, 
  nlyrs = nlayers, 
  extent = extent_obj, 
  crs = crs_obj
)

# Set all raster values to 12.0
values(time_proj_res_crop) <- 12.0

# Inspect the resulting raster
print(time_proj_res_crop)

# Inspect the resulting raster
print(time_proj_res_crop)
# End option 1


# Option 2: Hourly time of day raster 
# If you selected option 2 for subsetting the ERA5 covariates (step 4)

#######################
# Creating a daily time of day raster
############################
# We need 0-23 ordered correctly based on your time zone
# Extract raster parameters from a reference raster
nrows <- nrow(temp_celsius_proj_res_crop)
ncols <- ncol(temp_celsius_proj_res_crop)
nlayers <- nlyr(temp_celsius_proj_res_crop)
extent_obj <- ext(temp_celsius_proj_res_crop)
crs_obj <- crs(temp_celsius_proj_res_crop)

# Create a vector of time values (0.00 to 23.00, repeating across layers)
# The first layer will be 00:00 UTC, change the starting time based on your time zone 
# For example, Yukon time is 7 hours behind UTC, 00:00 - 07:00 = first layer starts at 17:00
time_values <- rep(c(17:23, 0:16), length.out = nlayers)

# Create a matrix with dimensions matching a single layer
single_layer <- matrix(0, nrow = nrows, ncol = ncols)

# Create a new raster by repeating the single layer and assigning time values
time_proj_res_crop <- rast(
  nrows = nrows, 
  ncols = ncols, 
  nlyrs = nlayers, 
  extent = extent_obj, 
  crs = crs_obj
)

# Populate all layers at once
values(time_proj_res_crop) <- rep(time_values, each = nrows * ncols)

# Inspect the resulting raster
print(time_proj_res_crop)
# End option 2




######################################
# Standardize time 
######################################
# We need to standardize the time of all rasters now

# Option 1: Daily time of day raster 
# If you selected option 1 previously
# Assign daily time to the expanded raster
time_reference <- seq(index_start_date,
                      index_end_date,
                      by = "day")
#End option 1


#Option 2: Hourly time of day raster 
#If you selected option 2 previously
# Assign hourly time to the expanded raster
time_reference <- seq(index_start_date_hourly,
                      index_end_date_hourly,
                      by = "hour")
#End option 2


# Set the time reference to all covariate rasters
terra::time(temp_celsius_proj_res_crop) <- time_reference
terra::time(wind_speed_proj_res_crop) <- time_reference
terra::time(precip_proj_res_crop) <- time_reference
terra::time(gdd_proj_res_crop) <- time_reference
terra::time(dssm_proj_res_crop) <- time_reference
terra::time(tpi_proj_res_crop) <- time_reference
terra::time(soil_proj_res_crop) <- time_reference
terra::time(time_proj_res_crop) <- time_reference


# Check mosquito covariates
print(temp_celsius_proj_res_crop)
print(wind_speed_proj_res_crop)
print(soil_proj_res_crop)
print(precip_proj_res_crop)
print(dssm_proj_res_crop)
print(tpi_proj_res_crop)

# Check oestrid fly covariates
print(temp_celsius_proj_res_crop)
print(gdd_proj_res_crop)
print(time_proj_res_crop)


# Clean the environment to save space
objects_to_keep <- c(
  "temp_celsius_proj_res_crop", 
  "wind_speed_proj_res_crop", 
  "precip_proj_res_crop", 
  "gdd_proj_res_crop", 
  "dssm_proj_res_crop", 
  "soil_proj_res_crop", 
  "tpi_proj_res_crop",
  "time_proj_res_crop"
)

# Remove all objects except the ones in the list
rm(list = setdiff(ls(), objects_to_keep))

# Garbage collection to free memory
gc()














##################################################################################################################
                                            #Step 6 Predict Using CHT model
##################################################################################################################
# Load models
gam_oest <- readRDS("C:/your/directory/containing/GAMoestFINAL.rds")
gam_mosq <- readRDS("C:/your/directory/containing/GAMmosqFINAL.rds")

###################################
# Predict mosquito harassment
###################################
# Create a list of separate rasters with unique names for gam_mosq
predictors_mosq <- list(
  ERAtemp = temp_celsius_proj_res_crop,
  ERAwind = wind_speed_proj_res_crop,
  ERAsoil = soil_proj_res_crop,
  ERA5kTPI = tpi_proj_res_crop,
  ERA2wkprecip = precip_proj_res_crop,
  DSSM = dssm_proj_res_crop
)

# Combine predictors into a matrix (layer by layer)
predictor_values_mosq <- as.data.frame(cbind(
  ERAtemp = as.vector(terra::values(predictors_mosq$ERAtemp)),
  ERAwind = as.vector(terra::values(predictors_mosq$ERAwind)),
  ERAsoil = as.vector(terra::values(predictors_mosq$ERAsoil)),
  ERA5kTPI = as.vector(terra::values(predictors_mosq$ERA5kTPI)),
  ERA2wkprecip = as.vector(terra::values(predictors_mosq$ERA2wkprecip)),
  DSSM = as.vector(terra::values(predictors_mosq$DSSM))
))

# Predict using the model
predicted_values_mosq <- predict(gam_mosq, newdata = predictor_values_mosq, type = "response")
predicted_values_mosq_index <- predicted_values_mosq[, 1] * 0 + 
  predicted_values_mosq[, 2] * 1/3 + 
  predicted_values_mosq[, 3] * 2/3 + 
  predicted_values_mosq[, 4] * 1

# Create a new raster to store predictions
prediction_raster_mosq <- terra::rast(predictors_mosq$ERAtemp) 

# Populate blank  with NAs to be filled in with predictions
terra::values(prediction_raster_mosq) <- NA 

# Assign predicted values back to the raster (match indices)
all_indices <- seq_len(ncell(prediction_raster_mosq))
valid_indices <- which(!is.na(terra::values(predictors_mosq$ERAtemp)))
terra::values(prediction_raster_mosq)[valid_indices] <- predicted_values_mosq_index

# View the result
print(prediction_raster_mosq)
min(prediction_raster_mosq)
max(prediction_raster_mosq)

############################################################################
# Predict oestrid harassment
############################################################################

# Create a list of separate rasters with unique names
predictors_oest <- list(
  ERAtemp = temp_celsius_proj_res_crop,
  ERAgdd = gdd_proj_res_crop,
  TimeNumeric = time_proj_res_crop
)

# Combine predictors into a matrix (layer by layer)
predictor_values_oest <- as.data.frame(cbind(
  ERAtemp = as.vector(terra::values(predictors_oest$ERAtemp)),
  ERAgdd = as.vector(terra::values(predictors_oest$ERAgdd)),
  TimeNumeric = as.vector(terra::values(predictors_oest$TimeNumeric))
))

# Predict using the model
predicted_values_oest <- predict(gam_oest, newdata = predictor_values_oest, type = "response")

# Square root transform the oestrid predictions (See Hein et al. 2025)
predicted_values_oest_sq <- sqrt(predicted_values_oest)


# Create a new raster to store predictions
prediction_raster_oest <- terra::rast(predictors_oest$ERAtemp) 

# Populate blank  with NAs to be filled in with predictions
terra::values(prediction_raster_oest) <- NA 

# Assign predicted values back to the raster (match indices)
all_indices <- seq_len(ncell(prediction_raster_oest))
valid_indices <- which(!is.na(terra::values(predictors_oest$ERAtemp)))
terra::values(prediction_raster_oest)[valid_indices] <- predicted_values_oest_sq

# View  the result
print(prediction_raster_oest)
min(prediction_raster_oest)
max(prediction_raster_oest)


# Save the resulting predicted rasters
writeRaster(prediction_raster_oest, "oest_pred_2018.nc", overwrite = FALSE)
writeRaster(prediction_raster_mosq, "mosq_pred_2018.nc", overwrite = FALSE)














##################################################################################################################
                                        # Step 7 Visualize the results
##################################################################################################################
# Lets visualize the resulting raster

#################################
# Visualizing Harassment over time
#################################
# We can find the average harassment across all pixels for each layer (point in time) and plot it

# Function to calculate averages and plot
plot_raster_averages <- function(raster_data, title) {
  # Extract time information
  time_values <- time(raster_data)
  
  # Calculate average values for each layer
  layer_averages <- global(raster_data, fun = "mean", na.rm = TRUE)
  
  # Create a data frame with time and average values
  plot_data <- data.frame(
    Time = as.POSIXct(time_values),  # Convert time to POSIXct for ggplot
    AvgValue = layer_averages$mean   # Corresponding average values
  )
  
  # Plot the averages over time
  ggplot(plot_data, aes(x = Time, y = AvgValue)) +
    geom_line(color = "blue") +  # Line plot
    geom_point(color = "red") +  # Points for each time step
    labs(
      title = title,
      x = "Date",
      y = "Daily Index Value"
    ) +
    theme_minimal()
}

# Plot for raster_data_mosq
plot_mosq <- plot_raster_averages(prediction_raster_mosq, "Daily Harassment Over Time (Mosquito)")

# Plot for raster_data_oest
plot_oest <- plot_raster_averages(prediction_raster_oest, "Daily Harassment Over Time (Oestrid Fly)")

# Print the plots
print(plot_mosq)
print(plot_oest)


############################################################
# Visualizing Harassment across space and time
############################################################
# We can plot both mosquito and oestrid fly harassment across space and time, and compare, using a shiny application 

# Use the in-memory prediction rasters
raster_data_mosq <- prediction_raster_mosq
raster_data_oest <- prediction_raster_oest

# Define the new extent (adjusting ymin to 68.245)
new_extent <- ext(-140, -138.0492, 68.5, 69.05083)

# Crop the rasters to the new extent
raster_data_mosq <- crop(raster_data_mosq, new_extent)
raster_data_oest <- crop(raster_data_oest, new_extent)

# Extract time property
time_mosq <- terra::time(raster_data_mosq)
time_oest <- terra::time(raster_data_oest)

# Define UI
ui <- fluidPage(
  titlePanel("CHT Mosquito and Oestrid Fly Raster"),
  sidebarLayout(
    sidebarPanel(
      h4("Navigation"),
      actionButton("prev_layer", "Previous Layer"),
      actionButton("next_layer", "Next Layer"),
      br(),
      h4("Current Date and Time:"),
      verbatimTextOutput("current_time")
    ),
    mainPanel(
      fluidRow(
        column(6, leafletOutput(outputId = "map_mosq", height = "600px")),
        column(6, leafletOutput(outputId = "map_oest", height = "600px"))
      )
    )
  )
)

# Define Server
server <- function(input, output, session) {
  # Initialize the reactive value for the layer number
  layer <- reactiveVal(1)
  max_layer <- min(terra::nlyr(raster_data_mosq), terra::nlyr(raster_data_oest))
  
  # Observe the "Next" button and increase the layer
  observeEvent(input$next_layer, {
    if (layer() < max_layer) {
      layer(layer() + 1)
    }
  })
  
  # Observe the "Previous" button and decrease the layer
  observeEvent(input$prev_layer, {
    if (layer() > 1) {
      layer(layer() - 1)
    }
  })
  
  # Display the current date and time
  output$current_time <- renderText({
    paste("Mosquito Prediction:", time_mosq[layer()],
          "\nOestrid Fly Prediction:", time_oest[layer()])
  })
  
  # Reactive functions to extract the selected layer for each raster
  selected_layer_mosq <- reactive({
    raster_data_mosq[[layer()]]
  })
  
  selected_layer_oest <- reactive({
    raster_data_oest[[layer()]]
  })
  
  # Render Leaflet Map for Mosquito Predictions
  output$map_mosq <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addRasterImage(selected_layer_mosq(), colors = terrain.colors(100), opacity = 0.7) %>%
      addLegend(
        pal = colorNumeric(terrain.colors(100), values(raster_data_mosq), na.color = "transparent"),
        values = values(raster_data_mosq),
        title = "Mosquito Prediction",
        opacity = 0.7
      )
  })
  
  # Render Leaflet Map for Oestr Prediction
  output$map_oest <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addRasterImage(selected_layer_oest(), colors = terrain.colors(100), opacity = 0.7) %>%
      addLegend(
        pal = colorNumeric(terrain.colors(100), values(raster_data_oest), na.color = "transparent"),
        values = values(raster_data_oest),
        title = "Oestrid Fly Prediction",
        opacity = 0.7
      )
  })
  
  # Observe layer changes and update the Mosquito Map
  observe({
    leafletProxy("map_mosq") %>%
      clearImages() %>%
      addRasterImage(selected_layer_mosq(), colors = terrain.colors(100), opacity = 0.7)
  })
  
  # Observe layer changes and update the Oestr Map
  observe({
    leafletProxy("map_oest") %>%
      clearImages() %>%
      addRasterImage(selected_layer_oest(), colors = terrain.colors(100), opacity = 0.7)
  })
}

# Run the app
shinyApp(ui = ui, server = server)








