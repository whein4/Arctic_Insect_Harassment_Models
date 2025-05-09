##################################################################################
##################################################################################
##################################################################################


# This R code is designed to hindcast insect harassment
# using the CHT insect indices introduced in Hein et al. 2025.
# Please refer to the documentation provided on Github
# https://github.com/whein4/Arctic_Insect_Harassment_Models/blob/main/README.md

# Contact: william.hein@mail.mcgill.ca


##################################################################################
##################################################################################
##################################################################################


# This code is a looped version of the PredictingCHT_basecode.R to allow for efficient hindcasting processing
# Note that this code is currently written to predict daily insect harassment, rather then hourly. For hourly CHT loop, refer to PredictingCHT_basecode.R
# We recommend that you use PredictingCHT_basecode.R to ensure your processing single years correctly prior to using this looped version
# Please reference PredictingCHT_basecode.R for code description

#install necessary libraries
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




##################################################################################
# Loop
##################################################################################

# Run on R version R version 4.4.2

# Select years to process (years <- as.character(seq(1950, 2024, by = 1)))
# Set your ERA5 key (https://cds.climate.copernicus.eu/profile)
# Select working directory (setwd())
# Select predictor's output directory (output_dir())
# Select final raster's output directory (FINAL_output_dir()), near the bottom of the script
# Select area of interest (area <- c(70.15, -146, 67.5, -135.3))
# Set the location for the TPI Raster (tpi_raster <- rast("summerrangeTPI.tif"))'
# Set location for the mosquito and oestrid fly models ( gam_oest <- readRDS("C:/Users/willi/OneDrive/Desktop/MASTERS/Project_data/Final Dataframes/Models/GAMoestFINAL.rds") & gam_mosq <- readRDS("C:/Users/willi/OneDrive/Desktop/MASTERS/Project_data/Final Dataframes/Models/GAMmosqFINAL.rds"))



# Define the years to process
years <- as.character(seq(1986, 1986, by = 1))


# Loop through each year
for (year in years) {
  year_str <- as.character(year)
  
  # Set the API key
  wf_set_key(key = "f967427f-5bf8-4a72-b787-67d4c82fee55")
  
  # Define output directories
  setwd("C:/your/working/directory")
  output_dir <- "C:/your/output/directory"
  
  # Define fixed parameters
  months <- c("05", "06", "07", "08")
  snowmonth <- c("05", "06", "07")
  precipmonth <- c("04", "05", "06", "07", "08")
  days <- sprintf("%02d", 1:31)
  times <- sprintf("%02d:00", 0:23)
  area <- c(70.15, -146, 67.5, -135.3)  
  
  
  # Define the start and end dates
  start_date <-  as.POSIXct(paste0(year, "-05-01 00:00"), tz = "UTC")
  index_start_date <- as.POSIXct(paste0(year, "-06-01 12:00"), tz = "UTC")
  index_end_date <- as.POSIXct(paste0(year, "-08-31 12:00"), tz = "UTC")
  end_date <- as.POSIXct(paste0(year, "-08-31 23:00"), tz = "UTC")
  snow_melt_end <- as.POSIXct(paste0(year, "-07-31 23:00"), tz = "UTC")
  snow_start_date <- as.POSIXct(paste0(year, "-05-01 00:00"), tz = "UTC")
  precipitation_start <- as.POSIXct(paste0(year, "-04-01 00:00"), tz = "UTC")
  
  
  
  ##########################################################
  #download ERA5 data
  ##########################################################
  ################
  #temp
  ###############
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
  
  ################
  #wind-u
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
  #wind-v
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
  #soil moisture
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
  #snow depth
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
  #precipitation
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
  
  
  
  
  
  
  # Check that the rasters loaded correctly
  print(temp_raw)
  print(snow_raw)
  print(wind_v_raw)
  print(wind_u_raw)
  print(precip_raw)
  print(soil_raw)
  
  
  #######################################################################################
  #manipulate ERA5 rasters
  ######################################################################################
  
  
  ######################################
  #temp to celsius
  ######################################
  
  # Convert temperature from Kelvin to Celsius
  temp_celsius <- temp_raw - 273.15
  
  ######################################
  #temp to cGDD
  ######################################
  
  # Aggregate the hourly raster to daily minimum
  temp_min_daily <- tapp(temp_celsius, "days", fun = min)
  
  # Aggregate the hourly raster to daily maximum
  temp_max_daily <- tapp(temp_celsius, "days", fun = max)
  
  # Compute the daily average temperature
  gdd_raster <- (temp_min_daily + temp_max_daily) / 2
  
  # Set all negative values to zero
  gdd_raster[gdd_raster < 0] <- 0
  
  # Inspect the resulting GDD raster
  print(gdd_raster)
  
  # Calculate cumulative GDD
  gdd_cumulative <- app(gdd_raster, function(x) {
    cumsum(x)  # Compute cumulative sum across layers for each pixel
  })
  
  # Inspect the resulting cumulative GDD raster
  print(gdd_cumulative)
  
  ######################################
  #wind u/v to wind speed
  ######################################
  
  # Calculate the square root of the sum of squares
  wind_speed <- sqrt(wind_u_raw^2 + wind_v_raw^2)
  
  # Inspect the resulting raster
  print(wind_speed)
  
  ######################################
  #soil to summer avg soil
  ######################################
  # Calculate the average across all layers
  soil_avg <- mean(soil_raw, na.rm = TRUE)
  
  # Inspect the resulting raster
  print(soil_avg)
  
  ######################################
  #snow depth to days since snow melt
  ######################################
  # Aggregate the snow_raw raster to daily minima
  snow_daily_min <- tapp(snow_raw, "days", fun = min)
  
  # Find the layer index where snow depth < 0.2 for the first time
  first_melt_layer <- app(snow_daily_min, function(x) {
    melt_index <- which(x < 0.2)[1]  # Find the first occurrence of snow depth < 0.2
    if (!is.na(melt_index)) {
      melt_index  # Return the layer index
    } else {
      NA  # No snow melt for this pixel
    }
  })
  
  # Inspect the result
  print(first_melt_layer)
  
  # Define the start and end dates
  start_date <- as.Date(paste0(year, "-05-31"))
  end_date <- as.Date(paste0(year, "-08-31"))
  all_dates <- seq(start_date, end_date, by = "day")
  num_days <- length(all_dates)
  
  # Create an empty raster stack with the same dimensions and resolution as snow_daily_min
  index_stack <- rast(snow_daily_min, nlyr = num_days)
  
  # Assign sequential indices (1, 2, ..., num_days) to the raster stack
  values(index_stack) <- rep(1:num_days, each = ncell(index_stack))
  
  # Assign date-based names to the layers
  index_stack <- setNames(index_stack, paste0("d_", all_dates))
  
  
  # Expand first_melt_layer to match the number of layers in the index_stack
  expanded_first_melt <- rast(lapply(1:num_days, function(i) first_melt_layer))
  
  # Ensure the expanded raster has the same layer names as the index stack
  expanded_first_melt <- setNames(expanded_first_melt, paste0("d_", all_dates))
  
  # Subtract the expanded first melt raster from the index stack
  days_since_snow_melt <- index_stack - expanded_first_melt
  
  # Inspect the result
  print(days_since_snow_melt)
  
  ######################################
  #precip  to 2 week sum of precip
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
  
  #######################################################################################
  #load TPI
  ######################################################################################
  tpi_raster <- rast("TPI_5km_raster.tif")
  #Plot 5km TPI raster
  plot(tpi_raster)
  
  
  #########################################################################################################
  #standardizing rasters
  #########################################################################################################
  
  ######################################
  #subset 2 week sum of precip
  ######################################
  precip_2wk <- subset(precip_2wk, 31:153)
  
  ######################################
  #set crs
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
  #set # layers
  ######################################
  #we only need 12:00 from June 1 - August 31 from hourly rasters
  # Define the sequence of layers to subset
  layer_indices <- seq(764, 2948, by = 24)
  
  # Subset each raster individually
  temp_celsius_proj_daily <- subset(temp_celsius_proj, layer_indices)
  wind_speed_proj_daily <- subset(wind_speed_proj, layer_indices)
  
  #subset daily rasters to ignore May
  precip_proj_daily <- subset(precip_proj, 32:123)
  gdd_proj_daily <- subset(gdd_proj, 32:123)
  dssm_proj_daily <- subset(dssm_proj, 32:123)
  
  #for single layer rasters:
  expand_single_layer <- function(single_layer, n_layers) {
    # Duplicate the single layer n_layers times
    expanded_raster <- terra::rast(replicate(n_layers, single_layer))
    
    return(expanded_raster)
  }
  
  # Expand soil_avg to 2952 layers
  soil_proj_daily <- expand_single_layer(soil_proj, 92)
  
  # Expand tpi to 2952 layers
  tpi_proj_daily <- expand_single_layer(tpi_proj, 92)
  
  
  print(temp_celsius_proj_daily)
  print(wind_speed_proj_daily)
  print(precip_proj_daily)
  print(gdd_proj_daily)
  print(dssm_proj_daily)
  print(soil_proj_daily)
  print(tpi_proj_daily)
  
  ######################################
  # Set Extent 
  ######################################
  # List of rasters
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
  
  
  ######################################
  # Set Resolution After Cropping
  ######################################
  # Use one of the rasters as the reference
  reference_raster <- tpi_proj_crop
  
  # Resample all rasters to align with the reference raster
  temp_celsius_proj_res_crop <- terra::resample(temp_celsius_proj_crop, reference_raster)
  wind_speed_proj_res_crop <- terra::resample(wind_speed_proj_crop, reference_raster)
  soil_proj_res_crop <- terra::resample(soil_proj_crop, reference_raster)
  precip_proj_res_crop <- terra::resample(precip_proj_crop, reference_raster)
  dssm_proj_res_crop <- terra::resample(dssm_proj_crop, reference_raster)
  gdd_proj_res_crop <- terra::resample(gdd_proj_crop, reference_raster)
  tpi_proj_res_crop <- tpi_proj_crop  # No resampling needed
  
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
  
  
  ############################################################################
  #create time raster
  ############################################################################
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
  
  ######################################
  #set time
  ######################################
  
  # Assign hourly time to the expanded raster
  time_reference <- seq(index_start_date,
                        index_end_date,
                        by = "day")
  terra::time(temp_celsius_proj_res_crop) <- time_reference
  terra::time(wind_speed_proj_res_crop) <- time_reference
  terra::time(precip_proj_res_crop) <- time_reference
  terra::time(gdd_proj_res_crop) <- time_reference
  terra::time(dssm_proj_res_crop) <- time_reference
  terra::time(tpi_proj_res_crop) <- time_reference
  terra::time(soil_proj_res_crop) <- time_reference
  terra::time(time_proj_res_crop) <- time_reference
  
  
  #check mosquito predictors
  print(temp_celsius_proj_res_crop)
  print(wind_speed_proj_res_crop)
  print(soil_proj_res_crop)
  print(precip_proj_res_crop)
  print(dssm_proj_res_crop)
  print(tpi_proj_res_crop)
  
  #check oestrid fly predictors
  print(temp_celsius_proj_res_crop)
  print(gdd_proj_res_crop)
  print(time_proj_res_crop)
  
  
  #Clean the environment
  objects_to_keep <- c(
    "temp_celsius_proj_res_crop", 
    "wind_speed_proj_res_crop", 
    "precip_proj_res_crop", 
    "gdd_proj_res_crop", 
    "dssm_proj_res_crop", 
    "soil_proj_res_crop", 
    "tpi_proj_res_crop",
    "time_proj_res_crop",
    "years",
    "pb",
    "FINAL_output_dir",
    "output_dir",
    "year_str"
  )
  
  # Remove all objects except the ones in the list
  rm(list = setdiff(ls(), objects_to_keep))
  
  # Garbage collection to free memory
  gc()
  
  ##################################################################################################################
  #Predict Using CHT Rasters
  ##################################################################################################################
  
  # Load models
  gam_oest <- readRDS("C:/Users/willi/OneDrive/Desktop/MASTERS/Project_data/Final Dataframes/Models/GAMoestFINAL.rds")
  gam_mosq <- readRDS("C:/Users/willi/OneDrive/Desktop/MASTERS/Project_data/Final Dataframes/Models/GAMmosqFINAL.rds")
  
  ############################################################################
  #predict mosquito harassment
  ############################################################################
  
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
  prediction_raster_mosq <- terra::rast(predictors_mosq$ERAtemp) # Use one of the predictors as a template
  terra::values(prediction_raster_mosq) <- NA # Clear the values
  
  # Assign predicted values back to the raster (match indices)
  all_indices <- seq_len(ncell(prediction_raster_mosq))
  valid_indices <- which(!is.na(terra::values(predictors_mosq$ERAtemp)))
  terra::values(prediction_raster_mosq)[valid_indices] <- predicted_values_mosq_index
  
  # View or save the result
  print(prediction_raster_mosq)
  min(prediction_raster_mosq)
  max(prediction_raster_mosq)
  
  ############################################################################
  #predict oestrid harassment
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
  predicted_values_oest_sq <- sqrt(predicted_values_oest)
  
  # Create a new raster to store predictions
  prediction_raster_oest <- terra::rast(predictors_oest$ERAtemp) # Use one of the predictors as a template
  terra::values(prediction_raster_oest) <- NA # Clear the values
  
  # Assign predicted values back to the raster (match indices)
  all_indices <- seq_len(ncell(prediction_raster_oest))
  valid_indices <- which(!is.na(terra::values(predictors_oest$ERAtemp)))
  terra::values(prediction_raster_oest)[valid_indices] <- predicted_values_oest_sq
  
  # View or save the result
  print(prediction_raster_oest)
  
  ################################################################################
  #Save the predicted rasters in a new fresh file
  ################################################################################
  # Define the output directory
  FINAL_output_dir <- "C:/your/FINALoutput/directory"
  
  # Save rasters
  output_file_mosq <- paste0(FINAL_output_dir, "/mosq_pred_", year_str, ".nc")
  writeCDF(prediction_raster_mosq, output_file_mosq, overwrite = TRUE)
  
  output_file_oest <- paste0(FINAL_output_dir, "/oest_pred_", year_str, ".nc")
  writeCDF(prediction_raster_oest, output_file_oest, overwrite = TRUE)
  
  # Ensure key variables are retained
  objects_to_keep <- c(
    "years", 
    "pb", 
    "FINAL_output_dir", 
    "output_dir"
  )
  
  # Remove only temporary objects
  rm(list = setdiff(ls(), objects_to_keep))
  
  # Garbage collection to free memory
  gc()
}










