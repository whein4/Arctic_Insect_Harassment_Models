# Arctic_Insect_Harassment_Models
Repository containing mosquito and oestrid fly harassment models and instructions for their use in calculating the CHT insect indices described in Hein et al., 2025.

Attached in this repository are the the .rds files containing the mosquito and oestrid fly models descridbed in Hein et al., 2025. These models can be used to predict insect harassment by mosquitoes and oestrid fly harassment in arctic ecosystems remotely using a trusted DSM and open source ERA5 land weather model variables. We provide code to predict a single year's hourly or daily insect harassment using the CHT indices in PredictingCHT_basecode.R as well as code to predict multiple years for hindcasting analysis (predicting daily) in PredictingCHT_loop.R.

The process is outlined below:

 Step 1: Downloading ERA5 rasters 
   - Download the required predictors for both the CHT mosquito and oestrid fly models from ERA5.
       - Mosquito Model predictors and ERA5 source variable:
            temperature = 2m_temperature
            wind speed = 10m_u_component_of_wind and 10m_v_component_of_wind
            2 week sum of precipitation = total_precipitation
            days since snow melt = snow_depth
            soil moisture = volumetric_soil_water_layer_1
        - Oestrid Model:
            temperature = 2m_temperature
            cumulative growing degree day = 2m_temperature
   - Pull directly from the ERA5 online data store (open access; cds.climate.copernicus.eu) into R
 Step 2: Manipulate ERA5 rasters
  - Manipulate the raw ERA5 rasters into CHT predictors
    - Temperature converted from kelvin to celsius
    - Derive cumulative growing degree day from temperature
    - Derive wind speed from U and V components of wind speed
    - Derive average soil moisture across the summer from soil moisture
    - Derive days since snow melt from snow depth
    - Derive 2-week precipitation from precipication
 Step 3: Download DSM and Calculate topographic position index (TPI)
  - You will need to compute TPI from ALOS DSM (code for this is provided in PredictingCHT_basecode.R)
    - You will need to download DSM tiles manually via https://www.eorc.jaxa.jp/ALOS/en/dataset/aw3d30/aw3d30_e.htm
    - Then you will merge the DSM tiles into a single DSM layer
    - Then you will resample to align resolutions with other predictors and calculate TPI using a 5 km neighborhood
    - Once TPI_5km_raster.tif has been created, you will not need recalculate for every year, we provide code to pull this .tif directly for future use
 Step 4: Standardize rasters
  - You will then standardize the resolution, extent, time, and layer number of raster predictors to facilitate standardized predictions
 Step 5: Create a time of day raster
  - To predict oestrid fly harassment you will need a time of day raster, this can simply be created by mirroring one of the standardized rasters and filling with time values
 Step 6: Predict Using CHT model
  - You will then stack the respective mosquito and oestrid fly predictors and predict, creating dataframes which are then manipulated back into a raster format
 Step 7: Visualize the results
  - In PredictingCHT_basecode.R we provide example code to visualize the resulting predicted CHT mosquito and oestrid fly multi-dimensional rasters.


For additional information on working with multidimensional rasters please refer to:
Introduction to Working with Raster Data in R: https://www.neonscience.org/resources/learning-hub/tutorials/introduction-working-raster-data-r
R for Spatial Statistics: https://gsp.humboldt.edu/olm/R/03_03_RasterFiles_Terra.html
Introduction to Geospatial Raster and Vector Data with R: https://ucsb-dreamlab.github.io/r-raster-vector-geospatial/05-raster-multi-band-in-r.html
Advanced tools for raster data in R: https://mhallwor.github.io/_pages/advanced_Rasters
Spatial analysis with R: https://rfrelat.github.io/Spatial2_MultiExamples.html

For questions about CHT insect indices and their prediction using the code provided in this repository contact me at william.hein@mail.mcgill.ca
















       
   
