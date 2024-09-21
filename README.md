# Arctic_Insect_Harassment_Models
Repository containing mosquito and oestrid fly harassment GAM models described in Hein et al. 2024.

Attached in this repository are the the .rds files containing the mosquito and oestrid fly models descridbed in "citation". These models can be used to predict insect harassment by mosquitoes and oestrid fly harassment in arctic ecosystems remotely using a trusted DSM and open source ERA5 land weather model variables. We describe the process to do so below:

Example 1: Predicting harassment onto a dataframe of GPS locations of caribou.

**Sample Covariates** 
1. Download covariates
   **Mosquito Covariates**
     - 30m DSM, we chose ALOS World 3D-30m (https://www.eorc.jaxa.jp/ALOS/en/dataset/aw3d30/aw3d30_e.htm)
         - To reduce download size, limit the raster dimensions to your study area
     - ERA5 hourly Land Weather model (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form)
         - Download 2m temperature, 10m u-component of wind and 10m v-component of wind to spatially and temporally cover your study's area and time period for each year.
         - Download volumetric soil water layer 1 to spatially cover your study area and temporally cover the summers (we used 05/01 to 08/31) of interest for each year/
         - Download snow depth to spatially cover your study area and temporally cover the time period in which snow melt occurs for each year (we used 05/01 to 07/15 but will vary     
           depending on location specific climate)
         - Download precipitation to spatially cover your study area and temporally cover at minimum two weeks before your study's time frame and up to the end of your study period for   
           each year.
         - NOTE: To avoid large raster files, download only the time periods of interest that you expect mosquitoes to be present. For example, a study from 2015-2018 you dont need to                download covariates for february of each year.
     **Oestrid fly Covariates**
     - ERA5 Land Weather model (https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-land?tab=form)
        - Download 2m temperature to spatially cover the study area and temporally cover the time period from prior to spring temperatures reaching above 0 degrees celsius to the end of             your study period for each year (we used 05/01 to 0/31).
2. Upload Covariates as single layer (DSM) or multidimensional rasters (ERA5 weather covariates)
3. Merge Multidimensional rasters
     -For each ERA5 hourly covariates they will grouped by month when you download them and you will need to merge them using the Merge Multidimensional Raster tool in the Multidimension        Tools toolbox. ArcGIS will recognize the timestamps for the raster layers and will automatically merge the monthly rasters into one raster, don't worry if you have gaps in the   
      timeframe, the merge tool retains the gaps as blanks and should not cause problems later.
     - You should now have 6 hourly multidimensional rasters spanning the time frame of your study, one for each of the ERA5 covariates outlined above (2m temperature, 10m u-component of         wind, 10 v-component of wind, volumetric soil water layer 1, snow depth, and precipitation). Check the time component in the properties to ensure you have the correct temporal             coverage.
4. Manipulate Rasters
     - Soil Moisture: using the aggregate multidimensional raster tool in the multidimension tools toolbox, calculate the mean soil moisture for each summer (we used 05/01 trhough 08/31), 
       this will be your soil moisture covariate for prediciting mosquito harassment.
     - Snow Depth: Using the aggregate multidimensional
     - Growing Degree Day: First using the minus tool from the spatial analyst toolbox, subtract 273.15 from the temperature raster to convert the units from kelvin to celsius. Then        using the aggregate multidimensional raster tool in the multidimension tools toolbox, create two new rasters from the temperature raster. The first will be aggregating over 24             hours using an aggregation definition of minimum to find the daily minimum temperature. The the second will be aggregating over 24 hours using an aggregation definition of maximum         to find the daily maximum temperature. Then you will calculate growing degree day with a base temperature of 0, this is done by subtracting the minimum daily temperature from the          maximum daily temperature and dividing by 2. This can be done using the python code provided below:

               from arcpy.sa import Raster, RasterCalculator
       
                # Define the paths to your rasters
                raster1_path = "DailyMaxTemperature.crf"
                raster2_path = "DailyMinTemperature.crf"
                
                # Create Raster objects from the file paths
                raster1 = Raster(raster1_path)
                raster2 = Raster(raster2_path)
                
                # Perform raster calculation
                out_rc_multi_raster = RasterCalculator([raster1, raster2], ["x", "y"], "(x + y) / 2")
                
                # Define the output path in your project directory
                output_path = "C:/Users/yourdirectory/gdd.crf"
                
                # Save the output raster
                out_rc_multi_raster.save(output_path)

       This will create a daily growing degree day raster, however we do only care about days above 0 so we need transform all negative gdd days to 0. You can do this using the python   
       code provided below:

             import arcpy
            from arcpy.sa import Con, Raster
            
            # Replace "GDD" with the actual name of your raster
            gdd = Raster("gdd.crf")
            
            # Use the Raster Calculator to set negative values to 0
            result_raster = Con(gdd < 0, 0, gdd)
            
            # Save the result back to the original raster
            result_raster.save("GDD")
       
       Now we have a growing degree day raster we can work with.
    -Topographic Position Index (TPI): To calculate TPI from the DSM you will first need to broaden the resolution of the DSM. You can do so using the aggregate tool from the spatial          analyst tools toolbox in arcGIS, setting the cell factor by 33.33 and setting the aggregation technique to mean (this will create a raster that changes the resolution from 30m to 1km).
    Next, using the focal statistics tool from the spatial analyst tools toolbox, input the 1 km DSM aand set the neighborhood to a circle with a 5 cell radius and calulcate the mean      
    elevation to calculate a 5 km neighborhood around each cell of the 1km DSM. Then use the raster calculator to subtract the 1km DSM from the 5km focal DSM. Now  we have a TPI raster 
    with a resolution of 1km and a neighborhood of 5km.
5. Sample covariates onto dataframe
    - To start ensure your dataframe has a unique ID, longitude, latitude, date time, date time minus 2 weeks, and date time start (a set month and day of each respective data points year
      you selected prior to spring temperatures reaching above 0 degrees celsius in step 1, in our study area we used 05/01/2024, 05/01/2023... etc.). We recommend setting all
      date time columns to UTC.
    - Using the add data function in arcGIS pro, add the dataframe as X Y point data. Then in the properties -> time section of this layer select "filter layer based on attribute values         and select "each feature has a single time field" then the date time column as the time field.
    - Using the sample tool from the spatial analyst tools toolbox, you can now sample the spatial and temporal covariates. For temperature, and both U and V components of wind, use the 
      sample tool from the spatial analyst tools toolbox and set the input raster as one of the aforementioned covariates, then your X Y dataframe as the input location raster or       
      features. Select to process as multidimensional and the unique ID field as your Object ID column. Set the dimension as StdTime, and start field as the date time column in your
      dataset, leave the statistics type blank. Hit run and the tool will create a table with the object ID column and sampled raster value at the time and location of each dataset. This
      can then be exported as a .csv and appended to your original dataset in a column named "temperature", "wind_v", and "wind_u" respectively. Calculate the square root of the summed          squared wind_u and wind_V column (=SQRT(A2^2+2B2^2) in excel)) to create a wind speed column labeled "wind".
    - To sample the two week precipitation covariate simply run the same sampling parameters described above after inputting the precipitation raster setting the end field or value as you
      date time minus 2 weeks column and change the statistic type to sum. Append the resulting dataframe to your original dataset as a column called "2wkprecip". Note that depending
      on the power of your computer and size of your X Y dataframe this may take a from a few minutes to a few hours.
    - To sample cumulative GDD again run the sampling parameters described above after inputting the GDD raster setting the end field or value as you date time start
      column and change the statistic type to sum. Append the resulting dataframe to your original dataset as a column called "cumGDD".
    - To sample soil moisture, again run the same sampling parameters described for the temperature and wind covariates after inputting the aggragated soil moisture raster and setting.
      Append the resulting dataframe to your original dataset as a column called "soilmoisture".
    - To sample soil moisture, again run the same sampling parameters described for the temperature and wind covariates after inputting the snow depth raster and setting.
      Append the resulting dataframe to your original dataset as a column called "soilmoisture".
      










       
   
