# Arctic_Insect_Harassment_Models
Repository containing mosquito and oestrid fly harassment GAM models described in Hein et al. 2024.

Attached in this repository are the the .rds files containing the mosquito and oestrid fly models descridbed in "citation". These models can be used to predict insect harassment by mosquitoes and oestrid fly harassment in arctic ecosystems remotely using a trusted DSM and open source ERA5 land weather model variables. We describe the process to do so below:

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
     -For each ERA5 hourly covariates they will grouped by month when you download them and you will need to merge them using the Merge Multidimensional Raster tool in the Multidimension           Tools toolbox. ArcGIS will recognize the timestamps for the raster layers and will automatically merge the monthly rasters into one raster, don't worry if you have gaps in the   
      timeframe, the merge tool retains the gaps as blanks and should not cause problems later.
     - You should now have 6 hourly multidimensional rasters spanning the time frame of your study, one for each of the ERA5 covariates outlined above (2m temperature, 10m u-component of         wind, 10 v-component of wind, volumetric soil water layer 1, snow depth, and precipitation). Check the time component in the properties to ensure you have the correct temporal             coverage.
4. 










       
   
