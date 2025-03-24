# Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong (1973-2022)

# Under Development

Supplementary materials used in the following studies:

Kwong, I. H. Y. (in progress). *Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong Using a 50-Year Archive of Remote Sensing Imagery* [Doctoral thesis, Department of Geography and Resource Management, The Chinese University of Hong Kong].

**Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2025). Spatial variations in forest succession rates revealed from multi-temporal habitat maps using Landsat imagery in subtropical Hong Kong. European Geosciences Union (EGU) General Assembly 2025, Vienna, Austria, 27 Apr–2 May 2025. https://doi.org/10.5194/egusphere-egu25-2667**


## GIS mapping results:

All raster layers (GeoTiff format) have a pixel size of 30 m covering the 1117-km2 terrestrial area in Hong Kong in this study. The time period of 1973–2022 was divided into 10 five-year periods in the mapping process.

*   **HabitatMapHK_6class_yyyy-yyyy.tif**: Raster data showing the 6 habitat classes mapped in this study. Pixel values range from 1 to 6 representing woodland, shrubland, grassland, barren land, built-up area, and water respectively.
  
*   **HabitatMapHK_6class_ArcGISsymbology.lyrx**: Used to apply the suggested symbology in ArcGIS Pro, as shown in the figure below.

*   **ClassificationProbability_yyyy-yyyy.tif**: The probability values belonging to each class for every pixel. They were the intermediate products generated from the classification workflow and used to determine the final class with the highest probability and compute the forest index in this study. The sum of probabilities for all six classes equal to 1. A scale factor of 10000 was applied to the GeoTiff files for storage convenience. 

*   **HabitatMapHK_8class_yyyy-yyyy.tif**: Based on the 6-class outputs, two more classes are added in this product, including wetland (pixel value 7) and plantation (pixel value 8), to serve as inputs for the habitat quality model.

*   **HabitatQualityHK_yyyy-yyyy.tif**: Habitat quality maps produced in this study. The pixel value is a continuous variable ranging from 0 to 1, with 1 meaning the highest habitat quality.


## GIS supplementary data:

Unless otherwise specified, all datasets were collected and compiled from January to June 2024 and represent the conditions at that time.

### Environmental Raster:

*   **DistanceFromCoast.tif**: Geometric distance (m) from the coastline.

*   **Elevation.tif**: Terrain height (m) from the LiDAR-based digital terrain model.

*   **Hillfire_10periods.tif**: Hill fires occurred in each five-year period, based on burn-area products by Chan et al. (2023) and manual digitisation for early years.

*   **Insolation.tif**: Annual amount of incoming solar radiation (kWh/m2) computed using SAGA GIS.

*   **Landslide_10periods.tif**: Landslides occurred in each five-year period, based on the Enhanced Natural Terrain Landslide Inventory (Dias et al., 2009).

*   **Northness.tif**: Terrain aspect from 1 (due north) to -1 (due south) computed from the DTM.

*   **Precipitation.tif**: Annual precipitation (mm) from Hong Kong Observatory (2022).

*   **Slope.tif**: Steepness (°) of the ground surface computed from the DTM.

*   **SoilCEC.tif**: Cation exchange capacity (CEC) (mmol/kg) of topsoil from Luo et al. (2007)

*   **SoilOrganicMatter.tif**: Organic matter content (%) of topsoil from Luo et al. (2007)

*   **Temperature.tif**: Annual mean temperature (°C) from Morgan and Guénard (2019).

*   **TopographicWetnessIndex.tif**: Amount of water accumulation due to topographic effects computed from the DTM.

*   **Typhoon_10periods.tif**: Wind speed (km/h) estimated from WindNinja based on maximum hourly mean wind records associated with typhoon events in each five-year period.

*   **WindSpeed.tif**: Mean wind speed (km/h) estimated from WindNinja based on monthly prevailing wind records.

### Human Activities:

*   **BuiltupAreas_10periods_shp.zip**: Shapefile (polygons) of built-up areas, with attributes on the years of construction (estimated from topographic maps) and density (high and low). It was used as threat factors in habitat quality mapping and variables in habitat changes.

*   **CountryParksProtectedAreas_shp.zip**: Shapefile (polygons) of protected areas (Country Parks, Special Areas, etc.), with attributes on the years of designation and revision. It was used as protection factors in habitat quality mapping and variables in habitat changes.

*   **PollutionSource_shp.zip**: Shapefile (polygons) of pollution sources (landfills, power stations, and incineration plants), with attributes on the years of construction and closure. It was used as threat factors in habitat quality mapping.

*   **Roads_10periods_shp.zip**: Shapefile (polylines) of roads, with attributes on the years of construction (estimated from topographic maps) and type (main and secondary). It was used as threat factors in habitat quality mapping.

### Mapping Reference:

*   **ForestIndex_FieldCollectedReferenceData.csv**:

*   **HabitatMapHK_EstimatedArea.csv**:

*   **HabitatMapHK_FieldCollectedReferenceData.csv**:

*   **HabitatMapHK_OfficeInterpretedReferenceData.csv**:

*   **HabitatQualityHK_FieldSurveyedEcologicalValue2008.csv**:

*   **LandsatHK_CrossSensorCalibrationPoints.csv**:

*   **LandsatHK_ImageMetadata.csv**:

*   **Plantation_1975_1990_2008_2019.tif**:

*   **SpeciesObsHK_SpeciesChecklist.csv**:

*   **SpeciesObsHK_SynonymList.csv**:


## Analysis codes:

### Part 1: Mapping Vegetation Habitats from a Satellite Image Time-Series

*   **P1_01_SearchAndDownloadFromGEE.ipynb**:
*   **P1_02_Preprocess_part1.py**:
*   **P1_03_TopographicCorrection.R**:
*   **P1_04_CrossSensorCal.R**:
*   **P1_05_ImageComposite.R**:
*   **P1_06_ExtractPixelValue.R**:
*   **P1_07_TrainingDataStat.R**:
*   **P1_08_TrainRFModel.R**:
*   **P1_09_TestProcedures.R**:
*   **P1_10_ApplyModel.R**:
*   **P1_11_AreaCoverage.R**:
*   **P1_12_CompareFieldData.R**:
*   **P1_13_SurvivalAnalysis.R**:

### Part 2: Computing Habitat Quality Maps with Reference to Habitat Type Information

*   **P2_01_ModifyInputMap.R**:
*   **P2_02_SensitivityParameter.R**:
*   **P2_03_RunInVESTModel.R**:
*   **P2_04_CompareFieldSurvey.R**:
*   **P2_05_CompareSpeciesObservation.R**:
*   **P2_06_DegradationRarity.R**:

### Part 3: Modelling Habitat Changes in Relation to Natural and Human Factors

*   **P3_01_CreateEnvRaster.R**:
*   **P3_02_ForestIndex.R**:
*   **P3_03_ExtractPixelValue.R**:
*   **P3_04_EFA.R**:
*   **P3_05_SEM.R**:
*   **P3_06_ModelValidation.R**:
*   **P3_07_YearChange.R**:

---

*Last updated in March 2025*
