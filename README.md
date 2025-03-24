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

*   **DistanceFromCoast.tif**: 

*   **Elevation.tif**: 

*   **Hillfire_10periods.tif**: 

*   **Insolation.tif**: 

*   **Landslide_10periods.tif**: 

*   **Northness.tif**: 

*   **Precipitation.tif**: 

*   **Slope.tif**: 

*   **SoilCEC.tif**: 

*   **SoilOrganicMatter.tif**: 

*   **Temperature.tif**: 

*   **TopographicWetnessIndex.tif**: 

*   **Typhoon_10periods.tif**: 

*   **WindSpeed.tif**: 

### Human Activities:

*   Buildings_10periods_shp.zip: 

*   CountryParksProtectedAreas_shp.zip:

*   PollutionSource_shp.zip:

*   Roads_10periods_shp.zip:

### Mapping Reference:

*   ForestIndex_FieldCollectedReferenceData.csv:

*   HabitatMapHK_EstimatedArea.csv:

*   HabitatMapHK_FieldCollectedReferenceData.csv:

*   HabitatMapHK_OfficeInterpretedReferenceData.csv:

*   HabitatQualityHK_FieldSurveyedEcologicalValue2008.csv:

*   LandsatHK_CrossSensorCalibrationPoints.csv:

*   LandsatHK_ImageMetadata.csv:

*   Plantation_1975_1990_2008_2019.tif:

*   SpeciesObsHK_SpeciesChecklist.csv:

*   SpeciesObsHK_SynonymList.csv:


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
