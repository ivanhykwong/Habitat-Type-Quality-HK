# Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong (1973-2022)

# Under Development

Supplementary materials used in the following studies:

Kwong, I. H. Y. (in progress). *Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong Using a 50-Year Archive of Remote Sensing Imagery* [Doctoral thesis, Department of Geography and Resource Management, The Chinese University of Hong Kong].

**Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2025). Spatial variations in forest succession rates revealed from multi-temporal habitat maps using Landsat imagery in subtropical Hong Kong. European Geosciences Union (EGU) General Assembly 2025, Vienna, Austria, 27 Apr–2 May 2025. https://doi.org/10.5194/egusphere-egu25-2667**

---

*Disclaimer: All datasets described here are for reference only. No express or implied warranty or representation is given to the accuracy or completeness of the data or its appropriateness for use in any particular circumstances.*

## GIS mapping results:

All raster layers (GeoTiff format) have a pixel size of 30 m covering the 1117-km2 terrestrial area in Hong Kong in this study. The time period of 1973–2022 was divided into 10 five-year periods in the mapping process. 

*   **HabitatMapHK_6class_yyyy-yyyy.tif**: Raster data showing the 6 habitat classes mapped in this study. Pixel values range from 1 to 6 representing woodland, shrubland, grassland, barren land, built-up area, and water respectively.

*   **HabitatMapHK_EstimatedArea.csv**: Area coverage (km2) of different habitat classes, as well as their confidence intervals, as mapped in this study. 
  
*   **HabitatMapHK_6class_ArcGISsymbology.lyrx**: Used to apply the suggested symbology in ArcGIS Pro, as shown in the figure below.

*   **ClassificationProbability_yyyy-yyyy.tif**: The probability values belonging to each class for every pixel. They were the intermediate products generated from the classification workflow and used to determine the final class with the highest probability and compute the forest index in this study. The sum of probabilities for all six classes equal to 1. A scale factor of 10000 was applied to the GeoTiff files for storage convenience. 

*   **HabitatMapHK_8class_yyyy-yyyy.tif**: Based on the 6-class outputs, two more classes are added in this product, including wetland (pixel value 7) and plantation (pixel value 8), to serve as inputs for the habitat quality model.

*   **HabitatQualityHK_yyyy-yyyy.tif**: Habitat quality maps produced in this study. The pixel value is a continuous variable ranging from 0 to 1, with 1 meaning the highest habitat quality.


## GIS supplementary data:

All datasets were collected and compiled from January to June 2024 and represent the conditions at that time. 

### Environmental Raster:

*   **DistanceFromCoast.tif**: Geometric distance (m) from the coastline.

*   **Elevation.tif**: Terrain height (m) from the LiDAR-based digital terrain model.

*   **Hillfire_10periods.tif**: Hill fires occurred in each five-year period, based on burn-area products by [Chan et al. (2023)](https://www.mdpi.com/2072-4292/15/6/1489) and manual digitisation for early years.

*   **Insolation.tif**: Annual amount of incoming solar radiation (kWh/m2) computed using SAGA GIS.

*   **Landslide_10periods.tif**: Landslides occurred in each five-year period, based on the [Enhanced Natural Terrain Landslide Inventory](https://portal.csdi.gov.hk/geoportal/?datasetId=cedd_rcd_1636520697377_96152&lang=en) ([Dias et al., 2009](https://hkss.cedd.gov.hk/hkss/filemanager/common/publications-resources/list-of-technical-papers/511_Dias%20et%20al%20(2009)_The%20Enhanced%20Natural%20Terrain%20Landslide%20Inventory.pdf)).

*   **Northness.tif**: Terrain aspect from 1 (due north) to -1 (due south) computed from the DTM.

*   **Precipitation.tif**: Annual precipitation (mm) from [Hong Kong Observatory (2022)](https://www.hko.gov.hk/en/cis/climahk.htm).

*   **Slope.tif**: Steepness (°) of the ground surface computed from the DTM.

*   **SoilCEC.tif**: Cation exchange capacity (CEC) (mmol/kg) of topsoil from [Luo et al. (2007)](https://books.google.com/books?id=ivzlMQAACAAJ).

*   **SoilOrganicMatter.tif**: Organic matter content (%) of topsoil from [Luo et al. (2007)](https://books.google.com/books?id=ivzlMQAACAAJ).

*   **Temperature.tif**: Annual mean temperature (°C) from [Morgan and Guénard (2019)](https://essd.copernicus.org/articles/11/1083/2019/).

*   **TopographicWetnessIndex.tif**: Amount of water accumulation due to topographic effects computed from the DTM.

*   **Typhoon_10periods.tif**: Wind speed (km/h) estimated from WindNinja based on maximum hourly mean wind records associated with typhoon events in each five-year period.

*   **WindSpeed.tif**: Mean wind speed (km/h) estimated from WindNinja based on monthly prevailing wind records.

### Human Activities:

*   **BuiltupAreas_10periods_shp.zip**: Shapefile (polygons) of built-up areas, with attributes on the years of construction (estimated from topographic maps) and density (high and low). It was used as threat factors in habitat quality mapping and variables in habitat changes.

*   **CountryParksProtectedAreas_shp.zip**: Shapefile (polygons) of protected areas (Country Parks, Special Areas, etc.), with attributes on the years of designation and revision. It was used as protection factors in habitat quality mapping and variables in habitat changes.

*   **PollutionSource_shp.zip**: Shapefile (polygons) of pollution sources (landfills, power stations, and incineration plants), with attributes on the years of construction and closure. It was used as threat factors in habitat quality mapping.

*   **Roads_10periods_shp.zip**: Shapefile (polylines) of roads, with attributes on the years of construction (estimated from topographic maps) and type (main and secondary). It was used as threat factors in habitat quality mapping.

### Mapping Reference:

*   **ForestIndex_FieldCollectedReferenceData.csv**: Field survey records of habitat types which were used to evaluate the forest index variable in this study.

*   **HabitatMapHK_FieldCollectedReferenceData.csv**: Field survey records of habitat types which were used to assess the habitat mapping results in this study.

*   **HabitatMapHK_OfficeInterpretedReferenceData.csv**: Reference points where the habitat class in each period was determined through visual interpretation of the aerial photographs and other historical records. The points were used for both training and validation of the habitat maps in this study.

*   **HabitatQualityHK_FieldSurveyedEcologicalValue2008.csv**: Field survey records of ecological values in 2008 which were used to evaluate the habitat quality maps in this study.

*   **LandsatHK_CrossSensorCalibrationPoints.csv**: Selected points that were assumed to remain unchanged over time and used to cross-calibrate different Landsat sensors in this study.

*   **LandsatHK_ImageMetadata.csv**: Metadata of the Landsat imagery (1,100 downloaded scenes and 607 valid scenes after pre-processing) acquired and processed in this study.

*   **Plantation_1975_1990_2008_2019.tif**: Pixels that were identified as plantations on four existing maps in different years (1975, 1990, 2008, 2019), as represented by the four layers contained in this raster file respectively. These pixels were used to help extract plantation class on the habitat map (when producing habitat quality) and denote areas with plantation activities (when modelling habitat changes) in this study.

*   **SpeciesObsHK_SpeciesChecklist.csv**: A species checklist of 7 taxa in Hong Kong (Plants, Butterflies, Birds, Reptiles, Dragonflies, Amphibians, Mammals) compiled from [AFCD](https://www.afcd.gov.hk/english/conservation/hkbiodiversity/speciesgroup/speciesgroup.html), [Hong Kong Biodiversity Information Hub](https://bih.gov.hk/en/hong-kong-species/index.html), and other secondary sources. Species of conservation concern are identified based on local assessments ([Corlett et al., 2000](https://www.researchgate.net/publication/364669045_Hong_Kong_vascular_plants_-_distribution_and_status); [Fellowes et al., 2002](https://www.researchgate.net/publication/268517493_Wild_animals_to_watch_terrestrial_and_freshwater_fauna_of_conservation_concern_in_Hong_Kong)), environmental protection laws, and national and global assessments. The checklist was used to match with the iNaturalist dataset to compute biodiversity metrics at grid levels and evaluate habitat quality maps in this study.

*   **SpeciesObsHK_SynonymList.csv**: A list of species name synonyms for matching names identified in iNaturalist and other secondary sources with the species checklist. It was used to pre-process the iNaturalist observation data and unify the species names of different records in this study.


## Analysis scripts:

### Part 1: Mapping Vegetation Habitats from a Satellite Image Time-Series

*   **P1_01_SearchAndDownloadFromGEE.ipynb**: Query and download all available Landsat 1-9 imagery covering the study area using Google Earth Engine. Atmospheric correction is performed if necessary.
*   **P1_02_Preprocess_part1.py**: Some basic pre-processing steps after downloading the images from cloud platform to local computer, such as mosaicking adjacent scenes and reprojecting to local coordinate system.
*   **P1_03_TopographicCorrection.R**: SCS+C topographic correction based on terrain slope, aspect, sun azimuth and sun elevation angles.
*   **P1_04_CrossSensorCal.R**: Cross-calibration of different Landsat sensors based on pseudo-invariant features, followed by computing variables for image classification.
*   **P1_05_ImageComposite.R**: Create image composites (median and standard deviation statistics) by combining all imagery acquried in the same period.
*   **P1_06_ExtractPixelValue.R**: Extract pixel values at the locations of reference points.
*   **P1_07_TrainingDataStat.R**: Summarise the characteristics of pixel values (e.g., spectral reflectance) of each habitat class and Landsat sensor.
*   **P1_08_TrainRFModel.R**: Train the random forest model, fuse probability outputs from each image, evaluate the model accuracies with cross-validation, and create the final model for classifying the entire dataset.
*   **P1_09_TestProcedures.R**: Modify the classification procedures and re-run the random forest models to evaluate their impacts on the classification accuracies.
*   **P1_10_ApplyModel.R**: Apply the random forest model and fusion steps to all images to create the habitat map for each period.
*   **P1_11_AreaCoverage.R**: Obtain the area coverage of each class on the habitat map as well as the confidence interval of the area estimates.
*   **P1_12_CompareFieldData.R**: Assess the accuracies of the habitat maps by overlaying with field-collected points and LiDAR height information at different times.
*   **P1_13_SurvivalAnalysis.R**: Analyse the number of years required for transitioning between vegetation classes as well as the correlations between transtion times and environmental variables.

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
