# Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong (1973-2022)

## Supplementary materials used in the following studies:

Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (Manuscript submitted for publication). **Multi-temporal assessment of habitat quality in Hong Kong using the InVEST model and multi-taxa citizen science data**.

Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2026). **Integrating five decades of Landsat imagery for territory-wide habitat mapping and change detection in a subtropical metropolitan city**. *Remote Sensing Applications: Society and Environment, 41*, 101910. https://doi.org/10.1016/j.rsase.2026.101910

Kwong, I. H. Y. (2025). ***Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong Using a 50-Year Archive of Remote Sensing Imagery*** [Doctoral thesis, Department of Geography and Resource Management, The Chinese University of Hong Kong]. https://repository.lib.cuhk.edu.hk/en/item/cuhk-3584537

Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2025). **Spatial variations in forest succession rates revealed from multi-temporal habitat maps using Landsat imagery in subtropical Hong Kong**. European Geosciences Union (EGU) General Assembly 2025, Vienna, Austria, 27 Apr–2 May 2025. https://doi.org/10.5194/egusphere-egu25-2667. [Poster Presentation: https://presentations.copernicus.org/EGU25/EGU25-2667_presentation-h291057.pdf]

---

*Disclaimer: All datasets described here are for reference only. No express or implied warranty or representation is given to the accuracy or completeness of the data or its appropriateness for use in any particular circumstances.*

## GIS mapping results:

All raster layers (GeoTiff format) have a pixel size of 30 m covering the 1117-km2 terrestrial area in Hong Kong in this study (Hong Kong 1980 Grid coordinate system). The time period of 1973–2022 was divided into 10 five-year periods in the mapping process. 

*   **HabitatMapHK_6class_yyyy-yyyy.tif**: Raster data showing the 6 habitat classes mapped in this study. Pixel values range from 1 to 6 representing woodland, shrubland, grassland, barren land, built-up area, and water respectively.

*   **HabitatMapHK_EstimatedArea.csv**: Area coverage (km2) of different habitat classes, as well as their confidence intervals, as mapped in this study. 
  
*   **HabitatMapHK_6class_ArcGISsymbology.lyrx**: Used to apply the suggested symbology in ArcGIS Pro, as shown in the figure below.

![Kwong2025_HabitatTypeMap](https://github.com/user-attachments/assets/828442cb-c287-4689-a1d9-f6e253e19921)

*   **ClassificationProbability_yyyy-yyyy.tif**: The probability values belonging to each class for every pixel. They were the intermediate products generated from the classification workflow and used to determine the final class with the highest probability and compute the forest index in this study. The sum of probabilities for all six classes is equal to 1. A scale factor of 10000 was applied to the GeoTiff files for storage convenience. 

*   **HabitatMapHK_8class_yyyy-yyyy.tif**: Based on the 6-class outputs, two more classes are added in this product, including wetland (pixel value 7) and plantation (pixel value 8), to serve as inputs for the habitat quality model.

*   **HabitatQualityHK_yyyy-yyyy.tif**: Habitat quality maps produced in this study. The pixel value is a continuous variable ranging from 0 to 1, with 1 meaning the highest habitat quality.

![Kwong2025_HabitatQualityMap](https://github.com/user-attachments/assets/bdd4dad4-7ddb-4066-a546-79e1ee0aac06)

## GIS supplementary data:

All datasets were collected and compiled from January to June 2024 and represent the conditions at that time. 

### Environmental Raster:

*   **DistanceFromCoast.tif**: Geometric distance (m) from the coastline.

*   **Elevation.tif**: Terrain height (m) from a [LiDAR-based digital terrain model](https://portal.csdi.gov.hk/geoportal/?lang=en&datasetId=cedd_rcd_1629267205233_87895).

*   **Hillfire_10periods.tif**: Hill fires occurred in each five-year period, based on burn-area products by [Chan et al. (2023)](https://www.mdpi.com/2072-4292/15/6/1489) and manual digitisation for early years.

*   **Insolation.tif**: Annual amount of incoming solar radiation (kWh/m2) computed using [SAGA GIS](https://saga-gis.sourceforge.io/saga_tool_doc/9.2.0/ta_lighting_2.html).

*   **Landslide_10periods.tif**: Landslides occurred in each five-year period, based on the [Enhanced Natural Terrain Landslide Inventory](https://portal.csdi.gov.hk/geoportal/?datasetId=cedd_rcd_1636520697377_96152&lang=en) ([Dias et al., 2009](https://hkss.cedd.gov.hk/hkss/filemanager/common/publications-resources/list-of-technical-papers/511_Dias%20et%20al%20(2009)_The%20Enhanced%20Natural%20Terrain%20Landslide%20Inventory.pdf)).

*   **Northness.tif**: Terrain aspect from 1 (due north) to -1 (due south) computed from the DTM.

*   **Precipitation.tif**: Annual precipitation (mm) (average between 1991-2020) from [Hong Kong Observatory](https://www.hko.gov.hk/en/cis/climahk.htm).

*   **Slope.tif**: Steepness (°) of the ground surface computed from the DTM.

*   **SoilCEC.tif**: Cation exchange capacity (CEC) (mmol/kg) of topsoil from [Luo et al. (2007)](https://books.google.com/books?id=ivzlMQAACAAJ).

*   **SoilOrganicMatter.tif**: Organic matter content (%) of topsoil from [Luo et al. (2007)](https://books.google.com/books?id=ivzlMQAACAAJ).

*   **Temperature.tif**: Annual mean temperature (°C) from [Morgan and Guénard (2019)](https://essd.copernicus.org/articles/11/1083/2019/).

*   **TopographicWetnessIndex.tif**: Amount of water accumulation due to topographic effects computed using [SAGA GIS](https://saga-gis.sourceforge.io/saga_tool_doc/9.2.0/ta_hydrology_15.html).

*   **Typhoon_10periods.tif**: Wind speed (km/h) estimated from [WindNinja](https://ninjastorm.firelab.org/windninja/) based on maximum hourly mean wind records associated with typhoon events in each five-year period.

*   **WindSpeed.tif**: Mean wind speed (km/h) estimated from [WindNinja](https://ninjastorm.firelab.org/windninja/) based on monthly [prevailing wind records](https://www.weather.gov.hk/en/cis/regione.htm).

![Kwong2025_EnvironmentalVariables](https://github.com/user-attachments/assets/d80674b7-9d21-4206-8b96-0e62a1913c58)

### Human Activities:

*   **BuiltupAreas_10periods_shp.zip**: Shapefile (polygons) of built-up areas, with attributes on the years of construction (estimated from topographic maps) and density (high and low). It was used as a threat factor in habitat quality mapping and variables in habitat changes.

*   **CountryParksProtectedAreas_shp.zip**: Shapefile (polygons) of protected areas (Country Parks, Special Areas, etc.), with attributes on the years of designation and revision. It was used as a protection factor in habitat quality mapping and variables in habitat changes.

*   **PollutionSource_shp.zip**: Shapefile (polygons) of pollution sources (landfills, power stations, and incineration plants), with attributes on the years of construction and closure. It was used as a threat factor in habitat quality mapping.

*   **Roads_10periods_shp.zip**: Shapefile (polylines) of roads, with attributes on the years of construction (estimated from topographic maps) and type (main and secondary). It was used as a threat factor in habitat quality mapping.

![Kwong2025_HumanActivitiesVariables](https://github.com/user-attachments/assets/afbae3bb-5674-481e-8805-b2feb4784688)

### Mapping Reference:

*   **ForestIndex_FieldCollectedReferenceData.csv**: Field survey records of habitat types which were used to evaluate the forest index variable in this study.

*   **HabitatMapHK_FieldCollectedReferenceData.csv**: Field survey records of habitat types which were used to assess the habitat mapping results in this study.

*   **HabitatMapHK_OfficeInterpretedReferenceData.csv**: Reference points where the habitat class in each period was determined through visual interpretation of the aerial photographs and other historical records. The points were used for both training and validation of the habitat maps in this study.

*   **HabitatMapHK_PredictionsFromProcedures.csv**: Prediction results for each reference point in each period when different mapping procedures were implemented.

*   **HabitatQualityHK_FieldSurveyedEcologicalValue2008.csv**: Field survey records of ecological values in 2008 which were used to evaluate the habitat quality maps in this study.

*   **LandsatHK_CrossSensorCalibrationPoints.csv**: Selected points that were assumed to remain unchanged over time and used to cross-calibrate different Landsat sensors in this study.

*   **LandsatHK_ImageMetadata.csv**: Metadata of the Landsat imagery (1,100 downloaded scenes and 607 valid scenes after pre-processing) acquired and processed in this study.

*   **Plantation_1975_1990_2008_2019.tif**: Pixels that were identified as plantations on four existing maps in different years (1975, 1990, 2008, 2019), as represented by the four layers contained in this raster file respectively. These pixels were used to help extract plantation class on the habitat map (when producing habitat quality) and denote areas with plantation activities (when modelling habitat changes) in this study.

*   **SpeciesObsHK_SpeciesChecklist.csv**: A species checklist of 7 taxa in Hong Kong (Plants, Butterflies, Birds, Reptiles, Dragonflies, Amphibians, Mammals) compiled from [AFCD](https://www.afcd.gov.hk/english/conservation/hkbiodiversity/speciesgroup/speciesgroup.html), [Hong Kong Biodiversity Information Hub](https://bih.gov.hk/en/hong-kong-species/index.html), and other secondary sources. Species of conservation concern are identified based on local assessments ([Corlett et al., 2000](https://www.researchgate.net/publication/364669045_Hong_Kong_vascular_plants_-_distribution_and_status); [Fellowes et al., 2002](https://www.researchgate.net/publication/268517493_Wild_animals_to_watch_terrestrial_and_freshwater_fauna_of_conservation_concern_in_Hong_Kong)), environmental protection laws, and national and global assessments. The checklist was used to match with the [iNaturalist observation data](https://doi.org/10.15468/dl.eexe9q) to compute biodiversity metrics at grid levels and evaluate habitat quality maps in this study.

*   **SpeciesObsHK_SynonymList.csv**: A list of species name synonyms for matching names used in iNaturalist and other secondary sources with the species checklist. It was used to pre-process the [iNaturalist observation data](https://doi.org/10.15468/dl.eexe9q) and unify the species names from different records in this study.

![Kwong2025_iNaturalistPoints](https://github.com/user-attachments/assets/9e3db740-e58f-40e7-a3e3-d26d1c6501f4)

<table><thead>
  <tr>
    <th colspan="2" rowspan="2">Taxa </th>
    <th colspan="2">Checklist</th>
    <th colspan="4">iNaturalist dataset (2018-2023)</th>
  </tr>
  <tr>
    <th>Number of species</th>
    <th>Those of conservation concerns</th>
    <th>Number of species found</th>
    <th>Those of conservation concerns</th>
    <th>Number of observations</th>
    <th>Those of species of conservation concerns</th>
  </tr></thead>
<tbody>
  <tr>
    <td colspan="2">Plants</td>
    <td>3347</td>
    <td>1139</td>
    <td>1206</td>
    <td>256</td>
    <td>29858</td>
    <td>2085</td>
  </tr>
  <tr>
    <td rowspan="6">Animals</td>
    <td>Mammals</td>
    <td>55</td>
    <td>40</td>
    <td>33</td>
    <td>21</td>
    <td>1709</td>
    <td>920</td>
  </tr>
  <tr>
    <td>Birds</td>
    <td>572</td>
    <td>230</td>
    <td>374</td>
    <td>163</td>
    <td>28313</td>
    <td>7442</td>
  </tr>
  <tr>
    <td>Reptiles</td>
    <td>90</td>
    <td>53</td>
    <td>62</td>
    <td>31</td>
    <td>4547</td>
    <td>750</td>
  </tr>
  <tr>
    <td>Amphibians</td>
    <td>25</td>
    <td>13</td>
    <td>22</td>
    <td>10</td>
    <td>3247</td>
    <td>616</td>
  </tr>
  <tr>
    <td>Butterflies</td>
    <td>245</td>
    <td>73</td>
    <td>222</td>
    <td>63</td>
    <td>37176</td>
    <td>2519</td>
  </tr>
  <tr>
    <td>Dragonflies</td>
    <td>132</td>
    <td>53</td>
    <td>96</td>
    <td>43</td>
    <td>7709</td>
    <td>643</td>
  </tr>
  <tr>
    <td colspan="2">All animals</td>
    <td>1119</td>
    <td>462</td>
    <td>809</td>
    <td>331</td>
    <td>82701</td>
    <td>12890</td>
  </tr>
  <tr>
    <td colspan="2">All taxa</td>
    <td>4466</td>
    <td>1601</td>
    <td>2015</td>
    <td>587</td>
    <td>112559</td>
    <td>14975</td>
  </tr>
</tbody></table>

## Analysis scripts:

### Part 1: Mapping Vegetation Habitats from a Satellite Image Time-Series

*   **P1_01_SearchAndDownloadFromGEE.ipynb**: Query and download all available Landsat 1-9 imagery covering the study area using Google Earth Engine. Atmospheric correction is performed if necessary.
*   **P1_02_Preprocess_part1.py**: Some basic pre-processing steps after downloading the images from cloud platform to local computer, such as mosaicking adjacent scenes and reprojecting to local coordinate system.
*   **P1_03_TopographicCorrection.R**: SCS+C topographic correction based on terrain slope, aspect, sun azimuth and sun elevation angles.
*   **P1_04_CrossSensorCal.R**: Cross-calibration of different Landsat sensors based on pseudo-invariant features, followed by computing variables for image classification.
*   **P1_05_ImageComposite.R**: Create image composites (median and standard deviation statistics) by combining all imagery acquired in the same period.
*   **P1_06_ExtractPixelValue.R**: Extract pixel values at the locations of reference points.
*   **P1_07_SpatialBlock.R**: Define spatial blocks and folds for cross-validating model predictions.
*   **P1_08_TrainRFModel.R**: Train the Random Forest model, fuse probability outputs from each image, evaluate the model accuracies with cross-validation, and create the final model for classifying the entire dataset.
*   **P1_09_TestProcedures.R**: Modify the classification procedures and re-run the Random Forest models to evaluate their impacts on the classification accuracies.
*   **P1_10_ApplyModel.R**: Apply the Random Forest model and fusion steps to all images to create the habitat map for each period.
*   **P1_11_AreaCoverage.R**: Obtain the area coverage of each class on the habitat map as well as the confidence interval of the area estimates.
*   **P1_12_CompareFieldData.R**: Assess the accuracies of the habitat maps by overlaying with field-collected points and LiDAR height information at different times.
*   **P1_13_SurvivalAnalysis.R**: Analyse the number of years required for transitioning between vegetation classes as well as the correlations between transition times and environmental variables.

### Part 2: Computing Habitat Quality Maps with Reference to Habitat Type Information

*   **P2_01_ModifyInputMap.R**: Add two more classes (wetland and plantation) to the habitat map to serve as inputs for the habitat quality model.
*   **P2_02_SensitivityParameter.R**: Compute the sensitivity parameter based on the statistical relationship between habitat changes and nearby threat factors to serve as inputs for the habitat quality model.
*   **P2_03_RunInVESTModel.R**: Execute the InVEST habitat quality model with both default and optimal half saturation constant values.
*   **P2_04_CompareFieldSurvey.R**: Evaluate the habitat quality mapping results with field-surveyed ecological value data.
*   **P2_05_CompareSpeciesObservation.R**: Pre-process and match the iNaturalist observation data, compute biodiversity metrics at grid levels, and compare with the estimated habitat quality values.
*   **P2_06_DegradationRarity.R**: Create degradation and rarity outputs from the InVEST habitat quality model and summarise their patterns over time.

### Part 3: Modelling Habitat Changes in Relation to Natural and Human Factors

*   **P3_01_CreateEnvRaster.R**: Create raster layers representing the environmental variables to serve as model inputs in this part of the study.
*   **P3_02_ForestIndex.R**: Create forest index raster by combining the probability layers generated from the Random Forest model and evaluate it with habitat classes and field-surveyed points.
*   **P3_03_ExtractPixelValue.R**: Extract the pixel values of the environmental variables for each of the 10 periods at the locations of randomly selected points.
*   **P3_04_EFA.R**: Conduct exploratory factor analysis for the variables related to the natural environment and reduce them to three factors.
*   **P3_05_SEM.R**: Conduct structural equation modelling based on the hypothesised model and computed variables.
*   **P3_06_ModelValidation.R**: Test the assumptions of individual paths in the structural equation model and compute metrics with cross-validation.
*   **P3_07_YearChange.R**: Examine the model interaction with time and path coefficients that vary against time.

---

*Last updated in February 2026*
