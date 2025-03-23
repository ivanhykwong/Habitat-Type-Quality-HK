# Spatio-Temporal Changes in Habitat Type and Quality in Hong Kong (1973-2022)

# Under Development

Supplementary materials used in the following studies:

**Kwong, I. H. Y., Lai, D. Y. F., Wong, F. K. K., & Fung, T. (2025). Spatial variations in forest succession rates revealed from multi-temporal habitat maps using Landsat imagery in subtropical Hong Kong. European Geosciences Union (EGU) General Assembly 2025, Vienna, Austria, 27 Apr–2 May 2025. https://doi.org/10.5194/egusphere-egu25-2667**

---

GIS and reference data:

*   **HabitatMapHK_6class_yyyy-yyyy.tif**: Raster data (GeoTiff format) showing the 6 habitat classes mapped in this study. Pixel values range from 1 to 6 representing woodland, shrubland, grassland, barren land, built-up area, and water respectively.
  
*   **HabitatMapHK_6class_note.xlsx**: Explanatory notes and a suggested symbology for different pixel values in the GeoTiff files.

*   **HabitatMapHK_6class_ArcGISsymbology.lyrx**: Used to apply the suggested symbology in ArcGIS Pro, as shown in the figure below.

---

Analysis codes:

**Part 1: Mapping Vegetation Habitats from a Satellite Image Time-Series**

*   **P1_1_SearchAndDownloadFromGEE.ipynb**:
*   **P1_2_Preprocess_part1.py**:
*   **P1_1_SearchAndDownloadFromGEE.ipynb**:
*   **P1_2_Preprocess_part1.py**:
*   **P1_3_TopographicCorrection.R**:
*   **P1_4_CrossSensorCal.R**:
*   **P1_5_ImageComposite.R**:
*   **P1_6_ExtractPixelValue.R**:
*   **P1_7_TrainingDataStat.R**:
*   **P1_8_TrainRFModel.R**:
*   **P1_9_TestProcedures.R**:
*   **P1_10_ApplyModel.R**:
*   **P1_11_AreaCoverage.R**:
*   **P1_12_CompareFieldData.R**:
*   **P1_13_SurvivalAnalysis.R**:

**Part 2: Computing Habitat Quality Maps with Reference to Habitat Type Information**

*   **P2_1_ModifyInputMap.R**:
*   **P2_2_SensitivityParameter.R**:
*   **P2_3_RunInVESTModel.R**:
*   **P2_4_CompareFieldSurvey.R**:
*   **P2_5_CompareSpeciesObservation.R**:
*   **P2_6_DegradationRarity.R**:

**Part 3: Modelling Habitat Changes in Relation to Natural and Human Factors**

*   **P3_1_CreateEnvRaster.R**:
*   **P3_2_ForestIndex.R**:
*   **P3_3_ExtractPixelValue.R**:
*   **P3_4_EFA.R**:
*   **P3_5_SEM.R**:
*   **P3_6_ModelValidation.R**:
*   **P3_7_YearChange.R**:

---

*Last updated in March 2025*
