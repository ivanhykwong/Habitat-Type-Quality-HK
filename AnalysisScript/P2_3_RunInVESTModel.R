## save the following InVEST model arguments as a json file
# {
#   "args": {
#     "access_vector_path": "D:/HabitatQuality/parameter/protected_y.shp",
#     "half_saturation_constant": "0.05",
#     "lulc_bas_path": "",
#     "lulc_cur_path": "D:/HabitatQuality/parameter/habmap.tif",
#     "lulc_fut_path": "",
#     "results_suffix": "test1",
#     "sensitivity_table_path": "D:/HabitatQuality/parameter/sensitivity.csv",
#     "threats_table_path": "D:/HabitatQuality/parameter/threat_table.csv",
#     "workspace_dir": "D:/HabitatQuality/workspace"
#   },
#   "invest_version": "3.14.1",
#   "model_name": "natcap.invest.habitat_quality"
# }

# first trial using default Half-saturation Constant 0.05

library(terra)
runinvestmodel <- function(firsttrial){
  # Load habitat maps and threat source shapefile
  setwd("D:/HabitatQuality/")
  build_all <- vect("threat/Buildings.shp")
  highbuild_all <- build_all[build_all$Density=="High"]
  lowbuild_all <- build_all[build_all$Density=="Low"]
  road_all <- vect("threat/Roads.shp")
  mainroad_all <- road_all[road_all$Class=="Main"]
  secroad_all <- road_all[road_all$Class=="Secondary"]
  pollution_all <- vect("threat/Pollution.shp")
  protected_all <- vect("protection/ProtectedAreas.shp")
  yearlist <- seq(1975,2020,5)
  maplist <- list.files("inputmap", "8class.tif$", full.names=TRUE)
  
  for (i in 1:10){  # for each period
    # subset input map and threats for this period
    targetyear <- yearlist[i]
    habmap <- rast(maplist[i])
    habmap[is.na(habmap)] <- 9  # NA to class 9
    highbuild_y <- highbuild_all[highbuild_all$AppearYear<=targetyear]
    highbuild_y_ras <- rasterize(highbuild_y, habmap, background=0)
    lowbuild_y <- lowbuild_all[lowbuild_all$AppearYear<=targetyear]
    lowbuild_y_ras <- rasterize(lowbuild_y, habmap, background=0)
    mainroad_y <- mainroad_all[mainroad_all$AppearYear<=targetyear]
    mainroad_y_ras <- rasterize(mainroad_y, habmap, background=0)
    secroad_y <- secroad_all[secroad_all$AppearYear<=targetyear]
    secroad_y_ras <- rasterize(secroad_y, habmap, background=0)
    pollution_y <- pollution_all[pollution_all$AppearYear<=targetyear & pollution_all$EndYear>=targetyear]
    pollution_y_ras <- rasterize(pollution_y, habmap, background=0)
    protected_y <- protected_all[protected_all$AppearYear<=targetyear & protected_all$EndYear>=targetyear]
    
    # write input map and threats as raster files
    setwd("D:/HabitatQuality/parameter")
    writeRaster(habmap, "habmap.tif", overwrite=TRUE)
    writeRaster(highbuild_y_ras, "highbuild_y.tif", overwrite=TRUE)
    writeRaster(lowbuild_y_ras, "lowbuild_y.tif", overwrite=TRUE)
    writeRaster(mainroad_y_ras, "mainroad_y.tif", overwrite=TRUE)
    writeRaster(secroad_y_ras, "secroad_y.tif", overwrite=TRUE)
    writeRaster(pollution_y_ras, "pollution_y.tif", overwrite=TRUE)
    writeVector(protected_y, "protected_y.shp", overwrite=TRUE)
    
    # link to InVEST program exe
    investexe <- "C:/Users/xxx/AppData/Local/Programs/InVEST 3.14.1 Workbench/resources/invest/invest.exe"
    if (firsttrial==TRUE){
      # link to argument json file
      system2(investexe, args = "run hq -d D:/HabitatQuality/parameter/arguments.json")
    } else {
      # use another argument json file that stores the optimal half_saturation_constant
      system2(investexe, args = "run hq -d D:/HabitatQuality/parameter/arguments_k03881.json")
    }
    
    # save the InVEST output files according to years
    setwd("D:/HabitatQuality/workspace")
    deg_sum_c <- rast("deg_sum_c_test1.tif")
    quality_c <- rast("quality_c_test1.tif")
    deg_sum_c <- mask(deg_sum_c, habmap==9, maskvalues=1)
    writeRaster(deg_sum_c, paste0("deg_sum_c_y",targetyear,".tif"))
    quality_c <- mask(quality_c, habmap==9, maskvalues=1)
    writeRaster(quality_c, paste0("quality_c_y",targetyear,".tif"))
    file.remove("deg_sum_c_test1.tif")
    file.remove("quality_c_test1.tif")
  }
}
runinvestmodel(firsttrial=TRUE)  # run the first trial


# inspect the first trial results and get optimal Half-saturation Constant
# (half of the highest grid cell degradation value on the landscape)

deg_sum_c_list <- list.files("D:/HabitatQuality/workspace", "deg_sum_c_y.+.tif$", full.names=TRUE)
getrasmax <- function(rasname){
  ras <- rast(rasname)
  return(global(ras, "max", na.rm=TRUE)[[1]])
}
deg_sum_c_list <- lapply(deg_sum_c_list, getrasmax)
print(paste("k value:", max(unlist(deg_sum_c_list))/2)) # 0.3881 in this study


# run the InVEST model again after modify json file
runinvestmodel(firsttrial=FALSE)


# summary statistics for the habitat quality maps
library(reshape2)
hqmap_list <- list.files("D:/HabitatQuality/workspace", pattern="quality_c_y", full.names=TRUE)
hqmap <- lapply(hqmap_list, rast)
hqmap <- rast(hqmap)
rc <- classify(hqmap, c(0, 0.2, 0.4, 0.6, 0.8, 1))
rc_f <- freq(rc)
rc_f <- dcast(rc_f, layer~value, value.var="count")
rc_f[,2:6] <- rc_f[,2:6]/rowSums(rc_f[,2:6])

hqmap_stat1 <- cbind(global(hqmap, c("mean","sd"), na.rm=TRUE),
                     global(hqmap, quantile, probs=c(0.25,0.5,0.75), na.rm=TRUE))
hqmap_stat <- cbind(rc_f, hqmap_stat1)
write.csv(hqmap_stat, "D:/HabitatQuality/HabitatQualityMap_Stat.csv")

