## save the following InVEST model arguments as a json file
# {
#   "args": {
#     "access_vector_path": "D:/HabitatQuality/parameter/protected_y.shp",
#     "half_saturation_constant": "0.3881",
#     "lulc_bas_path": "D:/HabitatQuality/parameter/habmap_base.tif",
#     "lulc_cur_path": "D:/HabitatQuality/parameter/habmap_cur.tif",
#     "lulc_fut_path": "",
#     "results_suffix": "test1",
#     "sensitivity_table_path": "D:/HabitatQuality/parameter/sensitivity.csv",
#     "threats_table_path": "D:/HabitatQuality/parameter/threat_table_withbase.csv",
#     "workspace_dir": "D:/HabitatQuality/workspace_comparefirst"
#   },
#   "invest_version": "3.14.1",
#   "model_name": "natcap.invest.habitat_quality"
# }

library(terra)

# function to run invest model
runinvestmodel_comparefirst <- function(){
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
  
  for (i in 2:10){
    setwd("D:/HabitatQuality/parameter")
    targetyear <- yearlist[i] # current year
    
    # write input map and threats for current year
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
    
    writeRaster(habmap, "habmap_cur.tif", overwrite=TRUE)
    writeRaster(highbuild_y_ras, "highbuild_cur.tif", overwrite=TRUE)
    writeRaster(lowbuild_y_ras, "lowbuild_cur.tif", overwrite=TRUE)
    writeRaster(mainroad_y_ras, "mainroad_cur.tif", overwrite=TRUE)
    writeRaster(secroad_y_ras, "secroad_cur.tif", overwrite=TRUE)
    writeRaster(pollution_y_ras, "pollution_cur.tif", overwrite=TRUE)
    writeVector(protected_y, "protected_y.shp", overwrite=TRUE)
    
    # write input map and threats for first year
    targetyear <- yearlist[1] # base year
    habmap <- rast(maplist[1])
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
    
    writeRaster(habmap, "habmap_base.tif", overwrite=TRUE)
    writeRaster(highbuild_y_ras, "highbuild_base.tif", overwrite=TRUE)
    writeRaster(lowbuild_y_ras, "lowbuild_base.tif", overwrite=TRUE)
    writeRaster(mainroad_y_ras, "mainroad_base.tif", overwrite=TRUE)
    writeRaster(secroad_y_ras, "secroad_base.tif", overwrite=TRUE)
    writeRaster(pollution_y_ras, "pollution_base.tif", overwrite=TRUE)

    # link to InVEST program exe
    investexe <- "C:/Users/xxx/AppData/Local/Programs/InVEST 3.14.1 Workbench/resources/invest/invest.exe"
    system2(investexe, args="run hq -d D:/HabitatQuality/parameter/arguments_comparefirst.json")

    # save the InVEST output files according to years
    setwd("D:/HabitatQuality/workspace_comparefirst")
    rarity_c <- rast("rarity_c_test1.tif")
    rarity_c <- mask(rarity_c, habmap==9, maskvalues=1)
    writeRaster(rarity_c, paste0("rarity_c_y",yearlist[i],".tif"))
    
    # remove unnecessary files
    file.remove(c("deg_sum_c_test1.tif", "quality_c_test1.tif",
                  "deg_sum_b_test1.tif", "quality_b_test1.tif",
                  "rarity_c_test1.tif"))
    setwd("D:/HabitatQuality/parameter")
    protected_y_shp <- list.files(getwd(), "protected_y")
    file.remove(c("habmap_cur.tif", "highbuild_cur.tif", "lowbuild_cur.tif", 
                  "mainroad_cur.tif", "secroad_cur.tif", "pollution_cur.tif",
                  "habmap_base.tif", "highbuild_base.tif", "lowbuild_base.tif", 
                  "mainroad_base.tif", "secroad_base.tif", "pollution_base.tif",
                  protected_y_shp))
  }
}
runinvestmodel_comparefirst()


# summarise rarity

maplist <- list.files("D:/HabitatQuality/inputmap", "8class.tif$", full.names=TRUE)
maplist <- maplist[2:10]
rarlist <- list.files("D:/HabitatQuality/workspace_comparefirst", "rarity_c_y", full.names=TRUE)

for (i in 1:9){
  habmap <- rast(maplist[i])  # for every period
  rarmap <- rast(rarlist[i])
  dfm <- zonal(rarmap, habmap, "mean", na.rm=TRUE) # mean rarity of each class
  colnames(dfm) <- c("class","rar_mean")
  dfs <- zonal(rarmap, habmap, "sd", na.rm=TRUE) # sd of rarity of each class
  colnames(dfs) <- c("class","rar_sd")
  df1 <- merge(dfm, dfs, by="class")
  df1 <- rbind(df1, data.frame(class="all",
                               rar_mean=global(rarmap,"mean",na.rm=TRUE)[[1]],
                               rar_sd=global(rarmap,"sd",na.rm=TRUE)[[1]]))
  colnames(df1) <- c("class",paste0("rar_mean",i),paste0("rar_sd",i))
  if (i==1){df <- df1}
  else {df <- merge(df, df1, by="class")}
}
write.csv(df,"rarity.csv")


# summarise degradation

maplist <- list.files("D:/HabitatQuality/inputmap", "8class.tif$", full.names=TRUE)
deglist <- list.files("D:/HabitatQuality/workspace", "deg_sum_c_y", full.names=TRUE)

for (i in 1:10){
  habmap <- rast(maplist[i])
  degmap <- rast(deglist[i])
  dfm <- zonal(degmap, habmap, "mean", na.rm=TRUE)
  colnames(dfm) <- c("class","deg_mean")
  dfs <- zonal(degmap, habmap, "sd", na.rm=TRUE)
  colnames(dfs) <- c("class","deg_sd")
  df1 <- merge(dfm, dfs, by="class")
  df1 <- rbind(df1, data.frame(class="all",
                               deg_mean=global(degmap,"mean",na.rm=TRUE)[[1]],
                               deg_sd=global(degmap,"sd",na.rm=TRUE)[[1]]))
  colnames(df1) <- c("class",paste0("deg_mean",i),paste0("deg_sd",i))
  if (i==1){df <- df1}
  else {df <- merge(df, df1, by="class")}
}
write.csv(df,"degradation.csv")
