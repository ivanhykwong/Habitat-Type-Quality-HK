library(terra)

# select random points (100000 from valid pixel)

setwd("D:/HabitatChange")
vp <- rast("HabitatPattern/validpixels.tif")
soil1 <- crop(rast("GISData/NaturalEnv/SoilCEC.tif"),vp)
soil2 <- crop(rast("GISData/NaturalEnv/SoilOrganicMatter.tif"),vp)
vp <- vp*soil1*soil2*0+1
# global(vp,sum, na.rm=T)  # pixel count: 813633
set.seed(2023)
pt1 <- spatSample(vp, 100000, "random", as.df=TRUE, xy=TRUE, na.rm=TRUE)
yc <- paste0("c",seq(1975,2020,5))
for (i in 1:10){          # duplicate the same 100000 points for 10 periods
  pt2 <- pt1[,c("x","y")]
  pt2$year <- i
  pt2$year1 <- yc[i]
  if (i==1) {pt <- pt2}
  else {pt <- rbind(pt, pt2)}
}
pt$ID <- 1:nrow(pt)
pt_df <- pt
pt <- vect(pt, geom=c("x","y"), crs=crs(vp), keepgeom=TRUE)
writeVector(pt, "SEM/randompt.shp")

# extract values

setwd("D:/HabitatChange")
pt <- vect("SEM/randompt.shp")
pt_df <- as.data.frame(pt)
pt_df$yrsinb <- (pt_df$year-1)*5

# variable is single layer and constant throughout periods
addcons <- function(pt_df, varname, raspath){  
  ras <- rast(raspath)
  v <- extract(ras, pt, ID=FALSE)
  colnames(v) <- varname
  return(cbind(pt_df,v))
}
# variable has 10 layers representing each period
addvary10 <- function(pt_df, varname, raspath){
  ras <- rast(raspath)
  for (i in 1:10){
    ras_y <- ras[[i]]
    pt_y <- pt[pt$year==i,]
    v_y <- extract(ras_y, pt_y, ID=FALSE)
    colnames(v_y) <- varname
    if (i==1) {v <- v_y}
    else {v <- rbind(v, v_y)}
  }
  return(cbind(pt_df,v))
}
# variable has 10 layers, extract the value for previous period
addvarypre <- function(pt_df, varname, raspath){
  ras <- rast(raspath)
  v <- data.frame(rep(NA, nrow(pt[pt$year==1,]))) # create 10000 na rows for year 1
  colnames(v) <- varname
  for (i in 1:9){
    ras_y <- ras[[i]]
    pt_y <- pt[pt$year==(i+1),]  # match raster 1 and year 2
    v_y <- extract(ras_y, pt_y, ID=FALSE)
    colnames(v_y) <- varname
    v <- rbind(v, v_y)
  }
  return(cbind(pt_df,v))
}
# variable has 10 layers, extract the value for next period
addvarynext <- function(pt_df, varname, raspath){
  ras <- rast(raspath)
  v1 <- data.frame(rep(NA, nrow(pt[pt$year==10,]))) # create 10000 na rows for year 10
  colnames(v1) <- varname
  for (i in 1:9){
    ras_y <- ras[[i+1]]
    pt_y <- pt[pt$year==i,]  # match raster 2 and year 1
    v_y <- extract(ras_y, pt_y, ID=FALSE)
    colnames(v_y) <- varname
    if (i==1) {v <- v_y}
    else {v <- rbind(v, v_y)}
  }
  v <- rbind(v, v1)
  return(cbind(pt_df,v))
}

# List of variables and sources
pt_df_data <- pt_df |>
  addcons("Elev","GISData/NaturalEnv/dtm.tif") |>
  addcons("Slope","GISData/NaturalEnv/slope.tif") |>
  addcons("North","GISData/NaturalEnv/northness.tif") |>
  addcons("TWI","GISData/NaturalEnv/TWI.tif") |>
  addcons("Insol","GISData/NaturalEnv/insolation.tif") |>  
  addcons("Temp","GISData/NaturalEnv/temperature.tif") |>  
  addcons("Precip","GISData/NaturalEnv/precipitation.tif") |>
  addcons("Wind","GISData/NaturalEnv/windspeed.tif") |>
  addcons("DistCoast","GISData/NaturalEnv/DistfromCoast.tif") |>
  addcons("CEC","GISData/NaturalEnv/SoilCEC.tif") |>
  addcons("SOM","GISData/NaturalEnv/SoilOrganicMatter.tif") |>

  addvary10("ProxBuilt","GISData/HumanAct/DistfromBuildings_10periods.tif") |>
  addvary10("CounPark","GISData/HumanAct/CountryParksSpecialAreas_10periods.tif") |>
  addvary10("PlantCur","GISData/HumanAct/Plantation_10periods_current.tif") |>
  addvary10("PlantPre","GISData/HumanAct/Plantation_10periods_previous.tif") |>

  addvarypre("DistWood","GISData/WoodSpatial/DistfromWoodland_10periods.tif") |>
  addvarypre("AreaWood","GISData/WoodSpatial/AreaNearestWoodland_10periods.tif") |>  

  addvary10("FireCur","GISData/Disturbance/Hillfire_10periods.tif") |>
  addvary10("Landslide","GISData/Disturbance/Landslide_10periods.tif") |>
  addvary10("Typhoon","GISData/Disturbance/Typhoon_10periods.tif") |>
 
  addvary10("ForIndCur","HabitatPattern/ForestIndex_10periods.tif") |>
  addvary10("HabQuaCur","HabitatPattern/HabitatQuality_10periods.tif") |>
  addvarynext("ForIndNext","HabitatPattern/ForestIndex_10periods.tif") |>
  addvarynext("HabQuaNext","HabitatPattern/HabitatQuality_10periods.tif")

write.csv(pt_df_data, "SEM/pt_df_data.csv", row.names=FALSE)

# correlation and summary
pt_df_data_sub <- pt_df_data[pt_df_data$year%in%(2:9),6:30] 
cor(pt_df_data_sub)
summary(pt_df_data_sub)

