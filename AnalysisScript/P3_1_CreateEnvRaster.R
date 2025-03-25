library(sf)
library(terra)

#### Insolation ####

# Potential Incoming Solar Radiation Tool in SAGA GIS
# C:/Users/xxx/Desktop/saga-9.2.0_x64/saga_cmd.exe ta_lighting 2 
# -GRD_DEM dsm2010_30m_fill.tif -GRD_SVF "Sky View Factor.sg-grd-z" -PERIOD 2 
# -HOUR_RANGE_MIN 0 -HOUR_RANGE_MAX 24 -DAY 2020-01-01 -DAY_STOP 2020-01-31 
# -GRD_DIRECT DirectInsolation_Jan.sg-grd-z -GRD_DIFFUS DiffusInsolation_Jan.sg-grd-z
# -GRD_TOTAL TotalInsolation_Jan.sg-grd-z

# Compare to HKO solar radiation records (Table 6 in https://www.hko.gov.hk/en/cis/normal/1991_2020/normals.htm)

library(terra)

setwd("D:/HabitatChange/GISData/insolation")
raslist <- paste0("TotalInsolation_", month.abb, ".sg-grd-z")
raslist <- lapply(raslist, rast)
hko <- read.csv("SunshinePercentage_KingsPark.csv")
sunper <- hko$SunshinePercent[1:12]/100

# multiply sunshine percentage
raslist_m <- raslist
for (i in 1:12){
  raslist_m[[i]] <- raslist_m[[i]]*sunper[i]
}
ras_m <- rast(raslist_m)

# extract value for King's Park
kingspark <- data.frame(x=835821,y=819181)
est <- extract(ras_m, kingspark, method="bilinear")
est <- array(unlist(est))[2:13]

tru <- hko$Radiation[1:12]
monthday <- c(31,28,31,30,31,30,31,31,30,31,30,31)
tru <- tru*monthday/3.6 # 1 kWh = 3.6 MJ
cor(est, tru) # 0.9589472
a <- lm(tru~est)
summary(a) # R-squared:  0.9196
plot(est, tru)
sqrt(mean((est-tru)^2)) # rmse

# sum up to annual
annualras <- sum(ras_m) # AnnualInsolation
writeRaster(annualras, "D:/HabitatChange/GISData/NaturalEnv/insolation.tif") 
plot(annualras)


#### Wind ####

# check wind speed
# https://www.hko.gov.hk/tc/cis/regione.htm (mean wind speed)
# https://www.hko.gov.hk/tc/cis/stn.htm (location)

ras <- rast("D:/HabitatChange/GISData/NaturalEnv/windspeed.tif")
plot(ras)
df <- read.csv("HKONWS.csv")  # csv file with 16 wind stations (Easting,Northing) and their wind speed
pt <- vect(df, geom=c("Easting", "Northing"), crs=crs(ras))
est <- extract(ras, pt, method="bilinear")[,2]
tru <- df[,"WindSpeed"]
cor.test(est, tru) # r=0.907884 n=16 df=14, p-value=1.179e-06
a <- lm(tru~est)
summary(a) # R-squared:  0.8243
plot(est, tru)
sqrt(mean((est-tru)^2)) # rmse


#### Distance from buildings ####

buildshp <- vect("D:/HabitatQuality/threat/Buildings.shp")
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  bs <- buildshp[buildshp$AppearYear<=y,]
  bs_ras <- rasterize(bs, tempras)
  bs_dist <- distance(bs_ras)
  names(bs_dist) <- paste0("c",y)
  if (i==1) {ras <- bs_dist}
  else {ras <- c(ras, bs_dist)}
}
writeRaster(ras, "D:/HabitatChange/GISData/HumanAct/DistfromBuildings_10periods.tif")


#### Country park ####

shp <- vect("D:/HabitatQuality/protection/ProtectedAreas.shp")
shp <- shp[shp$Type=="Country Parks",]
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  cp <- shp[(shp$AppearYear<=y)&(shp$EndYear>=y),]
  cp_ras <- rasterize(cp, tempras)
  cp_ras[is.na(cp_ras)] <- 0
  names(cp_ras) <- paste0("c",y)
  if (i==1) {ras <- cp_ras}
  else {ras <- c(ras, cp_ras)}
}
writeRaster(ras, "D:/HabitatChange/GISData/HumanAct/CountryParksSpecialAreas_10periods.tif")


#### Plantation #### 

# extract from 1975, 1990, 2010, 2020 map
setwd("D:/HabitatChange/GISData/Plantation")
p1975 <- vect("Plantation_1975.shp")
p1990 <- vect("Plantation_1990.shp")
p2010 <- rast("Plantation_2008.tif")
p2020 <- rast("Plantation_2019.tif")
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")

p1975_ras <- rasterize(p1975, tempras, background=0)
p1990_ras <- rasterize(p1990, tempras, background=0)
p2010_ras <- terra::resample(p2010, tempras, method="sum")
p2010_ras <- p2010_ras>15
p2020_ras <- terra::resample(p2020, tempras, method="mode")

# combine plantation pixels in all 4 periods
plantras <- p1975_ras + p1990_ras + p2010_ras + p2020_ras
plantras[is.na(plantras)] <- 0
plantras[plantras>=1] <- 1
names(plantras) <- "plantation"
writeRaster(plantras, "D:/HabitatChange/GISData/HumanAct/Plantation_union4periods.tif")

# overlay with habitat map to find plantation year
mapdir <- "D:/landsat/outputmap/"
maplist <- list.files(mapdir, ".tif$", full.names=TRUE)
maplist <- lapply(maplist, rast)
habmap <- rast(maplist)
plantras1 <- crop(plantras, habmap, mask=TRUE)
habmap <- ifel(plantras1==1 & habmap%in%c(1,2), 1, 0) # plantation==TRUE & habitatmap==wood or shrub
plantyear <- which.lyr(habmap) # first layer in the input that is TRUE
plantyear <- plantyear*5+1970
writeRaster(plantyear, "D:/HabitatChange/GISData/Plantation/plantyear.tif")

# create raster layer showing plantation pixels in current year
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  plantyear_y <- plantyear==y   # current year
  plantyear_y[is.na(plantyear_y)] <- 0
  names(plantyear_y) <- paste0("c",y)
  if (i==1) {ras <- plantyear_y}
  else {ras <- c(ras, plantyear_y)}
}
writeRaster(ras, "D:/HabitatChange/GISData/HumanAct/Plantation_10periods_current.tif")

# create raster layer showing plantation pixels in previous years
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  plantyear_y <- plantyear<y   # previous year
  plantyear_y[is.na(plantyear_y)] <- 0
  names(plantyear_y) <- paste0("c",y)
  if (i==1) {ras <- plantyear_y}
  else {ras <- c(ras, plantyear_y)}
}
writeRaster(ras, "D:/HabitatChange/GISData/HumanAct/Plantation_10periods_previous.tif")


#### Distance from woodland ####

mapdir <- "D:/landsat/outputmap/"
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")*0
maplist <- list.files(mapdir, ".tif$", full.names=TRUE)
maplist <- lapply(maplist, rast)
get.dist.from.wood <- function(ras){
  wood <- ifel(ras==1, 1, NA)
  wood[is.na(wood)] <- 0
  wood <- focal(wood, 3, fun="modal") # 3x3 majority filter
  wood[wood==0] <- NA
  woodline <- as.lines(as.polygons(wood)) # from boundary line
  distwood <- distance(tempras, woodline)
  return(distwood)
}
maplist <- lapply(maplist, get.dist.from.wood)
ras <- rast(maplist)
names(ras) <- paste0("c",seq(1975,2020,5))
writeRaster(ras, "D:/HabitatChange/GISData/WoodSpatial/DistfromWoodland_10periods.tif")


#### Area of nearest wood patch ####

mapdir <- "D:/landsat/outputmap/"
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")*0
maplist <- list.files(mapdir, ".tif$", full.names=TRUE)
maplist <- lapply(maplist, rast)
get.nearest.size <- function(ras){
  wood <- ifel(ras==1, 1, NA)
  wood[is.na(wood)] <- 0
  wood <- focal(wood, 3, fun="modal") # 3x3 majority filter
  wood[wood==0] <- NA
  y <- patches(wood)
  rz <- zonal(cellSize(y, unit="ha"), y, sum, as.raster=TRUE) # get area
  while (global(is.na(rz),"sum")[[1]] > 0) {
    rz <- focal(rz, 3, fun="mean", na.policy="only", na.rm=TRUE) # assign to nearest
  }
  return(rz)
}
maplist <- lapply(maplist, get.nearest.size)
ras <- rast(maplist)
names(ras) <- paste0("c",seq(1975,2020,5))
writeRaster(ras, "D:/HabitatChange/GISData/WoodSpatial/AreaNearestWoodland_10periods.tif")


#### Hill fire ####

setwd("D:/HabitatChange/GISData/Hillfire")
shp <- st_read("Fire_polygons_hmbl_nbr_veg_hk_bound_final.shp") # downloaded from Figshare of Chan et al. 2023
shp$year <- as.character(shp$fire_date)
shp$year <- substring(shp$year, 1, 4)
shp$year <- as.integer(shp$year)
write_sf(shp, "hillfire_1987-2020.shp")

v1 <- vect("hillfire_1973-1986.shp") # manually digitized from Landsat imagery
v2 <- vect("hillfire_1987-2020.shp")
firevect <- rbind(v1,v2)
firevect$yearclass <- round(firevect$year/5)*5
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  fv <- firevect[firevect$yearclass==y,]
  fv_ras <- rasterize(fv, tempras)
  fv_ras <- ifel(is.na(fv_ras),0,1)
  names(fv_ras) <- paste0("c",y)
  if (i==1) {ras <- fv_ras}
  else {ras <- c(ras, fv_ras)}
}
writeRaster(ras, "D:/HabitatChange/GISData/Disturbance/Hillfire_10periods.tif")


#### Landslide ####

setwd("D:/HabitatChange/GISData/Landslide")
shp <- vect("EnhancedNaturalTerrainLandslideInventory_FGDB/ENTLI (Up to Year 2019).gdb", layer="ENTLI_Trail")
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")

shp <- shp[,c("M_WIDTH","YEAR_1")]
shp$YEAR_1 <- as.integer(shp$YEAR_1)
shp$yearclass <- round(shp$YEAR_1/5)*5
yc <- seq(1975,2020,5)
for (i in 1:length(yc)){
  y <- yc[i]
  shp_y <- shp[shp$yearclass==y]
  shp_y <- buffer(shp_y, width=shp_y$M_WIDTH/2)  # create buffer using half of reported landslide width
  ras_y <- rasterize(shp_y, tempras, touches=TRUE)
  ras_y[is.na(ras_y)] <- 0
  if (i==1) {ras <- ras_y}
  else {ras <- c(ras, ras_y)}
}
writeRaster(ras, "D:/HabitatChange/GISData/Disturbance/Landslide_10periods.tif")


#### Habitat quality ####

hqmap_list <- list.files("D:/HabitatQuality/workspace", 
                         pattern="quality_c_y", full.names=TRUE)  # HQ maps created in P2
hqmap <- lapply(hqmap_list, rast)
hqmap <- rast(hqmap)  # combine all 10 maps
names(hqmap) <- paste0("c",seq(1975,2020,5))
writeRaster(hqmap, "D:/HabitatChange/HabitatPattern/HabitatQuality_10periods.tif")


#### valid pixels (class 1-4 in all periods) ####

mapdir <- "D:/landsat/outputmap/"
maplist <- list.files(mapdir, ".tif$", full.names=TRUE)
maplist <- lapply(maplist, rast)
ras <- rast(maplist)  # combine all 10 habitat maps
names(ras) <- paste0("c",seq(1975,2020,5))
ras1 <- ifel(ras<=4,1,0)  # class <= 4
ras1 <- min(ras1)    # min value across 10 layers
ras1[ras1==0] <- NA  # set 0 (class > 4 for any 1 year) to NA
writeRaster(ras1, "D:/HabitatChange/HabitatPattern/validpixels.tif")

