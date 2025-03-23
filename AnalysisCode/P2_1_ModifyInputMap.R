library(terra)
library(reshape2)

# Add wetland class

setwd("D:/HabitatQuality")
mapdir <- "D:/landsat/outputmap/"
maplist <- list.files(mapdir, ".tif$")
imgdir <- "D:/landsat/composite/"
imglist <- list.files(imgdir, ".tif$")
outfolder <- "D:/HabitatQuality/inputmap"

for (i in 1:10){
  habmap_name <- maplist[i]
  habmap <- rast(file.path(mapdir,habmap_name))
  img <- rast(file.path(imgdir,imglist[i]))
  crs(img) <- crs(habmap)
  
  ndvi <- (img[["NIR"]]-img[["R"]])/(img[["NIR"]]+img[["R"]])
  ndvi <- crop(ndvi, habmap, mask=TRUE)
  ndwi <- (img[["G"]]-img[["NIR"]])/(img[["G"]]+img[["NIR"]])
  ndwi <- crop(ndwi, habmap, mask=TRUE)
  
  # water rectify to built-up area
  # get "inland" water
  habmap_w <- ifel(is.na(habmap), 6, habmap)
  habmap_w <- ifel(habmap_w!=6, NA, habmap_w)
  habmap_w <- patches(habmap_w)
  habmap_w <- ifel(is.na(habmap_w), habmap, ifel(habmap_w==1,6,NA))
  habmap_w <- focal(habmap_w,w=5,fun="mean",na.policy="omit",na.rm=TRUE)
  # fill "inland" water using focal
  while (global(is.na(habmap_w),"sum")[[1]] > 0){
    habmap_w <- focal(habmap_w,focalMat(habmap_w, 30, "circle"),
                      fun="mean",na.policy="only",na.rm=TRUE)
  }
  habmap[habmap==6 & ndwi<0 & habmap_w==5] <- 5
  
  # remove water on reclaimed land to NA
  reclaim <- vect("D:/HabitatQuality/referencedata/reclaimedland.shp")
  reclaim <- rasterize(reclaim, habmap, background=0)
  habmap <- ifel(habmap==6 & reclaim==1, NA, habmap)
  
  # grassland rectify to wetland
  # low-lying brackish wetland: close to water, far from dry built-up
  elev <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")
  elev <- crop(elev, habmap, mask=TRUE)
  # water pixel > built-up pixel in surrounding 120m
  # reclass water to 1, built-up to -1, other to 0, focal sum -> positive
  habmap_wb <- habmap
  habmap_wb[habmap_wb <= 4] <- 0
  habmap_wb[habmap_wb == 6] <- 1
  habmap_wb[habmap_wb == 5] <- -1
  habmap_wb_f <- focal(habmap_wb, focalMat(habmap_wb, 120, "circle"),
                       fun="sum", na.policy="omit", na.rm=TRUE)
  habmap[habmap==3 & elev<5 & habmap_wb_f>0] <- 7
  
  # paddy wetland
  # paddy in 1966, check TWI (wetness of location), check current NDVI (vegetated)
  twi <- rast("D:/HabitatChange/GISData/NaturalEnv/TWI.tif")
  twi <- crop(twi, habmap, mask=TRUE)
  paddy <- vect("D:/HabitatQuality/referencedata/Paddy_1966.shp")
  paddy_ras <- mask(habmap, paddy, touches=FALSE)
  paddy_ras <- ifel(is.na(paddy_ras),0,1)
  habmap[habmap==3 & paddy_ras==1 & twi>8 & ndvi>0.6] <- 7
  
  outname <- gsub(".tif", "_7class.tif", habmap_name)
  outname <- gsub("outmap", "habmap", outname)
  writeRaster(habmap, file.path(outfolder,outname))
}

# check wetland accuracy

setwd("D:/HabitatQuality")
habmap2010 <- rast("inputmap/habmap_2008-2012_7class.tif")
habmap2020 <- rast("inputmap/habmap_2018-2022_7class.tif")
pt2010 <- read.csv("wetland/pt_2010_wetland.csv")
pt2020 <- read.csv("wetland/pt_2020_wetland.csv")

checkwetland <- function(ras, df, bufdist=60){
  # create point from df
  pt <- vect(df, geom=c("Rect_X","Rect_Y"), crs=crs(ras))
  pt_buff <- buffer(pt, bufdist)
  # extract value from ras
  pt_buff_ex <- extract(ras, pt_buff, fun=table)
  pt_buff_ex <- melt(pt_buff_ex, id.vars="ID", variable.name="mapclass",value.name="count")
  mapdf <- pt_buff_ex[pt_buff_ex$count>0,]
  # merge with refclass
  mapdf$refclass <- 7
  mapdf$classequal <- mapdf$mapclass == mapdf$refclass
  mapdf <- mapdf[order(-mapdf$classequal, -mapdf$count),]
  mapdf <- mapdf[!duplicated(mapdf$ID),]
  mapdf$ID <- as.numeric(mapdf$ID)
  mapdf <- mapdf[order(mapdf$ID),]
  print(sum(mapdf$classequal)/length(mapdf$classequal))
  return(mapdf)
}
mapdf2010 <- checkwetland(habmap2010, pt2010)
mapdf2020 <- checkwetland(habmap2020, pt2020)


# Add plantation class

# union 3 plantation layers from previous habitat mapping (1990,2010,2020)
# indicate locations where plantation activities have been conducted
# and plantation forests have successfully established
# plantation hardly be changed even after decades
# intersect with plantation layer & habitat map == woodland or shrubland -> plantation

# Prepare plantation layer ----

setwd("D:/HabitatChange/GISData/Plantation")
p1990 <- vect("Plantation_1990.shp")
p2010 <- rast("Plantation_2008.tif")
p2020 <- rast("Plantation_2019.tif")
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")

p1990_ras <- rasterize(p1990, tempras, background=0)
p2010_ras <- terra::resample(p2010, tempras, method="sum")
p2010_ras <- p2010_ras>15
p2020_ras <- terra::resample(p2020, tempras, method="mode")
plantras <- p1990_ras + p2010_ras + p2020_ras
plantras[is.na(plantras)] <- 0
plantras[plantras>=1] <- 1
names(plantras) <- "plantation"
writeRaster(plantras, "Plantation_union3periods.tif")


# modify habitat map ----

habmapdir <- "D:/HabitatQuality/inputmap"
setwd(habmapdir)
habmaplist <- list.files(habmapdir, "7class.tif$")

for (i in 1:10){
  habmapname <- habmaplist[i]
  habmap <- rast(habmapname)
  plantras1 <- crop(plantras, habmap, mask=TRUE)
  newmap <- ifel(plantras1==1 & habmap%in%c(1,2), 8, habmap)  # plantation==TRUE & habitatmap==wood or shrub
  writeRaster(newmap, gsub("7class", "8class", habmapname))
}
