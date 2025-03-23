library(raster)
library(sf)
library(dplyr)

imgfolder <- "D:/landsat/preprocess/variables"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "tif$")

trainpt <- read.csv("D:/landsat/train/trainpt.csv")
trainpt_sf <- st_as_sf(trainpt, coords=c("x", "y"), crs = crs(raster(imglist[1])))

extractdata <- function(imgname){
  img <- brick(imgname)
  # check sensor and assign bands
  sensorname <- substr(imgname, 20, 21)
  addvar <- c("ndvi","gndvi","glcm_mean","glcm_con")
  if (sensorname == "MS"){
    names(img) <- append(c("G","R","RE","NIR"), addvar)
  } else if (sensorname == "TM" | sensorname == "ET"){
    names(img) <- append(c("B","G","R","NIR","SW1","SW2"), addvar)
  } else if (sensorname == "OL"){
    names(img) <- append(c("C","B","G","R","NIR","SW1","SW2"), addvar)
  } 
  imgdatadf <- extract(img, trainpt_sf, df=TRUE)
  imgdatadf <- cbind(imgdatadf, "PointID" = trainpt_sf$PointID)
  # extract image year class
  imgyear <- as.numeric(substr(imgname, 9, 12))
  yearclass <- paste0("c", round(imgyear/5)*5)
  imgdatadf$imgname <- imgname
  imgdatadf$year <- imgyear
  imgdatadf$yearclass <- yearclass
  # extract image month
  imgmonth <- as.numeric(substr(imgname, 13, 14))
  if ((imgmonth <= 3)|(imgmonth >= 10)){
    monthclass <- "dry"
  } else {monthclass <- "wet"}
  imgdatadf$month <- imgmonth
  imgdatadf$monthclass <- monthclass
  print(paste0("Finished: ", imgname))
  return(na.omit(imgdatadf))
}

imgdata_all_list <- lapply(imglist, extractdata)
imgdata_all <- bind_rows(imgdata_all_list)

# Extract terrain (slope and northness)
sloperas <- raster("D:/landsat/dtm/slope_30m.tif")
northras <- raster("D:/landsat/dtm/northness_30m.tif")
df_terrain <- data.frame("PointID" = trainpt_sf$PointID,
                         "Slope" = extract(sloperas, trainpt_sf),
                         "Northness" = extract(northras, trainpt_sf))
imgdata_all_t <- merge(imgdata_all, df_terrain, by="PointID")

write.csv(imgdata_all_t, "D:/landsat/train/imgdata_all_t.csv")


# extract from images without calibration
imgfolder <- "D:/landsat/preprocess/variables"
# repeat the steps above (lines 6-52)
write.csv(imgdata_all_t, "D:/landsat/train/imgdata_all_t_nocal.csv")


# extract from composite
library("glcm")

imgfolder <- "D:/landsat/preprocess/composite"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "median")

trainpt <- read.csv("D:/landsat/train/trainpt.csv")
trainpt_sf <- st_as_sf(trainpt, coords=c("x", "y"), crs = crs(raster(imglist[1])))

focalallband <- function(img){
  blist <- list()
  for (i in 1:nlayers(img)){
    x <- img[[i]]
    x_focal <- focal(x, w=matrix(1,3,3), fun=mean, na.rm=TRUE, NAonly=TRUE)
    blist <- append(blist, x_focal)
  }
  img_focal <- brick(blist)
  return(img_focal)
}

computevar <- function(namedimg){
  gband <- namedimg[["G"]]
  rband <- namedimg[["R"]]
  nirband <- namedimg[["NIR"]]
  ndvi <- (nirband-rband)/(nirband+rband)
  gndvi <- (nirband-gband)/(nirband+gband)
  glcm_mean_con <- glcm(nirband, window = c(5, 5),
                        shift=list(c(0,1), c(1,1), c(1,0), c(1,-1)),
                        statistics = c("mean", "contrast"))
  s <- stack(list(namedimg, ndvi, gndvi, glcm_mean_con))
  names(s) <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con")
  return(s)
}

extractdata <- function(imgname){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("G","R","NIR")
  img <- computevar(img)
  imgdatadf <- extract(img, trainpt_sf, df=TRUE)
  imgdatadf <- cbind(imgdatadf, "PointID" = trainpt_sf$PointID)
  # extract image year class
  startchar <- nchar(imgname)-12
  endchar <- nchar(imgname)-9
  imgyear <- as.numeric(substr(imgname, startchar, endchar))+2
  yearclass <- paste0("c", round(imgyear/5)*5)
  imgdatadf$imgname <- imgname
  imgdatadf$year <- imgyear
  imgdatadf$yearclass <- yearclass
  print(paste0("Finished: ", imgname))
  return(na.omit(imgdatadf))
}

imgdata_all_list <- lapply(imglist, extractdata)
imgdata_all <- bind_rows(imgdata_all_list)


extractsd <- function(imgname){
  img <- brick(imgname)
  names(img) <- c("Gsd","Rsd","NIRsd","ndvisd","gndvisd")
  imgdatadf <- extract(img, trainpt_sf, df=TRUE)
  imgdatadf <- cbind(imgdatadf, "PointID" = trainpt_sf$PointID)
  # extract image year class
  startchar <- nchar(imgname)-12
  endchar <- nchar(imgname)-9
  imgyear <- as.numeric(substr(imgname, startchar, endchar))+2
  yearclass <- paste0("c", round(imgyear/5)*5)
  imgdatadf$imgname <- imgname
  imgdatadf$year <- imgyear
  imgdatadf$yearclass <- yearclass
  print(paste0("Finished: ", imgname))
  return(imgdatadf)
}
imglist <- list.files(imgfolder, pattern = "sd")
imgdata_all_list_sd <- lapply(imglist, extractsd)
imgdata_all_sd <- bind_rows(imgdata_all_list_sd)
imgdata_all_sd[is.na(imgdata_all_sd)] <- 0
imgdata_all_s <- merge(imgdata_all, imgdata_all_sd, by=c("PointID","yearclass"))

# Extract terrain (slope and northness)
sloperas <- raster("D:/landsat/dtm/slope_30m.tif")
northras <- raster("D:/landsat/dtm/northness_30m.tif")
df_terrain <- data.frame("PointID" = trainpt_sf$PointID,
                         "Slope" = extract(sloperas, trainpt_sf),
                         "Northness" = extract(northras, trainpt_sf))
imgdata_all_mediancom <- merge(imgdata_all_s, df_terrain, by="PointID")
head(imgdata_all_mediancom)
write.csv(imgdata_all_mediancom, "D:/landsat/train/imgdata_all_mediancom.csv")


