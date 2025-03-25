library(terra)
library(lubridate)

imgfolder <- "D:/landsat/preprocess/CrossSensorCal"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "tif$")

filteryear <- function(imglist, startyear, endyear){
  imglist_year <- lapply(imglist, substr, 9, 12)
  yearlist <- startyear:endyear
  for (i in 1:length(yearlist)){
    containyear1 <- imglist_year == yearlist[i]
    if (i==1){containyear <- containyear1}
    containyear <- containyear | containyear1
  }
  return(imglist[containyear])
}

KeepCommonBand <- function(imgname){
  img <- rast(imgname)
  sensorname <- substr(imgname, 20, 21)
  if (sensorname == "MS"){
    names(img) <- c("G","R","RE","NIR")
  } else if (sensorname == "TM" | sensorname == "ET"){
    names(img) <- c("B","G","R","NIR","SW1","SW2")
  } else if (sensorname == "OL"){
    names(img) <- c("C","B","G","R","NIR","SW1","SW2")
  }
  return(img[[c("G","R","NIR")]])
}

# median image

startyear1 <- c(1973,1978,1983,1988,1993,1998,2003,2008,2013,2018)
endyear1 <- c(1977,1982,1987,1992,1997,2002,2007,2012,2017,2022)
for (i in 1:10){
  startyear <- startyear1[i]
  endyear <- endyear1[i]
  imglist1 <- filteryear(imglist, startyear, endyear)
  imglist1 <- lapply(imglist1, KeepCommonBand)
  rsrc <- sprc(imglist1)
  outname <- paste0("D:/landsat/preprocess/composite/median_",startyear,"-",endyear,".tif")
  mosaic(rsrc, fun="median", filename=outname)
}

# Standard deviation image

computevi <- function(namedimg){
  gband <- namedimg[["G"]]
  rband <- namedimg[["R"]]
  nirband <- namedimg[["NIR"]]
  ndvi <- (nirband-rband)/(nirband+rband)
  names(ndvi) <- "ndvi"
  gndvi <- (nirband-gband)/(nirband+gband)
  names(gndvi) <- "gndvi"
  return(c(namedimg, ndvi, gndvi))
}

startyear1 <- c(1973,1978,1983,1988,1993,1998,2003,2008,2013,2018)
endyear1 <- c(1977,1982,1987,1992,1997,2002,2007,2012,2017,2022)
for (i in 1:10){
  startyear <- startyear1[i]
  endyear <- endyear1[i]
  imglist1 <- filteryear(imglist, startyear, endyear)
  imglist1 <- lapply(imglist1, KeepCommonBand)
  imglist1 <- lapply(imglist1, computevi)
  rastset <- sds(imglist1)
  outname <- paste0("D:/landsat/preprocess/composite/sd_",startyear,"-",endyear,".tif")
  app(rastset, sd, na.rm=TRUE, filename=outname)
}

