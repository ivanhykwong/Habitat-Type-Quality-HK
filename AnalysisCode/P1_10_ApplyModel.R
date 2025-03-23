library(raster)
library(terra)
library(sf)
library(dplyr)
library(data.table)
library(reshape2)
library(caret)
library(randomForest)
library(foreach)
library(doParallel)
library(kit)
library(groupdata2)
library(e1071)
library(ranger)
library(matrixStats)

imgfolder <- "D:/landsat/preproces/variables"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "tif$")

sloperas <- raster("D:/landsat/dtm/slope_30m.tif")
northras <- raster("D:/landsat/dtm/north_30m.tif")
predictors <- c("G", "R", "NIR", "ndvi","gndvi","glcm_mean","glcm_con", "Slope", "Northness")

predictlandsat <- function(imgname, predictors, rf_model){
  img <- stack(list(imgname, sloperas, northras))
  # check sensor and assign bands
  sensorname <- substr(imgname, 20, 21)
  addvar <- c("ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")
  if (sensorname == "MS"){
    names(img) <- append(c("G","R","RE","NIR"), addvar)
  } else if (sensorname == "TM" | sensorname == "ET"){
    names(img) <- append(c("B","G","R","NIR","SW1","SW2"), addvar)
  } else if (sensorname == "OL"){
    names(img) <- append(c("C","B","G","R","NIR","SW1","SW2"), addvar)
  }
  
  # extract image month
  imgmonth <- as.numeric(substr(imgname, 13, 14))
  if ((imgmonth <= 3)|(imgmonth >= 10)){
    monthclass <- "dry"
  } else {monthclass <- "wet"}

  # predictors
  p_yes <- predictors
  if("Slope" %in% p_yes == F){p_yes <- append(p_yes,"Slope")}
  if("Northness" %in% p_yes == F){p_yes <- append(p_yes,"Northness")}
  p_no <- predictors
  if("Slope" %in% p_no == T){p_no <- p_no[!p_no == "Slope"]}
  if("Northness" %in% p_no == T){p_no <- p_no[!p_no == "Northness"]}
  
  # choose model
  if (monthclass == "wet"){
    modeluse_t <- rf_model$wet_withtopo$finalModel
    modeluse_nt <- rf_model$wet_notopo$finalModel
  } else {
    modeluse_t <- rf_model$dry_withtopo$finalModel
    modeluse_nt <- rf_model$dry_notopo$finalModel
  }
  
  # apply model
  if (sensorname == "OL"){
    outname <- gsub("_hk1980.tif","_prob_t.tif",imgname)             # check L8 name
  } else {
    outname <- gsub("_calvar.tif","_prob_t.tif",imgname)
  }
  raster::predict(img[[p_yes]], model = modeluse_t, type="prob",     # with topo
                 filename = file.path(outfolder,outname), 
                 index=1:6, na.rm=TRUE)                               
  if (sensorname == "OL"){
    outname <- gsub("_hk1980.tif","_prob_nt.tif",imgname)            # check L8 name
  } else {
    outname <- gsub("_calvar.tif","_prob_nt.tif",imgname)
  }
  raster::predict(img[[p_no]], model = modeluse_nt, type="prob",     # without topo 
                 filename = file.path(outfolder,outname),
                 index=1:6, na.rm=TRUE)                              
  print(paste0("Finish: ", imgname))
}

outfolder <- "D:/landsat/predict"
rf_model <- readRDS("D:/landsat/train/rf_model.rds")
# parallel
cl <- parallel::makeCluster(3)
doParallel::registerDoParallel(cl)
foreach (i=1:length(imglist), .packages=c("raster","randomForest")) %dopar% {
  imgname <- imglist[i]
  predictlandsat(imgname, predictors, rf_model)
}
parallel::stopCluster(cl)


# Sum images within same yearclass

imgfolder <- "D:/landsat/predict"
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

extractband <- function(img, bnum){
  return(img[[bnum]])
}

calcsum <- function(imglist, startyear, endyear){
  imglist1 <- filteryear(imglist, startyear, endyear)
  raslist <- lapply(imglist1, rast)
  lyrsum <- function(bnum){
    lyr1 <- lapply(raslist, extractband, bnum = bnum)
    return(sum(rast(lyr1),na.rm=TRUE))
  }
  sumlist <- list()
  for (i in 1:6){
    s <- lyrsum(i)
    sumlist <- append(sumlist, list(s))
  }
  sumras <- rast(sumlist)
  return(sumras)
}

fillna_prob <- function(ras){
  while (global(is.na(ras),sum)[[1]] > 0){
    ras <- focal(ras, w=3, fun="mean", na.rm=TRUE, na.policy="only")
  }
  return(ras)
}
fillna_prob_allband <- function(ras){
  for (i in 1:nlyr(ras)){
    ras[[i]] <- fillna_prob(ras[[i]])
  }
  return(ras)
}

startyear1 <- c(1973,1978,1983,1988,1993,1998,2003,2008,2013,2018)
endyear1 <- c(1977,1982,1987,1992,1997,2002,2007,2012,2017,2022)
probmaplist <- list()
for (i in 1:10){
  startyear <- startyear1[i]
  endyear <- endyear1[i]
  outputprobmap <- calcsum(imglist, startyear, endyear)
  outputprobmap <- fillna_prob_allband(outputprobmap)
  probmaplist[[i]] <- outputprobmap
}

# Temporal smoothing and get highest prob class
relprob <- function(ras){
  totalprob <- sum(ras)
  return(ras/totalprob)
}
setwd("D:/landsat/outputmap")
probmaplist <- lapply(probmaplist, relprob)
for (i in 1:10){
  startyear <- startyear1[i]
  endyear <- endyear1[i]
  if (i==1){
    smoothprobmap <- (2*probmaplist[[i]] + probmaplist[[i+1]])/3
  } else if (i==10) {
    smoothprobmap <- (2*probmaplist[[i]] + probmaplist[[i-1]])/3
  } else {
    smoothprobmap <- (2*probmaplist[[i]] + probmaplist[[i-1]] + probmaplist[[i+1]])/4
  }
  writeRaster(smoothprobmap, paste0("prob/probmap_",startyear,"-",endyear,".tif"))
  outputmap <- which.max(smoothprobmap)
  writeRaster(outputmap, paste0("outmap_",startyear,"-",endyear,".tif"))
}



