#### Run lines 1-255 in "Landsat_8_TrainRFModel" first to load libraries and required functions

# Fusion method (pixel-level composite)

setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_mediancom.csv")
trainpt <- read.csv("trainpt.csv")

trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                            variable.name = "yearclass", value.name = "Class")
trainpt_m$Class <- factor(trainpt_m$Class, levels = c(1:6))
levels(trainpt_m$Class) <- c("c1", "c2", "c3", "c4", "c5", "c6")
traindata <- merge(imgdf, trainpt_m, by = c("PointID","yearclass"))
traindata$row_num <- seq.int(nrow(traindata))

predictors_withsd <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con",
                       "Gsd","Rsd","NIRsd","ndvisd","gndvisd","Slope","Northness")
predictors_nosd <- c("G","R","NIR","ndvi","gndvi",
                     "glcm_mean","glcm_con","Slope","Northness")

cv_pred <- buildmodel(traindata, predictors_withsd, separatewetdry = FALSE)
saveRDS(cv_pred, "D:/landsat/train/cv_pred_mediancom_withsd.rds")
cv_pred <- buildmodel(traindata, predictors_nosd, separatewetdry = FALSE)
saveRDS(cv_pred, "D:/landsat/train/cv_pred_mediancom_nosd.rds")

pred_df <- bind_rows(cv_pred[c(3,7)])
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))

pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)


# Period-specific model

setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_t.csv")
trainpt <- read.csv("trainpt.csv")

trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                            variable.name = "yearclass", value.name = "Class")
trainpt_m$Class <- factor(trainpt_m$Class, levels = c(1:6))
levels(trainpt_m$Class) <- c("c1", "c2", "c3", "c4", "c5", "c6")
traindata <- merge(imgdf, trainpt_m, by = c("PointID","yearclass"))
traindata$row_num <- seq.int(nrow(traindata))

predictors <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")

yc_list <- unique(traindata$yearclass)
for (yc in yc_list){
  td1 <- traindata[traindata$yearclass==yc,]
  cv_pred <- buildmodel(td1, predictors)
  saveRDS(cv_pred, paste0("D:/landsat/train/cv_pred_sepperiod_",yc,".rds"))
}
sepyc <- list.files("D:/landsat/train/", pattern="cv_pred_sepperiod", full.names=TRUE)
getdf <- function(rdsname){
  f <- readRDS(rdsname)
  df <- bind_rows(f[c(3,7,11,15)])
  return(df)
}
sepyc <- lapply(sepyc, getdf)
pred_df <- bind_rows(sepyc)
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))
pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult_sepperiod.csv")


# Sensor-specific model

traindata$sensor <- substr(traindata$imgname, 20, 21)
sen_list <- unique(traindata$sensor)
for (sen in sen_list){
  td1 <- traindata[traindata$sensor==sen,]
  addvar <- c("ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")
  if (sen == "MS"){ # predictors add remaining bands
    predictors <- append(c("G","R","RE","NIR"), addvar)
  } else if (sen == "TM" | sen == "ET"){
    predictors <- append(c("B","G","R","NIR","SW1","SW2"), addvar)
  } else if (sen == "OL"){
    predictors <- append(c("C","B","G","R","NIR","SW1","SW2"), addvar)
  }
  cv_pred <- buildmodel(td1, predictors)
  saveRDS(cv_pred, paste0("D:/landsat/train/cv_pred_sepsensor_",sen,".rds"))
}
sepsen <- list.files("D:/landsat/train/", pattern="cv_pred_sepsensor", full.names=TRUE)
getdf <- function(rdsname){
  f <- readRDS(rdsname)
  df <- bind_rows(f[c(3,7,11,15)])
  return(df)
}
sepsen <- lapply(sepsen, getdf)
pred_df <- bind_rows(sepsen)
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))
pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult_sepsensor.csv")


# Image-specific model

img_list <- unique(traindata$imgname)
for (img in img_list){
  td1 <- traindata[traindata$imgname==img,]
  sen <- substr(img, 20, 21)
  addvar <- c("ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")
  if (sen == "MS"){ # predictors add remaining bands
    predictors <- append(c("G","R","RE","NIR"), addvar)
  } else if (sen == "TM" | sen == "ET"){
    predictors <- append(c("B","G","R","NIR","SW1","SW2"), addvar)
  } else if (sen == "OL"){
    predictors <- append(c("C","B","G","R","NIR","SW1","SW2"), addvar)
  }
  cv_pred <- buildmodel(td1, predictors, separatewetdry = FALSE)
  saveRDS(cv_pred, paste0("D:/landsat/train/cv_pred_sepimage_",substr(img, 1, 28),".rds"))
}
sepsen <- list.files("D:/landsat/train/", pattern="cv_pred_sepimage", full.names=TRUE)
getdf <- function(rdsname){
  f <- readRDS(rdsname)
  df <- bind_rows(f[c(3,7)])
  return(df)
}
sepsen <- lapply(sepsen, getdf)
pred_df <- bind_rows(sepsen)
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))
pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult_sepimage.csv")


# Without cross-sensor calibration

setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_t_nocal.csv")
trainpt <- read.csv("trainpt.csv")

trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                            variable.name = "yearclass", value.name = "Class")
trainpt_m$Class <- factor(trainpt_m$Class, levels = c(1:6))
levels(trainpt_m$Class) <- c("c1", "c2", "c3", "c4", "c5", "c6")
traindata <- merge(imgdf, trainpt_m, by = c("PointID","yearclass"))
traindata$row_num <- seq.int(nrow(traindata))

predictors <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")

cv_pred <- buildmodel(traindata, predictors)
saveRDS(cv_pred, "D:/landsat/train/cv_pred_nocal.rds")

pred_df <- bind_rows(cv_pred[c(3,7,11,15)])
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))

pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult_nocal.csv")


# No temporal smoothing

cv_pred <- readRDS("D:/landsat/train/cv_pred.rds")
pred_df <- bind_rows(cv_pred[c(3,7,11,15)])
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))

pred_df_vote <- voteresult(pred_df, temporalsmooth = FALSE)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult_notempsmooth.csv")


# Reduce number of images

library(terra)
library(sf)
library(dplyr)

imgdir <- "D:/landsat/preprocess/variables"
setwd(imgdir)
imglist <- list.files(imgdir, pattern="*.tif$")
landaoi <- st_read("D:/landsat/ReferenceData/Land.shp")

df <- data.frame(imgname=imglist)
df$yearclass <- NA
df$pathrow <- NA
df$validpixel <- 0
for (imgname in imglist){
  ras <- rast(imgname)[[1]]
  rascrop <- crop(ras, landaoi, mask=TRUE)
  rasfreq <- freq(!is.na(rascrop))
  validcount <- rasfreq[rasfreq$value==1,]$count
  imgyear <- as.numeric(substr(imgname, 9, 12))
  yearclass <- paste0("c", round(imgyear/5)*5)
  df$yearclass[df$imgname==imgname] <- yearclass
  df$pathrow[df$imgname==imgname] <- substr(imgname,23,28)
  df$validpixel[df$imgname==imgname] <- validcount
}

df$path <- substr(df$pathrow,1,3)
percent_df <- df %>%
  group_by(yearclass, path) %>%
  mutate(percentile=percent_rank(validpixel)*100)
write.csv(percent_df, "D:/landsat/train/img_validcount.csv")

for (p in c(seq(90,10,-10),5,1)){
  percent_df_subset <- percent_df[percent_df$percentile >= (100-p),]
  percent_df_subset <- percent_df_subset[,c("imgname","percentile")]
  td1 <- merge(traindata, percent_df_subset, by="imgname", all=FALSE)
  if (sum(td1$monthclass=="wet")>0 & sum(td1$monthclass=="dry")>0){
    cv_pred <- buildmodel(td1, predictors, separatewetdry = TRUE)
  } else {
    cv_pred <- buildmodel(td1, predictors, separatewetdry = FALSE)
  }
  pred_df <- bind_rows(cv_pred[c(3,7,11,15)]) 
  pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))
  pred_df_vote <- voteresult(pred_df)
  vresult <- doaccuracy(pred_df_vote)
  vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
  write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),
            paste0("vresult_reduceto",p,".csv"))
}

