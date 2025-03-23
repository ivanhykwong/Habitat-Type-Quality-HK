library(raster)
library(sf)
library(dplyr)

calpt <- st_read("D:/landsat/referencedata/UnchangePoint.shp")
calpt_df <- st_drop_geometry(calpt)
calpt_df$ID <- 1:nrow(calpt_df)

# MSS

imgfolder <- "D:/landsat/preprocess/topocor"
setwd(imgfolder)
mss <- list.files(imgfolder, pattern = "MS.+.tif$")
mss <- lapply(mss, brick)

mss_df <- lapply(mss, mss_extract)
mss_df <- bind_rows(mss_df)
mss_df <- mss_df %>% group_by(ID) %>%
  summarise(G_mssm = median(G), R_mssm = median(R), 
            NIR_B3_mssm = median(NIR_B3), NIR_mssm = median(NIR),
            Class = median(Class))
mss_df$Sensor <- "mss"
write.csv(mss_df, "D:/landsat/CrossSensorCal/mss_df.csv")

# TM

imgfolder <- "D:/landsat/preprocess/topocor"
setwd(imgfolder)
tm <- list.files(imgfolder, pattern = "TM.+.tif$")
tm <- lapply(tm, brick)

tm_extract <- function(img){
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_df <- extract(img, calpt, df = TRUE)
  img_df <- merge(img_df, calpt_df, by = "ID")
  img_df <- na.omit(img_df)
  return(img_df)
}

tm_df <- lapply(tm, tm_extract)
tm_df <- bind_rows(tm_df)
tm_df <- tm_df %>% group_by(ID) %>%
  summarise(B_tm = median(B), G_tm = median(G), R_tm = median(R), 
            NIR_tm = median(NIR), SW1_tm = median(SW1), SW2_tm = median(SW2),
            Class = median(Class))
tm_df$Sensor <- "TM"
write.csv(tm_df, "D:/landsat/CrossSensorCal/tm_df.csv")

# ETM

imgfolder <- "D:/landsat/preprocess/topocor"
setwd(imgfolder)
et <- list.files(imgfolder, pattern = "ET.+.tif$")
et <- lapply(et, brick)

et_extract <- function(img){
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_df <- extract(img, calpt, df = TRUE)
  img_df <- merge(img_df, calpt_df, by = "ID")
  img_df <- na.omit(img_df)
  return(img_df)
}

et_df <- lapply(et, et_extract)
et_df <- bind_rows(et_df)
et_df <- et_df %>% group_by(ID) %>%
  summarise(B_et = median(B), G_et = median(G), R_et = median(R), 
            NIR_et = median(NIR), SW1_et = median(SW1), SW2_et = median(SW2),
            Class = median(Class))
et_df$Sensor <- "et"
write.csv(et_df, "D:/landsat/CrossSensorCal/et_df.csv")

# OLI

imgfolder <- "D:/landsat/preprocess/topocor"
setwd(imgfolder)
ol <- list.files(imgfolder, pattern = "OL.+.tif$")
ol <- lapply(ol, brick)

ol_extract <- function(img){
  names(img) <- c("C", "B", "G", "R", "NIR", "SW1", "SW2")
  img_df <- extract(img, calpt, df = TRUE)
  img_df <- merge(img_df, calpt_df, by = "ID")
  img_df <- na.omit(img_df)
  return(img_df)
}

ol_df <- lapply(ol, ol_extract)
ol_df <- bind_rows(ol_df)
ol_df <- ol_df %>% group_by(ID) %>%
  summarise(C_ol = median(C), B_ol = median(B), G_ol = median(G), R_ol = median(R), 
            NIR_ol = median(NIR), SW1_ol = median(SW1), SW2_ol = median(SW2),
            Class = median(Class))
ol_df$Sensor <- "ol"
write.csv(ol_df, "D:/landsat/CrossSensorCal/ol_df.csv")

# Merge

setwd("D:/landsat/CrossSensorCal")
mss_df <- read.csv("mss_df.csv")[-1]
tm_df <- read.csv("tm_df.csv")[-1]
et_df <- read.csv("et_df.csv")[-1]
ol_df <- read.csv("ol_df.csv")[-1]

df <- merge(mss_df, tm_df[1:(ncol(tm_df)-2)], by = "ID")
df <- merge(df, et_df[1:(ncol(et_df)-2)], by = "ID")
df <- merge(df, ol_df[1:(ncol(ol_df)-2)], by = "ID")

write.csv(df, "CrossSensorCal.csv")

# Calculate LM for each band, return equations

lm_cal <- function(x, y){
  a <- lm(df[[y]]~df[[x]])
  rsquare <- summary(a)$r.squared
  intercept <- coefficients(a)[[1]]
  x_coef <- coefficients(a)[[2]]
  x_cal <- df[[x]]*x_coef+intercept
  x_cal_mae <- mean(abs(x_cal-df[[y]]))
  lm_df <- data.frame(x=x, y=y, rsquare=rsquare, mae=x_cal_mae,
                      intercept=intercept, x_coef=x_coef)
  return(lm_df)
}

lm_result_list <- list(lm_cal("G_mssm", "G_ol"), lm_cal("R_mssm", "R_ol"), lm_cal("NIR_mssm", "NIR_ol"),
                       lm_cal("B_tm", "B_ol"), lm_cal("G_tm", "G_ol"), lm_cal("R_tm", "R_ol"),
                       lm_cal("NIR_tm", "NIR_ol"), lm_cal("SW1_tm", "SW1_ol"), lm_cal("SW2_tm", "SW2_ol"),
                       lm_cal("B_et", "B_ol"), lm_cal("G_et", "G_ol"), lm_cal("R_et", "R_ol"),
                       lm_cal("NIR_et", "NIR_ol"), lm_cal("SW1_et", "SW1_ol"), lm_cal("SW2_et", "SW2_ol"))
lm_result <- bind_rows(lm_result_list)
write.csv(lm_result, "lm_result.csv")


# Apply equations to images

setwd("D:/landsat/CrossSensorCal")
lm_result <- read.csv("lm_result.csv")

CalBand <- function(imgbrick, bandindex, x_name){
  ras <- imgbrick[[bandindex]]
  ras_cal <- ras*(lm_result[lm_result$x==x_name, "x_coef"])+
    (lm_result[lm_result$x==x_name, "intercept"])
  ras_cal[ras_cal > 1] <- 1
  ras_cal[ras_cal < 0] <- 0
  return(ras_cal)
}

mss_cal_fun <- function(imgname){
  imgbrick <- brick(imgname)
  b1_cal <- CalBand(imgbrick, 1, "G_mssm")
  b2_cal <- CalBand(imgbrick, 2, "R_mssm")
  b3_cal <- imgbrick[[3]]
  b4_cal <- CalBand(imgbrick, 4, "NIR_mssm")
  brick_cal <- brick(b1_cal, b2_cal, b3_cal, b4_cal)
  imgname_base <- gsub("_hk1980.tif", "", basename(imgname))
  outfile <- paste0("D:/landsat/preprocess/CrossSensorCal/",imgname_base,"_cal.tif")
  writeRaster(brick_cal, outfile)
}

imgfolder <- "D:/landsat/preprocess/topocor"
mss <- list.files(imgfolder, pattern = "MS.+.tif$", full.names = TRUE)
lapply(mss, mss_cal_fun)

tm_cal_fun <- function(imgname){
  imgbrick <- brick(imgname)
  b1_cal <- CalBand(imgbrick, 1, "B_tm")
  b2_cal <- CalBand(imgbrick, 2, "G_tm")
  b3_cal <- CalBand(imgbrick, 3, "R_tm")
  b4_cal <- CalBand(imgbrick, 4, "NIR_tm")
  b5_cal <- CalBand(imgbrick, 5, "SW1_tm")
  b6_cal <- CalBand(imgbrick, 6, "SW2_tm")
  brick_cal <- brick(b1_cal, b2_cal, b3_cal, b4_cal, b5_cal, b6_cal)
  imgname_base <- gsub("_hk1980.tif", "", basename(imgname))
  outfile <- paste0("D:/landsat/preprocess/CrossSensorCal/",imgname_base,"_cal.tif")
  writeRaster(brick_cal, outfile)
  print(paste0("Finish: ", imgname_base))
}

imgfolder <- "D:/landsat/preprocess/topocor"
tm <- list.files(imgfolder, pattern = "TM.+.tif$", full.names = TRUE)
lapply(tm, tm_cal_fun)

et_cal_fun <- function(imgname){
  imgbrick <- brick(imgname)
  b1_cal <- CalBand(imgbrick, 1, "B_et")
  b2_cal <- CalBand(imgbrick, 2, "G_et")
  b3_cal <- CalBand(imgbrick, 3, "R_et")
  b4_cal <- CalBand(imgbrick, 4, "NIR_et")
  b5_cal <- CalBand(imgbrick, 5, "SW1_et")
  b6_cal <- CalBand(imgbrick, 6, "SW2_et")
  brick_cal <- brick(b1_cal, b2_cal, b3_cal, b4_cal, b5_cal, b6_cal)
  imgname_base <- gsub("_hk1980.tif", "", basename(imgname))
  outfile <- paste0("D:/landsat/preprocess/CrossSensorCal/",imgname_base,"_cal.tif")
  writeRaster(brick_cal, outfile)
  print(paste0("Finish: ", imgname_base))
}

imgfolder <- "D:/landsat/preprocess/topocor"
et <- list.files(imgfolder, pattern = "ET.+.tif$", full.names = TRUE)
lapply(et, et_cal_fun)


imgfolder <- "D:/landsat/preprocess/topocor"
ol <- list.files(imgfolder, pattern = "OL.+.tif$", full.names = TRUE)
for (imgname in ol){
  imgbrick <- brick(imgname)
  imgname_base <- gsub("_hk1980.tif", "", basename(imgname))
  outfile <- paste0("D:/landsat/preprocess/CrossSensorCal/",imgname_base,"_cal.tif")
  writeRaster(imgbrick, outfile)
  print(paste0("Finish: ", imgname_base))
}


# Compute Variables

library("glcm")

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
  return(stack(list(namedimg, ndvi, gndvi, glcm_mean_con)))
}

imgfolder <- "D:/landsat/preprocess/CrossSensorCal"
outfolder <- "D:/landsat/preprocess/variables"
setwd(imgfolder)
ms <- list.files(imgfolder, pattern = "MS.+.tif$")
for (imgname in ms){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("G", "R", "RE", "NIR")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_cal.tif","_calvar.tif",imgname)))
  print(paste0("Finish: ", imgname))
}

tm <- list.files(imgfolder, pattern = "TM.+.tif$")
for (imgname in tm){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_cal.tif","_calvar.tif",imgname)))
  print(paste0("Finish: ", imgname))
}

et <- list.files(imgfolder, pattern = "ET.+.tif$")
for (imgname in et){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_cal.tif","_calvar.tif",imgname)))
  print(paste0("Finish: ", imgname))
}


imgfolder <- "D:/landsat/preprocess/topocor"
outfolder <- "D:/landsat/preprocess/variables"
setwd(imgfolder)
ol <- list.files(imgfolder, pattern = "OL.+.tif$")
for (imgname in ol){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("C", "B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_cal.tif","_calvar.tif",imgname)))
  print(paste0("Finish: ", imgname))
}


# Compute Variables without calibration

library("glcm")

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
  return(stack(list(namedimg, ndvi, gndvi, glcm_mean_con)))
}

imgfolder <- "D:/landsat/preprocess/topocor"
outfolder <- "D:/landsat/preprocess/variables_nocal"
setwd(imgfolder)
ms <- list.files(imgfolder, pattern = "MS.+.tif$")
for (imgname in ms){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("G", "R", "RE", "NIR")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_hk1980.tif","_var.tif",imgname)))
  print(paste0("Finish: ", imgname))
}

tm <- list.files(imgfolder, pattern = "TM.+.tif$")
for (imgname in tm){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_hk1980.tif","_var.tif",imgname)))
  print(paste0("Finish: ", imgname))
}

et <- list.files(imgfolder, pattern = "ET.+.tif$")
for (imgname in et){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_hk1980.tif","_var.tif",imgname)))
  print(paste0("Finish: ", imgname))
}

imgfolder <- "D:/landsat/preprocess/topocor"
outfolder <- "D:/landsat/preprocess/variables_nocal"
setwd(imgfolder)
ol <- list.files(imgfolder, pattern = "OL.+.tif$")
for (imgname in ol){
  img <- brick(imgname)
  img <- focalallband(img)
  names(img) <- c("C", "B", "G", "R", "NIR", "SW1", "SW2")
  img_withvar <- computevar(img)
  writeRaster(img_withvar, file.path(outfolder,gsub("_hk1980.tif","_var.tif",imgname)))
  print(paste0("Finish: ", imgname))
}
