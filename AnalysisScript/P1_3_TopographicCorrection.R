library(raster)
library(insol)
dtm <- raster("D:/landsat/dtm/dtm_30m.tif")
dtm_slope <- raster("D:/landsat/dtm/slope_30m.tif")
dtm_aspect <- raster("D:/landsat/dtm/aspect_30m.tif")  # computed from ArcGIS

# SCS+C topographic correction
# https://ieeexplore.ieee.org/document/1499030
# modified based on topocorr function in landsat package

scsc <- function(ras, sunazimuth, sunelev){
  # compute terrain shaded area
  ras_shade <- insol::doshade(dtm, normalvector(90-sunelev, sunazimuth)) # (sun zenith, sun azimuth)
  # compute illumination
  slope_rad <- (pi/180) * dtm_slope
  aspect_rad <- (pi/180) * dtm_aspect
  sunzenith_rad <- (pi/180) * (90 - sunelev)
  sunazimuth_rad <- (pi/180) * sunazimuth
  IL <- cos(slope_rad)*cos(sunzenith_rad)+sin(slope_rad)*sin(sunzenith_rad)*cos(sunazimuth_rad - aspect_rad)
  IL[IL < 0] <- NA
  # mask shaded area (remove from data)
  ras <- mask(ras, ras_shade, maskvalue = 0)
  ras <- mask(ras, IL)
  # mask flat area (no need correction)
  flatarea <- dtm_slope < 5
  ras_mask <- mask(ras, flatarea, maskvalue = 1)
  # loop through bands
  blist <- list()
  ras_na <- ras[[1]]*0
  for (i in 1:nlayers(ras)){
    x_raw <- ras[[i]]
    x <- ras_mask[[i]]
    # c correction
    band.lm <- lm(as.vector(x) ~ as.vector(IL))
    C <- coefficients(band.lm)[[1]]/coefficients(band.lm)[[2]]
    # SCS: xout <- x * (cos(sunzenith) * cos(slope))/IL
    xout <- x * (cos(sunzenith_rad) * cos(slope_rad) + C)/(IL + C)
    # merge corrected and original areas
    x_cor <- cover(xout, x_raw)
    # clamp extreme values
    x_cor[(x_cor < -0.05) | (x_cor > 1)] <- NA
    x_cor[x_cor < 0] <- 0
    blist <- append(blist, x_cor)
    # apply masked pixels to all bands
    ras_na <- mask(ras_na, x_cor)
  }
  ras_cor <- brick(blist)
  ras_cor <- mask(ras_cor, ras_na)
  return(ras_cor)
}

imgfolder <- "D:/landsat/preprocess/rectangle"
outfolder <- "D:/landsat/preprocess/topocor"
angle_csv <- read.csv("D:/landsat/preprocess/dtm/df_metadata.csv")

setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "*.tif$")
n <- 1
for (f in imglist){
  print(paste0("Processing: ", n, " / ", length(imglist)))
  outfile <- file.path(outfolder,f)
  if (file.exists(outfile) == TRUE){
    n <- n+1
    next
  }
  ras <- brick(f)
  sunazimuth <- angle_csv[angle_csv$filename == f,]$SunAzimuth
  sunelev <- angle_csv[angle_csv$filename == f,]$SunElevation
  ras_cor <- scsc(ras, sunazimuth, sunelev)
  writeRaster(ras_cor, outfile)
  print(paste0("Exported: ", outfile))
  n <- n+1
}


