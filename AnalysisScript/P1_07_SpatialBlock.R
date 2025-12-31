library(sf)
library(terra)
library(raster)
library(dplyr)
library(glcm)
library(blockCV)
library(data.table)

imgfolder <- "D:/landsat/preprocess/composite"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "median")
rasters <- lapply(imglist, rast)

# create a mean raster from all years
nR <- length(rasters)
nBands <- nlyr(rasters[[1]])
mean_layers <- vector("list", nBands)
for (b in 1:nBands) {
  # extract band b from every raster and stack them
  band_stack <- rast(lapply(rasters, function(x) x[[b]]))
  # compute mean across the stack (pixel-wise)
  mean_layers[[b]] <- app(band_stack, mean, na.rm=TRUE)
}
mean_by_band <- rast(mean_layers)

# compute predictor variables from the mean raster
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

img <- brick(mean_by_band)
names(img) <- c("G","R","NIR")
img <- computevar(img)
img <- rast(img)
crs(img)  <- "epsg:2326"
sloperas <- rast("D:/landsat/dtm/slope_30m.tif")
northras <- rast("D:/landsat/dtm/northness_30m.tif")
img <- c(img,sloperas,northras)

aoi <- vect("D:/landsat/ReferenceData/Land.shp")
img <- crop(img, aoi, mask=TRUE)

# measure spatial autocorrelation from predictor raster (for choosing block size)
set.seed(2022)
sac1 <- cv_spatial_autocor(r = img)

# assign reference data to blocks (10 fold)
setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_t.csv")
trainpt <- read.csv("trainpt.csv")
trainpt_m <- reshape2::melt(trainpt[,c(1,3:14)], id.vars = c("PointID","x","y"),
                            variable.name = "yearclass", value.name = "Class")
pa_data <- sf::st_as_sf(trainpt_m, coords = c("x", "y"), crs = 2326)

set.seed(2022)
sb1 <- cv_spatial(pa_data, column="Class", r=img, k=10, size=sac1$range, 
                  selection="random", iteration=100, seed=2022)
saveRDS(sb1,"spatialblock.rds")

