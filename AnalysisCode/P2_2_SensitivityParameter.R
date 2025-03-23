library(terra)
library(data.table)
library(pcaPP)
library(psych)

# Load habitat maps and threat source shapefile

setwd("D:/HabitatQuality/")
build_all <- vect("threat/Buildings.shp")
highbuild_all <- build_all[build_all$Density=="High"]
lowbuild_all <- build_all[build_all$Density=="Low"]
road_all <- vect("threat/Roads.shp")
mainroad_all <- road_all[road_all$Class=="Main"]
secroad_all <- road_all[road_all$Class=="Secondary"]
pollution_all <- vect("threat/Pollution.shp")
yearlist <- seq(1975,2020,5)
maplist <- list.files("inputmap", "8class.tif$")
threatstable <- read.csv("parameter/threat_table.csv")

# Compute sensitivity

dt <- data.table()
for (i in 2:10){
  targetyear <- yearlist[i]
  habmap <- rast(file.path("inputmap",maplist[i]))       # habitat in this period
  prehabmap <- rast(file.path("inputmap",maplist[i-1]))  # habitat in previous period
  
  highbuild_y <- highbuild_all[highbuild_all$AppearYear<=targetyear]   # threat source in this period
  lowbuild_y <- lowbuild_all[lowbuild_all$AppearYear<=targetyear]
  mainroad_y <- mainroad_all[mainroad_all$AppearYear<=targetyear]
  secroad_y <- secroad_all[secroad_all$AppearYear<=targetyear]
  pollution_y <- pollution_all[pollution_all$AppearYear<=targetyear & pollution_all$EndYear>=targetyear]
  
  # buffer half of maximum impact distance
  highbuild_y_buf <- buffer(highbuild_y, threatstable[threatstable$THREAT=="highbuild","MAX_DIST"]*1000/2)
  highbuild_y_bufras <- mask(rasterize(highbuild_y_buf, habmap, background=0), habmap)
  lowbuild_y_buf <- buffer(lowbuild_y, threatstable[threatstable$THREAT=="lowbuild","MAX_DIST"]*1000/2)
  lowbuild_y_bufras <- mask(rasterize(lowbuild_y_buf, habmap, background=0), habmap)
  mainroad_y_buf <- buffer(mainroad_y, threatstable[threatstable$THREAT=="mainroad","MAX_DIST"]*1000/2)
  mainroad_y_bufras <- mask(rasterize(mainroad_y_buf, habmap, background=0), habmap)
  secroad_y_buf <- buffer(secroad_y, threatstable[threatstable$THREAT=="secroad","MAX_DIST"]*1000/2)
  secroad_y_bufras <- mask(rasterize(secroad_y_buf, habmap, background=0), habmap)
  pollution_y_buf <- buffer(pollution_y, threatstable[threatstable$THREAT=="pollution","MAX_DIST"]*1000/2)
  pollution_y_bufras <- mask(rasterize(pollution_y_buf, habmap, background=0), habmap)
  
  # create a data frame of habitat changes (this & previous period) and occurrence of threats
  ras <- c(habmap, prehabmap, highbuild_y_bufras, lowbuild_y_bufras,
           mainroad_y_bufras, secroad_y_bufras, pollution_y_bufras)
  names(ras) <- c("curhabitat","prehabitat","highbuild","lowbuild",
                  "mainroad","secroad","pollution")
  dt_y <- as.data.table(as.data.frame(ras))   
  dt_y <- na.omit(dt_y)
  dt <- rbind(dt, dt_y)
}

# Phi coefficients
df_sen <- data.frame()
for (h in c(1,2,3,4,5,6,7,8)) {  # for each habitat class
  df1 <- data.frame(habitat=h)
  dt_h <- dt[dt$curhabitat==h | dt$prehabitat==h,]                            # subset class in either period == h
  dt_h <- dt_h[, change := ifelse(curhabitat==prehabitat, 0, 1)]              # both period ==h (no change) -> 0, either one (have change) -> 1
  for (t in c("highbuild","lowbuild","mainroad","secroad","pollution")) {
    df1 <- cbind(df1, phic = phi(table(dt_h$change, dt_h[[t]]), digits = 5))  # phi, actually same result as kendall
    colnames(df1)[ncol(df1)] <- t
  }
  df_sen <- rbind(df_sen, df1)
}
write.csv(df_sen, "parameter/sensitivity.csv")


