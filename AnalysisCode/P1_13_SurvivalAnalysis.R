library(terra)
library(survival)
library(ggsurvfit)
library(ggplot2)

##### Create dense habitat class raster for all 45 years

# read raw prob tif
rasdir <- "D:/landsat/outputmap/prob"
setwd(rasdir)
probmapname <- list.files(rasdir, pattern="tif$")
probmaplist <- lapply(probmapname, rast)

# interpolate and create dense prob tif
setwd("D:/landsat/outputmap/dense")
for (i in 1:46){
  premap <- probmaplist[[floor(i/5+0.8)]]
  nexmap  <- probmaplist[[ceiling(i/5+0.8)]]
  preyearlist <- c(rep(c(0,1,2,3,4),9),0)
  nexyearlist <- c(rep(c(0,4,3,2,1),9),0)
  preyear <- preyearlist[i]
  nexyear <- nexyearlist[i]
  if (preyear==0) {outmap <- premap}
  else {outmap <- (premap*nexyear+nexmap*preyear)/5}
  outmap <- which.max(outmap)
  writeRaster(outmap, paste0("outmap",i,".tif"))
}


##### Survival Analysis and Median Number of years

setwd("D:/HabitatChange/SurvivalAnalysis")
hablist <- paste0("D:/landsat/outputmap/dense/outmap",1:46,".tif")
ras <- lapply(hablist, rast)
ras <- rast(ras)
vp <- rast("D:/HabitatChange/HabitatPattern/validpixels.tif")  # valid pixels created in P3_1
ras <- crop(ras, vp, mask=TRUE)

df <- as.data.frame(ras, xy=TRUE)
firstgrass <- function(x){
  if (3 %in% x == FALSE) {return(NA)}
  else {return(min(which(x==3)))}
}
firstshrubwood <- function(x){
  x[x==1] <- 2
  if (2 %in% x == FALSE) {return(NA)}
  else {return(min(which(x==2)))}
}
firstshrub <- function(x){
  if (2 %in% x == FALSE) {return(NA)}
  else {return(min(which(x==2)))}
}
firstwood <- function(x){
  if (1 %in% x == FALSE) {return(NA)}
  else {return(min(which(x==1)))}
}

df$firstgrass <- apply(df[,3:48], 1, firstgrass)  # year of first grass exist
df$firstshrubwood <- apply(df[,3:48], 1, firstshrubwood)  # year of first shrub/wood exist
df$firstshrub <- apply(df[,3:48], 1, firstshrub) # year of first shrub exist
df$firstwood <- apply(df[,3:48], 1, firstwood) # year of first wood exist

df_gs <- df[,c("x","y","firstgrass","firstshrubwood")]
df_sw <- df[,c("x","y","firstshrub","firstwood")]

df_gs <- df_gs[!is.na(df_gs$firstgrass),]   # filter grass exist
df_gs$event <- !is.na(df_gs$firstshrubwood)   # find censored (not change)
df_gs["firstshrubwood"][is.na(df_gs["firstshrubwood"])] <- 46  # convert censored to last year
df_gs$time <- df_gs$firstshrubwood - df_gs$firstgrass
df_gs <- df_gs[df_gs$time>0,]
df_gs$type <- "GrassToShrub"

df_sw <- df_sw[!is.na(df_sw$firstshrub),]   # filter shrub exist
df_sw$event <- !is.na(df_sw$firstwood)   # find censored (not change)
df_sw["firstwood"][is.na(df_sw["firstwood"])] <- 46  # convert censored to last year
df_sw$time <- df_sw$firstwood - df_sw$firstshrub
df_sw <- df_sw[df_sw$time>0,]
df_sw$type <- "ShrubToWood"

df_gw <- df[,c("x","y","firstgrass","firstwood")]
df_gw <- df_gw[!is.na(df_gw$firstgrass),]   # filter grass exist
df_gw <- df_gw[!is.na(df_gw$firstwood),]   # filter wood exist
df_gw$time <- df_gw$firstwood - df_gw$firstgrass
df_gw <- df_gw[df_gw$time>0,]  # grass in early years and wood in later years

df1 <- rbind(df_gs[,5:7],df_sw[,5:7])
s3 <- survfit2(Surv(time, event) ~ type, data = df1)
quantile(s3, probs = c(0.25, 0.5)) # GS: 7, 21; SW: 10, 29
ggsurvfit(s3)

p <- ggsurvfit(s3, linewidth=1, type = "risk") +
  labs(title = "Number of years for transition between vegetation classes",
       x = "Year", y = "Overall survival probability") + 
  scale_color_manual(labels = c("Transition from grassland\nto shrubland",
                                "Transition from shrubland\nto woodland"),
                     values = c('#55FF00', '#267300')) +
  add_quantile(y_value = 0.5, linetype = 2) +
  add_quantile(y_value = 0.25, linetype = 3) +
  scale_ggsurvfit(y_scales = list(limits=c(0,1))) +
  theme(plot.title = element_text(size=12, hjust=0.5), plot.title.position = "plot")

p
ggsave("survivalplot.pdf", width=12, height=9.5, units="cm")


##### Write grass-to-shrub-year and shrub-to-wood-year as raster files

setwd("D:/HabitatChange/SurvivalAnalysis")
vp <- rast("D:/HabitatChange/HabitatPattern/validpixels.tif")  # valid pixels created in P3_1

ras_gs <- vect(df_gs, geom=c("x","y"), crs=crs(vp))
ras_gs <- ras_gs[ras_gs$firstgrass<=20,]  # use only pixels start year <= 20
ras_gs <- rasterize(ras_gs, vp, "time")
writeRaster(ras_gs, "grasstoshrubyear.tif")

ras_sw <- vect(df_sw, geom=c("x","y"), crs=crs(vp))
ras_sw <- ras_sw[ras_sw$firstshrub<=20,]  # use only pixels start year <= 20
ras_sw <- rasterize(ras_sw, vp, "time")
writeRaster(ras_sw, "shrubtowoodyear.tif")


##### Environmental Variables

setwd("D:/HabitatChange/GISData")
Elev <- rast("NaturalEnv/dtm.tif")
Slope <- rast("NaturalEnv/slope.tif")
North <- rast("NaturalEnv/northness.tif")
TWI <- rast("NaturalEnv/TWI.tif")
Insol <- rast("NaturalEnv/insolation.tif")
Temp <- rast("NaturalEnv/temperature.tif")
Precip <- rast("NaturalEnv/precipitation.tif")
Wind <- rast("NaturalEnv/windspeed.tif")
DistCoast <- rast("NaturalEnv/DistfromCoast.tif")
CEC <- rast("NaturalEnv/SoilCEC.tif")
SOM <- rast("NaturalEnv/SoilOrganicMatter.tif")
ProxBuilt <- rast("HumanAct/DistfromBuildings_10periods.tif")[[10]]
CounPark <- rast("HumanAct/CountryParksSpecialAreas_10periods.tif")[[10]]
Plant <- rast("HumanAct/Plantation_union4periods.tif")
DistWood <- rast("WoodSpatial/DistfromWoodland_10periods.tif")[[1]]
AreaWood <- rast("WoodSpatial/AreaNearestWoodland_10periods.tif")[[1]]
Fire <- app(rast("Disturbance/Hillfire_10periods.tif"),sum)
Landslide <- app(rast("Disturbance/Landslide_10periods.tif"),sum)
Typhoon <- app(rast("Disturbance/Typhoon_10periods.tif"),mean)

envras <- c(Elev,Slope,North,TWI,Insol,Temp,Precip,
            Wind,DistCoast,CEC,SOM,
            ProxBuilt,CounPark,Plant,
            DistWood,AreaWood,Fire,Landslide,Typhoon)

names(ras_gs) <- "GrassShrub"
names(ras_sw) <- "ShrubWood"

ras <- c(ras_gs, ras_sw, crop(envras, ras_gs))
names(ras) <- c("GrassShrub", "ShrubWood", 
                "Elev","Slope","North","TWI","Insol","Temp","Precip",
                "Wind","DistCoast","CEC","SOM",
                "ProxBuilt","CounPark","Plant",
                "DistWood","AreaWood","Fire","Landslide","Typhoon")
df <- as.data.frame(ras)
df <- df[!(is.na(df$GrassShrub)&is.na(df$ShrubWood)),]
df <- df[!(is.na(df$CEC)|is.na(df$SOM)),]

setwd("D:/HabitatChange/SurvivalAnalysis")
write.csv(df,"successiontime_envras.csv",row.names=FALSE)


##### Correlation Analysis

library(pcaPP)
setwd("D:/HabitatChange/SurvivalAnalysis")
df <- read.csv("successiontime_envras.csv")
df_gs <- df[!is.na(df$GrassShrub), ]
df_sw <- df[!is.na(df$ShrubWood), ]

df_kendall <- data.frame()
for (i in 3:21){
  r1 <- cor.fk(df_gs$GrassShrub, df_gs[,i])
  r2 <- cor.fk(df_sw$ShrubWood, df_sw[,i])
  df_k1 <- data.frame(id=i, GrassShrub_k=r1, ShrubWood_k=r2)
  df_kendall <- rbind(df_kendall, df_k1)
}

df_kendall$variable <- colnames(df)[3:21]
write.csv(df_kendall, "df_kendall.csv")

lm_gs <- lm(GrassShrub~.,df_gs[,c(1,3:21)])
summary(lm_gs)
lm_sw <- lm(ShrubWood~.,df_sw[,c(2,3:21)])
summary(lm_sw)

