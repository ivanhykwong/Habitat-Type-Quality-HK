# Load data and variables

setwd("D:/HabitatChange/SEM")
dat <- read.csv("pt_df_data.csv")
dat <- dat[dat$year%in%(2:9),] # remove period 1 and 10

# EFA of variables related to natural environment

library(psych)
col <- c("Elev","Slope","North","TWI","Temp","Precip","Insol","Wind","DistCoast","CEC","SOM")
envdf <- dat[,col]
envdf[col] <- scale(envdf[col]) # standardise to mean 0 and sd 1
summary(envdf)

envdf_fa <- fa(envdf, nfactors=3, rotate="varimax", fm="ml")
envdf_fa          
envdf_fa$loadings

# write factor summary and loadings

sink(file = "EFA_result.txt")
envdf_fa$loadings # cutoff = 0.1
summary(envdf_fa)
envdf_fa
sink(file = NULL)

# append factor scores to original dataframe

envfa <- as.data.frame(envdf_fa$scores)
colnames(envfa) <- c("envfa1","envfa2","envfa3")
dat <- cbind(dat, envfa)
write.csv(dat, "pt_df_data_8periods_withfa.csv", row.names=FALSE)


# convert factor scores to raster layer

library(terra)

setwd("D:/HabitatChange/SEM")
dat <- read.csv("pt_df_data_8periods_withfa.csv")
dat <- dat[dat$year==2, c("x","y","envfa1","envfa2","envfa3")]

vp <- rast("D:/HabitatChange/HabitatPattern/validpixels.tif")
pt <- vect(dat, geom=c("x","y"), crs=crs(vp))

fa1 <- interpIDW(vp, pt, "envfa1",radius=1000, maxPoints=10)
fa1 <- mask(fa1, vp)
fa2 <- interpIDW(vp, pt, "envfa2",radius=1000, maxPoints=10)
fa2 <- mask(fa2, vp)
fa3 <- interpIDW(vp, pt, "envfa3",radius=1000, maxPoints=10)
fa3 <- mask(fa3, vp)
fa_ras <- c(fa1,fa2,fa3)

writeRaster(fa_ras, "envfa_raster.tif")

