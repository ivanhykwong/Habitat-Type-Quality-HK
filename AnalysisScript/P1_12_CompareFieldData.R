library(terra)
library(sf)
library(caret)
library(dplyr)
library(data.table)
library(reshape2)
library(ggplot2)

# Function to extract habitat class at field point locations

getmapclass <- function(ras, df, bufdist=60){
  validclass <- c("Woodland","Shrubland","Grassland","Barren land","Built-up area","Water")
  classdf <- data.frame(classnum=1:6,classnam=validclass)
  # create point from df
  colnames(df) <- c("x","y","class")
  df <- df[df$class %in% validclass,]
  df$id <- 1:nrow(df)
  refdf <- cbind(df[,c("id","x","y")], refclass=df[,"class"])
  pt <- vect(df, geom=c("x","y"), crs=crs(ras))
  pt_buff <- buffer(pt, bufdist)
  # extract value from ras
  pt_buff_ex <- extract(ras, pt_buff, fun=table)
  pt_buff_ex <- melt(pt_buff_ex, id.vars="ID", variable.name="classnum",value.name="count")
  pt_buff_ex <- pt_buff_ex[pt_buff_ex$count>0,]
  mapdf <- merge(pt_buff_ex, classdf, by="classnum")
  mapdf <- mapdf[,c("ID","count","classnam")]
  colnames(mapdf) <- c("id","count","mapclass")
  # merge with refclass
  mapdf <- merge(mapdf, refdf, by="id")
  mapdf$classequal <- mapdf$mapclass == mapdf$refclass
  mapdf <- mapdf[order(-mapdf$classequal, -mapdf$count),]
  mapdf <- mapdf[!duplicated(mapdf$id),]
  mapdf$id <- as.numeric(mapdf$id)
  mapdf <- mapdf[order(mapdf$id),]
  mapdf <- mapdf[,c("id","x","y","refclass","mapclass","count","classequal")]
  return(mapdf)
}

# Function to calculate accuracy (and CI) by comparing predicted and reference classes

getaccuracy <- function(mapdf, ras){
  validclass <- c("Woodland","Shrubland","Grassland","Barren land","Built-up area","Water")
  confmat <- table(factor(mapdf$mapclass, levels=validclass),
                   factor(mapdf$refclass, levels=validclass))
  OA <- sum(diag(confmat))/sum(confmat)
  PA <- diag(confmat)/colSums(confmat)
  UA <- diag(confmat)/rowSums(confmat)
  useclass <- colSums(confmat) != 0
  confmat <- confmat[useclass,useclass]
  PA_t <- PA[useclass]
  UA_t <- UA[useclass]
  
  nclass <- ncol(confmat)
  df1 <- freq(ras)
  df1$area <- df1$count * res(ras)[1] * res(ras)[2] / (1000^2)
  maparea <- as.numeric(df1$area)[useclass]
  conf <- 1.96
  A <- sum(maparea)
  W_i <- maparea / A
  n_i <- rowSums(confmat) 
  p <- W_i * confmat / n_i
  p[is.na(p)] <- 0
  PA <- diag(p) / colSums(p)
  UA <- diag(p) / rowSums(p)
  OA_CI <- conf * sqrt(sum(W_i ^ 2 * UA * (1 - UA) / (n_i - 1)))
  UA_CI <- conf * sqrt(UA * (1 - UA) / (n_i - 1)) 
  N_j <- sapply(1:nclass, function(x) sum(maparea / n_i * confmat[ , x]) )
  tmp <- sapply(1:nclass, function(x) sum(maparea[-x] ^ 2 * confmat[-x, x] / n_i[-x] * ( 1 - confmat[-x, x] / n_i[-x]) / (n_i[-x] - 1)) )
  PA_CI <- conf * sqrt(1 / N_j ^ 2 * (maparea ^ 2 * ( 1 - PA ) ^ 2 * UA * (1 - UA) / (n_i - 1) + PA ^ 2 * tmp))
  
  pred_df_a <- data.frame(class=names(PA), pa=PA_t, pa_ci=PA_CI, 
                          ua=UA_t, ua_ci=UA_CI)
  return(pred_df_a)
}

# Load points from 1995 and 2015 field surveys and compute accuracies

setwd("D:/landsat/fieldpoint")
df_1995_2015 <- read.csv("pt_1995_2015.csv")
df_1995 <- df_1995_2015[,c("Easting", "Northing", "Class1995")]
df_2015 <- df_1995_2015[,c("Easting", "Northing", "Class2015")]

ras_1995 <- rast("D:/landsat/outputmap/outmap_1993-1997.tif")
ras_2015 <- rast("D:/landsat/outputmap/outmap_2013-2017.tif")

mapdf_1995 <- getmapclass(ras_1995, df_1995, 120)
mapdf_1995$year <- 1995
ac_1995 <- getaccuracy(mapdf_1995, ras_1995)
ac_1995$year <- 1995

mapdf_2015 <- getmapclass(ras_2015, df_2015, 120)
mapdf_2015$year <- 2015
ac_2015 <- getaccuracy(mapdf_2015, ras_2015)
ac_2015$year <- 2015

# Load points from 2010 and 2020 field surveys and compute accuracies

df_2010 <- read.csv("pt_2010.csv")
df_2010 <- df_2010[,c("Rect_X", "Rect_Y", "RefClass")]
ras_2010 <- rast("D:/landsat/outputmap/outmap_2008-2012.tif")

mapdf_2010 <- getmapclass(ras_2010, df_2010)
mapdf_2010$year <- 2010
ac_2010 <- getaccuracy(mapdf_2010, ras_2010)
ac_2010$year <- 2010

df_2020 <- read.csv("pt_2020.csv")
df_2020 <- df_2020[,c("Rect_X", "Rect_Y", "RefClass")]
ras_2020 <- rast("M:/landsat/outputmap/outmap_2018-2022.tif")

mapdf_2020 <- getmapclass(ras_2020, df_2020)
mapdf_2020$year <- 2020
ac_2020 <- getaccuracy(mapdf_2020, ras_2020)
ac_2020$year <- 2020

mapdf_all <- rbind(mapdf_1995,mapdf_2010,mapdf_2015,mapdf_2020)
ac_all <- rbind(ac_1995,ac_2010,ac_2015,ac_2020)
write.csv(mapdf_all, "mapdf_all.csv")
write.csv(ac_all, "ac_all.csv")


# Compare with LiDAR height in 2010 and 2020

ras_2010 <- rast("D:/landsat/outputmap/outmap_2008-2012.tif")
ras_2020 <- rast("D:/landsat/outputmap/outmap_2018-2022.tif")
lidar_2010 <- rast("D:/landsat/referencedata/ndsm_2010.tif")
lidar_2020 <- rast("D:/landsat/referencedata/ndsm_2020.tif")

chmdf <- terra::as.data.frame(c(ras_2010, lidar_2010), na.rm=TRUE)
colnames(chmdf) <- c("class","height")
chmdf <- as.data.table(chmdf)
df_trim <- data.table()
for (i in 1:3){
  df1 <- chmdf[chmdf$class==i,]
  df1$percentile <- rank(df1$height, ties.method = "first")/length(df1$height)
  df1 <- df1[percentile %between% c(0.01, 0.99)]
  df_trim <- rbind(df_trim, df1)
}
df_trim$class <- as.factor(paste0("2010c",df_trim$class))

chmdf <- terra::as.data.frame(c(ras_2020, lidar_2020), na.rm=TRUE)
colnames(chmdf) <- c("class","height")
chmdf <- as.data.table(chmdf)
df_trim1 <- data.table()
for (i in 1:3){
  df1 <- chmdf[chmdf$class==i,]
  df1$percentile <- rank(df1$height, ties.method = "first")/length(df1$height)
  df1 <- df1[percentile %between% c(0.01, 0.99)]
  df_trim1 <- rbind(df_trim1, df1)
}
df_trim1$class <- as.factor(paste0("2020c",df_trim1$class))
df_trim <- rbind(df_trim, df_trim1)

p <- ggplot(df_trim, aes(x=class, y=height, fill=class)) + 
  stat_boxplot(geom="errorbar", width = 0.3) +
  geom_boxplot(color="black", outlier.color="black", outlier.shape = 1, 
               outlier.size=0.5, outlier.alpha=0.05) +
  scale_fill_manual(values=c("#267300", "#44CC00", "#84FF47", "#267300", "#44CC00", "#84FF47")) +
  scale_y_continuous(trans='sqrt', breaks=c(0,1,seq(2,20,by=2))) +
  theme_classic()
p

df_trim_stat <- df_trim %>% group_by(class) %>%
   summarise(min=min(height),q1=quantile(height,0.25),median=median(height),
             mean=mean(height),q3=quantile(height,0.75),max=max(height))
df_trim_stat

write.csv(df_trim,"HeightCompareLidar.csv")

