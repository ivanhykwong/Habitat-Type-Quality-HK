library(terra)

# read probability layers created from RF model
rasdir <- "D:/landsat/outputmap/prob"
setwd(rasdir)
probmapname <- list.files(rasdir, pattern="tif$")
probmaplist <- lapply(probmapname, rast)

# Extract valid pixels
setwd("D:/HabitatChange/HabitatPattern")
vp <- rast("validpixels.tif")

rastodf <- function(ras){
  ras <- crop(ras, vp, mask=TRUE)
  df <- as.data.frame(ras, xy=TRUE)   # convert to dataframe, each layer become a column (probability value of the class)
  colnames(df) <- c("x","y","wood","shrub","grass","barren","builtup","water")
  return(df)
}
df_list <- lapply(probmaplist, rastodf)

yc <- seq(1975,2020,5)
for (i in 1:length(yc)){  # add year
  df1 <- df_list[[i]]
  df1$year <- yc[i]
  df_list[[i]] <- df1
}
df <- do.call(rbind, df_list)

# Perform PCA and normalise to forest index
pc <- prcomp(df[,3:8], center=TRUE, scale.=FALSE)
summary(pc)  # Proportion of Variance
pc$rotation  # Loadings
pcscore <- cbind(df[,c("x","y","year")],PC1=pc$x[,1])  # data frame with x,y coordinates and PC1 values
globalmin <- min(pcscore$PC1)   # -0.6836
globalmax <- max(pcscore$PC1)   # 0.5789
pcscore$forestindex <- (pcscore$PC1-globalmin)/(globalmax-globalmin) # 0 to 1

# create forest index raster
tempras <- rast("D:/HabitatChange/GISData/NaturalEnv/dtm.tif")
y_list <- unique(pcscore$year)
for (i in 1:10){
  y <- y_list[i]
  pc_y <- pcscore[pcscore$year==y,]
  pc_y <- vect(pc_y, geom=c("x","y"), crs=crs(tempras))
  pc_ras_y <- rasterize(pc_y, tempras, field="forestindex")
  if (i==1) {pc_ras <- pc_ras_y}
  else {pc_ras <- c(pc_ras, pc_ras_y)}
}
plot(pc_ras)
names(pc_ras) <- paste0("c",seq(1975,2020,5))
writeRaster(pc_ras, "ForestIndex_10periods.tif")

# Compare to field survey points (year 2010 and 2020)
setwd("D:/HabitatChange/HabitatPattern")
ras <- rast("ForestIndex_10periods.tif")

pt <- read.csv("referencedata/FieldPoint_ForestIndex.csv")
pt <- vect(pt, geom=c("Rect_X", "Rect_Y"), crs=crs(ras))
pt2010 <- pt[pt$Year==2010,]
pt2010 <- extract(ras["c2010"], pt2010, bind=TRUE)
pt2010 <- as.data.frame(pt2010)[,c("Habitat","c2010")]
pt2010 <- na.omit(pt2010)
colnames(pt2010) <- c("Habitat","index")
pt2010$type <- "field2010"

pt2020 <- pt[pt$Year==2020,]
pt2020 <- extract(ras["c2020"], pt2020, bind=TRUE)
pt2020 <- as.data.frame(pt2020)[,c("Habitat","c2020")]
pt2020 <- na.omit(pt2020)
colnames(pt2020) <- c("Habitat","index")
pt2020$type <- "field2020"

# Compare to output class
deskpt <- vect("D:/HabitatChange/SEM/randompt.shp")  # random points created in P3_3
deskpt_df <- as.data.frame(deskpt)

# extract forest index
for (i in 1:10){
  ras_y <- ras[[i]]
  pt_y <- deskpt[deskpt$year==i,]
  v_y <- extract(ras_y, pt_y, ID=FALSE)
  colnames(v_y) <- "forestindex"
  if (i==1) {v <- v_y}
  else {v <- rbind(v, v_y)}
}
deskpt_df <- cbind(deskpt_df,v)

# extract habitat map class
maplist <- list.files("D:/landsat/outputmap/", ".tif$", full.names=TRUE)
maplist <- lapply(maplist, rast)
habras <- rast(maplist)
names(habras) <- paste0("c",seq(1975,2020,5))
for (i in 1:10){            
  ras_y <- habras[[i]]
  pt_y <- deskpt[deskpt$year==i,]
  v_y <- extract(ras_y, pt_y, ID=FALSE)
  colnames(v_y) <- "habitatclass"
  if (i==1) {v <- v_y}
  else {v <- rbind(v, v_y)}
}
deskpt_df <- cbind(deskpt_df,v)

# combine field points and desktop points into a single data frame for plotting
pt <- data.frame(Habitat=deskpt_df$habitatclass, index=deskpt_df$forestindex)
pt$Habitat <- factor(pt$Habitat, levels=c(3,2,1,0), labels=c("Woodland","Shrubland","Grassland","Barren land"))
deskpt <- na.omit(pt)
deskpt$type <- "desktop"
set.seed(1)
deskpt_samp <- deskpt[sample(nrow(deskpt), 1000),]
df <- rbind(deskpt_samp, pt2010, pt2020)
set.seed(1)
df <- df[sample(1:nrow(df)), ]

# create jitter plot
df$type <- factor(df$type, levels=c("field2020","field2010","desktop"))
df$Habitat <- factor(df$Habitat, levels=c("Woodland","Shrubland","Grassland","Barren land",
                                          "Lowland forest","Rural plantation","Plantation",
                                          "Mixed shrubland","Woody shrubland","Shrubby grassland",
                                          "Bare rock","Bare rock/ soil","Mixed barren land"))
HabitatColour <- c("#267300","#55ff00","#d3e645","#e1e1e1",
                   "#267300","#007d7d","#007d7d",
                   "#38a800","#38a800","#d0ff73",
                   "#e1e1e1","#e1e1e1","#ffebaf")

library(ggplot2)
p <- ggplot(df, aes(type, index)) +
  geom_jitter(aes(colour = Habitat), size=0.8, width=0.25, height=0) +
  scale_color_manual(values=HabitatColour) +
  scale_x_discrete(labels = c("Field-surveyed\npoints (2020-21)","Field-surveyed\npoints (2009)","Habitat maps\n(1973-2022)")) +
  xlab("") + ylab("Forest index") +
  coord_flip() +
  theme_classic() +theme(axis.text=element_text(color="black"))
p
ggsave("D:/HabitatChange/HabitatPattern/referencedata/CompareFieldAndDesktop.pdf",
       width = 16, height = 7, units = "cm")


