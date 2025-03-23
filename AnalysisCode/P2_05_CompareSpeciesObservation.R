# GBIF.org (22 January 2024) GBIF Occurrence Download 
# https://doi.org/10.15468/dl.eexe9q
# Records included: 322984 records from 1 published datasets
# Filter used: "Country is Hong Kong", "DatasetKey is iNaturalist Research-grade Observations"

# Filter observations (year 2018-2023, identify to species level, uncertainty<=500m)

setwd("D:/HabitatQuality/iNaturalist")
df <- read.table("0067764-231120084113126.csv", sep="\t", header=TRUE, row.names=NULL, quote="", fill = TRUE)
df <- df[,c("gbifID","year","verbatimScientificName","taxonRank","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters")]
df <- df[df$year >= 2018 & df$year <= 2023,]
df <- df[df$taxonRank %in% c("SPECIES","SUBSPECIES","VARIETY","FORM"),]
df <- df[!is.na(df$coordinateUncertaintyInMeters),]
df <- df[df$coordinateUncertaintyInMeters <= 500,]
df_gbif <- df

# Match with species checklist and conservation concern

specieslist <- read.csv("SynonymList.csv")
speciesconcern <- read.csv("SpeciesConservationConcern.csv")
df <- merge(df_gbif, specieslist, by.x="verbatimScientificName", by.y="Synonym")
df <- merge(df, speciesconcern, by="ChecklistName")
df <- df[,c("decimalLatitude", "decimalLongitude", "Taxa", "ChecklistName", "Concern")]
colnames(df) <- c("Lat","Lon","Taxa","Name","Concern")
write.csv(df, "iNaturalist_observations_filtered.csv")
table(speciesconcern$Taxa, speciesconcern$Concern)

# fishnet grid subset

library(terra)
setwd("D:/HabitatQuality/iNaturalist")
fishnet <- vect("FishnetGrid.shp")  # 1km x 1km grid polygon created in ArcGIS
habmap <- rast("D:/HabitatQuality/inputmap/habmap_2018-2022_8class.tif")
habmap <- ifel(is.na(habmap), 0, 1)
fishnet_prop <- zonal(habmap, fishnet, fun=mean, weights=TRUE)
colnames(fishnet_prop) <- "fishnet_prop"
fishnet <- cbind(fishnet, fishnet_prop)
fishnet_subset <- fishnet[fishnet$fishnet_prop >= 0.25,]
fishnet_subset$GridID <- fishnet_subset$Id
writeVector(fishnet_subset, "FishnetGrid_subset.shp")

# Match observation pts and grids (spatial join)

library(sf)
library(dplyr)
setwd("D:/HabitatQuality/iNaturalist")
fishnet <- st_read("FishnetGrid_subset.shp")
obsdf <- read.csv("iNaturalist_observations_filtered.csv")
obspt <- st_as_sf(obsdf, coords = c("Lon", "Lat"), crs="EPSG:4326")
obspt <- st_transform(obspt, st_crs(fishnet))
obspt <- st_join(obspt, fishnet)
obspt_df <- st_drop_geometry(obspt[,c("GridID","Taxa","Name","Concern")])
obspt_df <- obspt_df[!is.na(obspt_df$GridID),]

# summary table of species count

table(obspt_df$Taxa, obspt_df$Concern)
obspt_df %>% group_by(Taxa) %>% summarise(n_distinct(Name), n_distinct(Name[Concern=="Y"]))

# compute biodiversity variables

taxa <- c("All", "Animals", unique(obspt_df$Taxa))
f <- st_drop_geometry(fishnet)
for (t in taxa){
  if (t=="All") {obspt_df_taxa <- obspt_df}
  else if (t=="Animals") {obspt_df_taxa <- obspt_df[obspt_df$Taxa!="Plants",]}
  else {obspt_df_taxa <- obspt_df[obspt_df$Taxa==t,]}
  gridstat <- obspt_df_taxa %>% group_by(GridID) %>%
    summarise(pt.count = n(), sp.count = n_distinct(Name), 
              con.sp.count = n_distinct(Name[Concern=="Y"]),
              con.sp.prop = con.sp.count/sp.count, 
              con.pt.prop = sum(Concern=="Y")/pt.count)
  gridstat <- gridstat[gridstat$pt.count>=30,]  # remove grids with less than 30 points
  gridcolnames <- colnames(gridstat)[2:ncol(gridstat)]
  colnames(gridstat) <- c("GridID", paste0(t,".",gridcolnames))
  f <- merge(f, gridstat, by="GridID", all=TRUE)
}
head(f)

# Match HQ map and grids

library(terra)
setwd("D:/HabitatQuality/iNaturalist")
fishnet_v <- vect("FishnetGrid_subset.shp")
hqmap <- rast("D:/HabitatQuality/workspace/quality_c_y2020.tif")
hq_grid <- zonal(hqmap, fishnet_v, fun=mean, na.rm=TRUE, weights=TRUE)
colnames(hq_grid) <- "habitatquality"
hq_grid_df <- cbind(fishnet_v[,"GridID"], hq_grid)
m <- merge(hq_grid_df, f, by="GridID", all=TRUE)

# Correlation analysis

cormatrix <- cor(as.data.frame(m), method="pearson", use="pairwise.complete.obs")
write.csv(cormatrix[,2], "correlation.csv")

# p value

df <- as.data.frame(m)
dfcol <- colnames(df)
dfcol <- dfcol[6:length(dfcol)]
df_p <- data.frame()
for (i in dfcol){
  ct <- cor.test(df$habitatquality, df[,i], method="pearson", exact=FALSE)
  df_p1 <- data.frame(variable=i, p=ct$p.value)
  df_p <- rbind(df_p, df_p1)
}
df_p$sig <- df_p$p < 0.01
write.csv(df_p, "correlation_pvalue.csv")

# write shapefile

shp <- f
shp[is.na(shp)] <- -999
shp <- merge(hq_grid_df, shp, by="GridID", all.x=TRUE)
colm <- names(shp)
colm <- gsub("Dragonflies","Dr",colm)
colm <- gsub("Animals","An",colm)
colm <- gsub("Plants","Pl",colm)
colm <- gsub("Butterflies","Bu",colm)
colm <- gsub("Birds","Bi",colm)
colm <- gsub("Reptiles","Re",colm)
colm <- gsub("Amphibians","Am",colm)
colm <- gsub("Mammals","Ma",colm)
colm <- gsub("con","c",colm)
colm <- gsub(".","_",colm,fixed = TRUE)
names(shp) <- colm
writeVector(shp, "fishnet_value.shp")
write.csv(as.data.frame(shp), "fishnet_value_csv.csv")

# create scatterplot

library(dplyr)
library(reshape2)
library(ggplot2)

setwd("D:/HabitatQuality/iNaturalist")
dat <- read.csv("fishnet_value_csv.csv")
dat <- dat[dat$CountNA<72,]  # remove taxa with too few valid grids

# Count of species of conservation concern
df <- dat %>% select(c("habitatqua", contains("_c_sp_c")))
taxa <- c("All","Animals","Plants","Butterflies","Birds","Reptiles","Dragonflies","Amphibians","Mammals")
colnames(df) <- c("habitatqua",taxa)
df <- melt(df, id.vars="habitatqua", variable.name="taxa", value.name="spcount")
df <- df[!is.na(df$spcount),]
df$taxa <- factor(df$taxa, levels=c("All","Plants","Animals","Mammals","Birds",
                                    "Reptiles","Amphibians","Butterflies","Dragonflies")) # re-arrange sequence
levels(df$taxa) <- c("All taxa","Plants","All animals","Mammals","Birds",
                     "Reptiles","Amphibians","Butterflies","Dragonflies") # change display name
df["spcount"][df["spcount"] == 0] <- 0.5

p <- ggplot(data=df, aes(spcount, habitatqua)) +
  geom_smooth(method="lm", linewidth=0.8, color="gray40", fill="gray") +
  geom_point(size=0.8) + 
  facet_wrap(vars(taxa)) +
  scale_x_continuous(name="Count of species of conservation concerns", 
                     trans="log2",
                     breaks = c(0.5,1,2,4,8,16,32,64),
                     label = c(0,1,2,4,8,16,32,64)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name="Estimated habitat quality") +
  theme_bw()
p
ggsave("scatterplot_c_sp_c.pdf",width=20,height=15,unit="cm")
summary(df)
head(df)

# Proportion of points of species of conservation concern
df <- dat %>% select(c("habitatqua", contains("_c_pt_p")))
taxa <- c("All","Animals","Plants","Butterflies","Birds","Reptiles","Dragonflies","Amphibians","Mammals")
colnames(df) <- c("habitatqua",taxa)
df <- melt(df, id.vars="habitatqua", variable.name="taxa", value.name="ptprop")
df <- df[!is.na(df$ptprop),]
df$taxa <- factor(df$taxa, levels=c("All","Plants","Animals","Mammals","Birds",
                                    "Reptiles","Amphibians","Butterflies","Dragonflies")) # re-arrange sequence
levels(df$taxa) <- c("All taxa","Plants","All animals","Mammals","Birds",
                     "Reptiles","Amphibians","Butterflies","Dragonflies") # change display name

p <- ggplot(data=df, aes(ptprop, habitatqua)) +
  geom_smooth(method="lm", linewidth=0.8, color="gray40", fill="gray") +
  geom_point(size=0.8) + 
  facet_wrap(vars(taxa)) +
  scale_x_continuous(name="Proportion of observations of species of conservation concerns",
                     limits=c(0,0.9)) +
  coord_cartesian(ylim=c(0,1)) +
  scale_y_continuous(name="Estimated habitat quality") +
  theme_bw()
p
ggsave("scatterplot_c_pt_p.pdf",width=20,height=15,unit="cm")
summary(df)
head(df)

