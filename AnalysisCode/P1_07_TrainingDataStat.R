library(data.table)
library(dplyr)

# pixel-level fusion

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

df <- traindata %>% group_by(Class) %>%
  summarise_at(predictors_withsd, c(mean,sd), na.rm = TRUE)
write.csv(df, "trainstat_pixelfusion.csv")


# difference by sensor

setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_t.csv")
trainpt <- read.csv("trainpt.csv")

trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                            variable.name = "yearclass", value.name = "Class")
trainpt_m$Class <- factor(trainpt_m$Class, levels = c(1:6))
levels(trainpt_m$Class) <- c("c1", "c2", "c3", "c4", "c5", "c6")
traindata <- merge(imgdf, trainpt_m, by = c("PointID","yearclass"))
traindata$row_num <- seq.int(nrow(traindata))
traindata$sensor <- substr(traindata$imgname, 20, 21)

predictors <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")
predictors <- c("C","B","G","R","RE","NIR","SW1","SW2")

df <- traindata %>% group_by(sensor, Class) %>%
  summarise_at(predictors, c(mean,sd), na.rm = TRUE)
write.csv(df, "trainstat_sensor.csv")


