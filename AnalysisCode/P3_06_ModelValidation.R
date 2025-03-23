# standardise dataframe

setwd("D:/HabitatChange/SEM")
dat <- read.csv("pt_df_data_8periods_withfa.csv")

contcol <- c("yrsinb","ProxBuilt","DistWood","AreaWood","Typhoon",
             "ForIndCur","HabQuaCur","ForIndNext","HabQuaNext") # continuous variables
dat[contcol] <- scale(dat[contcol])

# manually test lm assumption for each element in modellist
# subset 1000 rows for test and plots

dat1 <- dat
set.seed(2024)
dat <- dat[sample(nrow(dat), 1000), ]
table(dat$year)

modellist <- list(
  lm(ProxBuilt~yrsinb+envfa1+envfa2+envfa3,dat),
  glm(CounPark~yrsinb+envfa1+envfa2+envfa3,'binomial',dat),
  glm(PlantCur~yrsinb+envfa1+envfa2+envfa3,'binomial',dat),
  glm(PlantPre~yrsinb+envfa1+envfa2+envfa3,'binomial',dat),
  
  glm(FireCur~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,'binomial',dat),
  glm(Landslide~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,'binomial',dat),  
  lm(Typhoon~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,dat),
  
  lm(ForIndCur~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  lm(HabQuaCur~yrsinb+ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  
  lm(ForIndNext~yrsinb+ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  lm(HabQuaNext~yrsinb+HabQuaCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat)
)

# collinearity
### all correlation < 0.5
### all VIF < 2
### remove HabQuaCur in ForIndNext~
### remove ForIndCur and ForIndNext in HabQuaNext~

library(car)
library(gvlma)
for (i in 1:11){
  print(i)
  model <- modellist[[i]]
  print(vif(model))
}

# lm plots

model <- modellist[[8]]
summary(model)
vif(model)
plot(model, 1)
durbinWatsonTest(model)  # test existence of autocorrelation in residuals; stat close to 2 and p > 0.05; cannot reject -> no autocorrelation
plot(model, 2)
plot(model, 3)
ncvTest(model)  # test for heteroscedasticity; p > .05, suggesting that our data is homoscedastic. Yay!
plot(model, 4)
plot(model, 5)  # not cross Cook’s distance line, so we’re okay
plot(model, c(1:3, 5))
gvlma(model)
par(mfrow=c(2,2))
plot(model)

# spatial autocorrelation
# https://forum.posit.co/t/morans-i-for-point-data/180588/2

library(ape)
library(dplyr)
library(terra)
model <- lm(ForIndCur~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat)
lmpred <- predict(model)
lmpred <- data.frame(x=dat$x, y=dat$y, t=dat$year, e=lmpred - dat$ForIndCur)  # residual in each observation
pt <- lmpred %>% group_by(x,y) %>% summarise(e=mean(e))  # summarise residuals by locations
pt$x <- as.numeric(pt$x)
pt$y <- as.numeric(pt$y)
vp <- rast("D:/HabitatChange/HabitatPattern/validpixels.tif")
vp <- aggregate(vp, fact=30)
pt <- vect(pt, geom=c("x","y"), crs=crs(vp))  # convert to spatial points and raster
writeVector(pt, "residual_pt.shp")
ras <- rasterize(pt, vp, "e", fun="mean")
writeRaster(ras, "residual_raster.tif")

pt <- as.data.frame(ras, xy=T)  # calculate Moran's I test
colnames(pt) <- c("x","y","e")
pt.dists <- as.matrix(dist(cbind(pt$x, pt$y)))
pt.dists.inv <- 1/pt.dists
diag(pt.dists.inv) <- 0
Moran.I(pt$e, pt.dists.inv)

# temporal autocorrelation

time_resi <- lmpred %>% group_by(t) %>% summarise(mean=mean(e), sd=sd(e))
plot(time_resi$t, time_resi$mean)
write.csv(time_resi, "residual_time.csv")


# Cross-validated R square using observation in other periods

dat <- dat1
m1 <- "ForIndCur~yrsinb+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon"
m2 <- "HabQuaCur~yrsinb+ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon"
m3 <- "ForIndNext~yrsinb+ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon"
m4 <- "HabQuaNext~yrsinb+HabQuaCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon"

ysamp <- 2:9
df1 <- data.frame()
df2 <- data.frame()
df3 <- data.frame()
df4 <- data.frame()
for (y in ysamp) {  # 8 fold based on 8 periods
  print(y)
  dftrain <- dat[dat$year != y,]  # train using other years
  dftest <- dat[dat$year == y,]   # test using one year
  mod1 <- lm(m1, dftrain)
  mod1_pred <- predict(mod1, dftest)
  df1 <- rbind(df1, data.frame(obs=dftest$ForIndCur, est=mod1_pred))
  mod2 <- lm(m2, dftrain)
  mod2_pred <- predict(mod2, dftest)
  df2 <- rbind(df2, data.frame(obs=dftest$HabQuaCur, est=mod2_pred))
  mod3 <- lm(m3, dftrain)
  mod3_pred <- predict(mod3, dftest)
  df3 <- rbind(df3, data.frame(obs=dftest$ForIndNext, est=mod3_pred))
  mod4 <- lm(m4, dftrain)
  mod4_pred <- predict(mod4, dftest)
  df4 <- rbind(df4, data.frame(obs=dftest$HabQuaNext, est=mod4_pred))
}

cor(df1$est, df1$obs)**2  # ForIndCur: 0.426
cor(df2$est, df2$obs)**2  # HabQuaCur: 0.704
cor(df3$est, df3$obs)**2  # ForIndNext: 0.954
cor(df4$est, df4$obs)**2  # HabQuaNext: 0.860


