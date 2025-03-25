# load SEM data

library(piecewiseSEM)
setwd("D:/HabitatChange/SEM")
dat <- read.csv("pt_df_data_8periods_withfa.csv")

contcol <- c("yrsinb","ProxBuilt","DistWood","AreaWood","Typhoon",
             "ForIndCur","HabQuaCur","ForIndNext","HabQuaNext") # continuous variables
dat[contcol] <- scale(dat[contcol])

# construct SEM model
# "year since baseline" variable not included in this model
model <- psem(
  lm(ProxBuilt~envfa1+envfa2+envfa3,dat),
  glm(CounPark~envfa1+envfa2+envfa3,'binomial',dat),
  glm(PlantCur~envfa1+envfa2+envfa3,'binomial',dat),
  glm(PlantPre~envfa1+envfa2+envfa3,'binomial',dat),
  
  glm(FireCur~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,'binomial',dat),
  glm(Landslide~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,'binomial',dat),  
  lm(Typhoon~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,dat),
  
  lm(ForIndCur~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  lm(HabQuaCur~ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  
  lm(ForIndNext~ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  lm(HabQuaNext~HabQuaCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  
  CounPark %~~% ProxBuilt,
  PlantCur %~~% ProxBuilt,
  PlantPre %~~% ProxBuilt,
  PlantCur %~~% CounPark,
  PlantPre %~~% CounPark,
  PlantPre %~~% PlantCur,
  AreaWood %~~% DistWood,
  FireCur %~~% Landslide,
  FireCur %~~% Typhoon,
  Landslide %~~% Typhoon
)

# create multi-group SEM using year

pmultigroup <- multigroup(model, group = "year")

# write model output

options(max.print=10000)
sink(file = "coef_sepyear_print.txt")
pmultigroup
sink(file = NULL)
options(max.print=1000)

# convert to data frame

pcoefs <- pmultigroup$group.coefs
for (i in 1:length(pcoefs)){
  cf <- pcoefs[[i]]
  constrained <- cf[,10]
  cf_est <- cf[,c("Response","Predictor","Std.Estimate")]
  colnames(cf_est) <- c("Response","Predictor",paste0("year",i+1))
  if (i==1) {
    cf_est1 <- cbind(constrained, stat="StdEstimate", cf_est)
  } else {
    cf_est1 <- merge(cf_est1, cf_est, by=c("Response","Predictor"))
  }
}
df <- cf_est1

# find monotonic trends

df$kendall <- NA
df$kendall_p <- NA
for (r in 1:nrow(df)){
  df1 <- as.numeric(df[r,5:12])
  kendall <- cor.test(df1, 2:9, method="kendall")
  df[r,13] <- kendall$estimate
  df[r,14] <- kendall$p.value
}
write.csv(df, "coef_sepyear.csv")

# Confidence interval of the estimates

library(dplyr)
df1 <- data.frame()
for (y in 2:9){
  print(y)
  print(Sys.time())
  d1 <- dat[dat$year==y, ]
  for (i in 1:1000){  # repeat 1000 times
    if (i %% 100 == 0) {print(i);print(Sys.time())}
    d <- d1[sample(nrow(d1), 10000), ]   # subsample 10000 rows
    model <- psem(
      lm(ProxBuilt~envfa1+envfa2+envfa3,d),
      glm(CounPark~envfa1+envfa2+envfa3,'binomial',d),
      glm(PlantCur~envfa1+envfa2+envfa3,'binomial',d),
      glm(PlantPre~envfa1+envfa2+envfa3,'binomial',d),
      glm(FireCur~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,'binomial',d),
      lm(Typhoon~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood,d),
      lm(ForIndCur~envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Typhoon,d),
      lm(HabQuaCur~ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Typhoon,d),
      lm(ForIndNext~ForIndCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Typhoon,d),
      lm(HabQuaNext~HabQuaCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Typhoon,d),
      CounPark %~~% ProxBuilt,
      PlantCur %~~% ProxBuilt,
      PlantPre %~~% ProxBuilt,
      PlantCur %~~% CounPark,
      PlantPre %~~% CounPark,
      PlantPre %~~% PlantCur,
      AreaWood %~~% DistWood,
      FireCur %~~% Typhoon
    )
    co <- coefs(model)
    c1 <- co[(co$Response=="FireCur")&(co$Predictor=="AreaWood"),"Std.Estimate"]
    c2 <- co[(co$Response=="ForIndCur")&(co$Predictor=="envfa2"),"Std.Estimate"]
    c3 <- co[(co$Response=="ForIndCur")&(co$Predictor=="CounPark"),"Std.Estimate"]
    c4 <- co[(co$Response=="ForIndCur")&(co$Predictor=="PlantPre"),"Std.Estimate"]
    c5 <- co[(co$Response=="HabQuaCur")&(co$Predictor=="ProxBuilt"),"Std.Estimate"]
    c6 <- co[(co$Response=="HabQuaCur")&(co$Predictor=="PlantPre"),"Std.Estimate"]
    df2 <- data.frame(year=y, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6)
    df1 <- rbind(df1, df2)
  }
}
write.csv(df1, "coef_sepyear_10000subsamp_1000times.csv")
df2 <- df1 %>% group_by(year) %>% summarise(
  c1a=quantile(c1,0.025), c1b=mean(c1), c1c=quantile(c1,0.975),
  c2a=quantile(c2,0.025), c2b=mean(c2), c2c=quantile(c2,0.975),  
  c3a=quantile(c3,0.025), c3b=mean(c3), c3c=quantile(c3,0.975),
  c4a=quantile(c4,0.025), c4b=mean(c4), c4c=quantile(c4,0.975),
  c5a=quantile(c5,0.025), c5b=mean(c5), c5c=quantile(c5,0.975),
  c6a=quantile(c6,0.025), c6b=mean(c6), c6c=quantile(c6,0.975),
)
write.csv(df2, "coef_sepyear_10000subsamp_1000times_95ci.csv")

# plot coefficients against time

library(ggplot2)
library(reshape2)
setwd("D:/HabitatChange/SEM")
df <- read.csv("coef_sepyear_select.csv")
df <- df[,c("pos","stat","X1980","X1985","X1990","X1995","X2000","X2005","X2010","X2015")]
df <- melt(df, id.vars=c("pos","stat"), variable.name="year", value.name="value")
df <- dcast(df, pos + year ~ stat)
df$year <- as.integer(substr(df$year,2,5))

p <- ggplot(data=df, aes(x=year)) +
  geom_line(aes(y=StdEstimate)) + 
  geom_ribbon(aes(ymin=Low, ymax=Up), alpha=0.15) +
  facet_wrap(vars(pos), scales="free_y", nrow=2) +
  scale_x_continuous(name="Year") +
  scale_y_continuous(name="Standardized estimate") +
  theme_bw()
p
ggsave("coef_sepyear_lineplot.pdf",width=20,height=11,unit="cm") 

