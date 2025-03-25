# standardise dataframe

setwd("D:/HabitatChange/SEM")
dat <- read.csv("pt_df_data_8periods_withfa.csv")

contcol <- c("yrsinb","ProxBuilt","DistWood","AreaWood","Typhoon",
         "ForIndCur","HabQuaCur","ForIndNext","HabQuaNext") # continuous variables
dat[contcol] <- scale(dat[contcol])

# construct SEM model
library(piecewiseSEM)

model <- psem(
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
  lm(HabQuaNext~yrsinb+HabQuaCur+envfa1+envfa2+envfa3+ProxBuilt+CounPark+PlantCur+PlantPre+DistWood+AreaWood+FireCur+Landslide+Typhoon,dat),
  
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

s <- summary(model)

# write model output

options(max.print=10000)
sink(file = "psem_output.txt")
s
sink(file = NULL)
options(max.print=1000)

coefs(model, standardize = "scale")
r2 <- s$R2[9:12,c("Response","R.squared")]

# compute direct and indirect effects

library(semEff)
library(snow)
library(boot)

model.boot <- bootEff(model, R=1)
model.semeff <- semEff(model.boot)

# write semEff output

options(max.print=10000)
sink(file = "modelsemeff.txt")
model.semeff
summary(model.semeff)
sink(file = NULL)
options(max.print=1000)

# convert semEff output to a dataframe

a <- model.semeff$Summary
df1 <- data.frame()
for (i in 2:length(a)){
  df <- a[[i]]
  df$Response <- names(a[i])
  df <- df[,c(14,1,2,4,6,8,10,11,13)]
  colnames(df) <- c("Response","Category","Predictor","Effect","Bias","Std.Err.","LowerCI","UpperCI","sig")	
  df <- df[2:nrow(df),]
  df <- df[substr(df$Predictor,1,1) != " ",]
  df1 <- rbind(df1,df)
}
write.csv(df1,"modelsemeff_coef.csv")

