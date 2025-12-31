library(terra)
library(dplyr)
library(reshape2)

# Get area of each class

imgfolder <- "D:/landsat/outputmap"
setwd(imgfolder)
imglist <- list.files(imgfolder, pattern = "tif$")

for (i in 1:length(imglist)){
  imgname <- imglist[i]
  ras <- rast(imgname)
  df1 <- freq(ras)
  df1$area <- df1$count * res(ras)[1] * res(ras)[2] / (1000^2)
  df1 <- df1[,c("value","area")]
  colnames(df1) <- c("class", substr(imgname,8,16))
  if (i==1){df <- df1}
  else {df <- merge(df, df1, by="class")}
}
write.csv(df, "D:/landsat/outputmap/predictedarea.csv")

# Get area CI

trainpt <- read.csv("D:/landsat/train/trainpt.csv")
cv_pred <- readRDS("D:/landsat/train/cv_pred.rds")
pred_df <- bind_rows(cv_pred[c(3,7,11,15)])
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))
voteresult <- function(pred_df, temporalsmooth = TRUE){
  if (inherits(pred_df, "list")){pred_df <- pred_df[[1]]}
  # combine vote
  pred_df_vote <- pred_df %>%
    group_by(PointID, yearclass, repeats) %>%
    summarise(s1=sum(c1),s2=sum(c2),s3=sum(c3),s4=sum(c4),s5=sum(c5),s6=sum(c6))
  pred_df_vote$probsum <- apply(pred_df_vote[4:9], 1, sum)
  for (i in 4:9){pred_df_vote[,i] <- pred_df_vote[,i]/pred_df_vote$probsum}
  pred_df_vote <- as.data.table(pred_df_vote[,1:9])
  # temporal smoothing
  if (temporalsmooth==TRUE){
    tempsmooth <- function(pd_df){
      pd_df_smooth <- data.table()
      for (i in unique(pd_df$PointID)){
        d1 <- pd_df[pd_df$PointID==i,]
        if (nrow(d1) < 10){
          allyear <- paste0("c",seq(1975,2020,5))
          missingyear <- allyear[!allyear %in% d1$yearclass]
          for (m in missingyear){
            d1 <- rbind(d1, data.table(PointID=i, yearclass=m, repeats=r,
                                       s1=0,s2=0,s3=0,s4=0,s5=0,s6=0))
          }
        }
        d1 <- d1[order(d1$yearclass),]
        d1_dec <- d1[order(d1$yearclass, decreasing = TRUE),]
        p1 <- apply(d1_dec[,4:9], 2, WMA, n=2, wts=c(1,2))[10,, drop=FALSE]
        p2 <- apply(d1[,4:9], 2, WMA, n=3, wts=c(1,2,1))[3:10,, drop=FALSE]
        p3 <- apply(d1[,4:9], 2, WMA, n=2, wts=c(1,2))[10,, drop=FALSE]
        d1 <- cbind(d1[,1:3], rbind(p1,p2,p3))
        pd_df_smooth <- rbind(pd_df_smooth, d1)
      }
      return(pd_df_smooth)
    }
    pred_df_vote_smooth <- data.table()
    for (r in unique(pred_df_vote$repeats)){
      pred_df_vote_s1 <- pred_df_vote[pred_df_vote$repeats==r,]
      pred_df_vote_s1 <- tempsmooth(pred_df_vote_s1)
      pred_df_vote_smooth <- rbind(pred_df_vote_smooth, pred_df_vote_s1)
    }
    pred_df_vote <- pred_df_vote_smooth
  }
  # Get class
  pred_df_vote$predclass <- apply(pred_df_vote[,4:9], 1, which.max)
  pred_df_vote$predclass <- factor(pred_df_vote$predclass, levels = c(1:6))
  levels(pred_df_vote$predclass) <- c("c1", "c2", "c3", "c4", "c5", "c6")
  # join ref class
  trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                              variable.name = "yearclass", value.name = "refclass")
  trainpt_m$refclass <- factor(trainpt_m$refclass, levels = c(1:6))
  levels(trainpt_m$refclass) <- c("c1", "c2", "c3", "c4", "c5", "c6")
  pred_df_vote <- merge(pred_df_vote, trainpt_m, by = c("PointID","yearclass"))
  return(pred_df_vote)
}
pred_df_vote <- voteresult(pred_df)

getareaCI <- function(pred_df_vote, yearclass1, imgname){
  # accuracy assessment
  # https://blogs.fu-berlin.de/reseda/area-adjusted-accuracies/
  # area adjusted accuracy assessment (as in Olofsson et al. 2014)
  refpt <- pred_df_vote[pred_df_vote$yearclass==yearclass1,]
  confmat <- table(refpt$predclass, refpt$refclass)/5
  nclass <- length(unique(refpt$refclass))
  ras <- rast(imgname)
  df1 <- freq(ras)
  df1$area <- df1$count * res(ras)[1] * res(ras)[2] / (1000^2)
  maparea <- as.numeric(df1$area)
  conf <- 1.96
  A <- sum(maparea)
  W_i <- maparea / A
  n_i <- rowSums(confmat) 
  p <- W_i * confmat / n_i
  p[is.na(p)] <- 0
  p_area <- colSums(p) * A
  p_area_CI <- conf * A * sqrt(colSums((W_i * p - p ^ 2) / (n_i - 1))) 
  return(p_area_CI)
}

yearclasslist <- paste0("c",seq(1975,2020,5))
areaCIdf <- data.frame(class=1:6)
for (i in 1:10){
  areaCI <- getareaCI(pred_df_vote, yearclasslist[i], imglist[i])
  areaCIdf <- cbind(areaCIdf, areaCI=areaCI)
}
write.csv(areaCIdf, "D:/landsat/outputmap/predictedarea_CI.csv")

