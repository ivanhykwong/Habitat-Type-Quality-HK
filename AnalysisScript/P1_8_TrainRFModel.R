library(raster)
library(sf)
library(dplyr)
library(data.table)
library(reshape2)
library(caret)
library(randomForest)
library(foreach)
library(doParallel)
library(kit)
library(groupdata2)
library(e1071)
library(ranger)
library(matrixStats)
library(TTR)
library(tuneRanger)

# Prepare data

setwd("D:/landsat/train")
imgdf <- fread("imgdata_all_t.csv")
trainpt <- read.csv("trainpt.csv")

trainpt_m <- reshape2::melt(trainpt[,c(1,5:14)], id.vars = c("PointID"),
                            variable.name = "yearclass", value.name = "Class")
trainpt_m$Class <- factor(trainpt_m$Class, levels = c(1:6))
levels(trainpt_m$Class) <- c("c1", "c2", "c3", "c4", "c5", "c6")
traindata <- merge(imgdf, trainpt_m, by = c("PointID","yearclass"))
traindata$row_num <- seq.int(nrow(traindata))

predictors <- c("G","R","NIR","ndvi","gndvi","glcm_mean","glcm_con","Slope","Northness")

# define function: trainpt id -> row numbers (for select points in each fold)

getrn <- function(fold1, td){
  for (n in fold1) {fold1[n] <- trainpt$PointID[n]}
  rn <- td[td$PointID %in% fold1,]$row_num
  return(rn)
} # foldindex_row <- lapply(foldindex, getrn, traindata)

# Function to build RF model and perform 10-fold cross validation

buildmodel <- function(traindata, predictors, separatewetdry=TRUE, ncore=4, seed=2022){
  set.seed(seed)
  seeds <- sample.int(n=2000, 100)
  # separate wet dry
  if (separatewetdry == TRUE){
    td_wet <- traindata[traindata$monthclass == "wet"]
    td_wet$row_num <- seq.int(nrow(td_wet))
    td_dry <- traindata[traindata$monthclass == "dry"]
    td_dry$row_num <- seq.int(nrow(td_dry))
  } else {
    td_all <- traindata
    td_all$row_num <- seq.int(nrow(td_all))
  }
  
  # predict prob
  rawtoprob <- function(rawpred){
    rawpred <- rawpred$predictions
    df <- data.frame(matrix(ncol = 6, nrow = nrow(rawpred)))
    for (n in 1:6){df[,n] <- apply(rawpred==n, 1, sum)}
    df <- df/sum(df[1,])
    colnames(df) <- c("c1", "c2", "c3", "c4", "c5", "c6")
    return(df)
  }

  # Define ranger function
  doranger <- function(td1, predictors1, modelname){
    # tune mtry
    td_t <- as.data.table(balance(td1, "min", cat_col="Class"))
    res <- data.frame()
    for (s in 1:10){
      set.seed(seeds[s])
      res1 <- tuneMtryFast(Class ~., data=td_t[,c(..predictors1,"Class")],
                           num.treesTry=200, stepFactor=1.5, num.threads = ncore)
      res1 <- as.data.frame(res1)
      res <- rbind(res, res1)      
    }
    res_mtry <- res %>% group_by(mtry) %>% summarize(error=mean(OOBError))
    res_mtry <- res_mtry[order(res_mtry$error),]
    bestmtry <- res_mtry$mtry[1]
    print(paste0("bestmtry: ",bestmtry))
    # cross validation
    set.seed(seeds[4])
    ptid_fold <- createMultiFolds(trainpt$PointID, k = 10, times = 5)

    doranger_f <- function(f){
      td2 <- td1[td1$PointID %in% ptid_fold[[f]],]
      pd2 <- td1[!(td1$PointID %in% ptid_fold[[f]]),]
      set.seed(seeds[5]+f)
      td2_b <- as.data.table(balance(td2, "max", cat_col="Class"))
      set.seed(seeds[6]+f)
      m <- ranger(x = td2_b[,..predictors1], y = td2_b[["Class"]], mtry = bestmtry,
                  seed = seeds[6]+f, num.threads = ncore, verbose = TRUE)
      set.seed(seeds[7]+f)
      m_pred <- predict(m, data = pd2, predict.all = TRUE, 
                        seed = seeds[7]+f, num.threads = ncore, verbose = TRUE)
      m_pred <- rawtoprob(m_pred)
      m_pred <- data.table(pd2[,c("PointID","yearclass")], 
                           repeats=names(ptid_fold[f]), m_pred)
      return(m_pred)
    }
    cv_pred <- data.table()
    for (f in 1:length(ptid_fold)){
      print(paste0("Start fold: ", f))
      m_pred <- doranger_f(f)
      cv_pred <- rbind(cv_pred, m_pred)
    }
    set.seed(seeds[8])
    td_b <- as.data.table(balance(td1, "max", cat_col="Class"))
    set.seed(seeds[9]+f)
    finalmodel <- ranger(x = td_b[,..predictors1], y = td_b[["Class"]], 
                         mtry = bestmtry, importance = "impurity", 
                         seed = seeds[9]+f, num.threads = ncore, verbose = TRUE)
    return(list(modelname=modelname, mtry=bestmtry,
                cv_pred=cv_pred, finalmodel=finalmodel))
  }

  # Do RF
  out <- list()
  p_yes <- predictors
  if("Slope" %in% p_yes == F){p_yes <- append(p_yes,"Slope")}
  if("Northness" %in% p_yes == F){p_yes <- append(p_yes,"Northness")}
  p_no <- predictors
  if("Slope" %in% p_no == T){p_no <- p_no[!p_no == "Slope"]}
  if("Northness" %in% p_no == T){p_no <- p_no[!p_no == "Northness"]}
  
  if (separatewetdry == TRUE){
    p_wet <- doranger(td_wet, p_yes, modelname="wet_withtopo")
    out <- append(out, p_wet)
    print("Finish 1/4")
    p_dry <- doranger(td_dry, p_yes, modelname="dry_withtopo")
    out <- append(out, p_dry)
    print("Finish 2/4")
    p_wet <- doranger(td_wet, p_no, modelname="wet_notopo")
    out <- append(out, p_wet)
    print("Finish 3/4")
    p_dry <- doranger(td_dry, p_no, modelname="dry_notopo")
    out <- append(out, p_dry)
  } else {
    p_all <- doranger(td_all, p_yes, modelname="all_withtopo")
    out <- append(out, p_all)
    print("Finish 1/2")
    p_all <- doranger(td_all, p_no, modelname="all_notopo")
    out <- append(out, p_all)
  }
  return(out)
}


# Functions to combine votes and compute accuracy
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

doaccuracy <- function(pred_df_vote){
  # accuracy assessment
  # https://blogs.fu-berlin.de/reseda/area-adjusted-accuracies/
  # area adjusted accuracy assessment (as in Olofsson et al. 2014)
  confmat <- table(pred_df_vote$predclass, pred_df_vote$refclass)/5
  nclass <- length(unique(pred_df_vote$refclass))
  mappixels <- table(pred_df_vote[pred_df_vote$PointID<=825,]$predclass)
  maparea <- as.numeric(mappixels/sum(mappixels))*1117
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
  OA <- sum(diag(confmat))/sum(confmat)
  PA <- diag(confmat)/colSums(confmat)
  UA <- diag(confmat)/rowSums(confmat)
  # kappa
  diagonal.counts <- diag(confmat)
  N <- sum(confmat)
  row.marginal.props <- rowSums(confmat)/N
  col.marginal.props <- colSums(confmat)/N
  Po <- sum(diagonal.counts)/N
  Pe <- sum(row.marginal.props*col.marginal.props)
  ka <- (Po - Pe)/(1 - Pe)
  
  pred_df_a <- data.frame(class=c("c1", "c2", "c3", "c4", "c5", "c6"),
                          pa=PA, pa_ci=PA_CI, ua=UA, ua_ci=UA_CI)
  pred_df_a <- rbind(pred_df_a, data.frame(class="oa_ka", pa=OA, pa_ci=OA_CI,
                                           ua=ka, ua_ci=0))
  return(pred_df_a)
}

doaccuracy_sepyear <- function(pred_df_vote){
  yearclass <- unique(pred_df_vote$yearclass)
  pred_df_a <- data.frame()
  for (y in yearclass){
    pred_df_y <- pred_df_vote[pred_df_vote$yearclass==y,]
    pred_df_a1 <- doaccuracy(pred_df_y)
    pred_df_a1$yearclass <- y
    pred_df_a <- rbind(pred_df_a, pred_df_a1)
  }
  return(pred_df_a)
}


# run
cv_pred <- buildmodel(traindata, predictors)
saveRDS(cv_pred, "D:/landsat/train/cv_pred.rds")

pred_df <- bind_rows(cv_pred[c(3,7,11,15)])
pred_df$repeats <- substr(pred_df$repeats, nchar(pred_df$repeats), nchar(pred_df$repeats))

pred_df_vote <- voteresult(pred_df)
vresult <- doaccuracy(pred_df_vote)
vresult_sepyear <- doaccuracy_sepyear(pred_df_vote)
write.csv(rbind(cbind(vresult,yearclass="All"),vresult_sepyear),"vresult.csv")


# Compute confusion matrix
confmat <- table(pred_df_vote$predclass, pred_df_vote$refclass)/5
write.csv(confmat, "confmat.csv")
yearclass <- unique(pred_df_vote$yearclass)
y = yearclass[1]
pred_df_y <- pred_df_vote[pred_df_vote$yearclass==y,]
table(pred_df_y$predclass, pred_df_y$refclass)/5


# Create final model (use randomforest)
buildfinmodel <- function(traindata, predictors, rangermodel){
  td_wet <- traindata[traindata$monthclass == "wet"]
  td_dry <- traindata[traindata$monthclass == "dry"]
  # Define rf function
  dorf <- function(td, predictors, mtry){
    set.seed(2021)
    td_b <- as.data.table(balance(td, "max", cat_col="Class"))
    set.seed(2022)
    grid_rf <- expand.grid(mtry = mtry)
    fit_rf <- trainControl(method = "none", seeds = 2022, classProbs = TRUE)
    set.seed(2023)
    m2 <- caret::train(td_b[,..predictors], td_b[["Class"]], method="rf", type="prob",
                       trControl = fit_rf, tuneGrid = grid_rf, 
                       importance = TRUE, verbose = TRUE)
    return(m2)
  }
  # Do RF
  out <- list()
  p_yes <- predictors
  if("Slope" %in% p_yes == F){p_yes <- append(p_yes,"Slope")}
  if("Northness" %in% p_yes == F){p_yes <- append(p_yes,"Northness")}
  p_no <- predictors
  if("Slope" %in% p_no == T){p_no <- p_no[!p_no == "Slope"]}
  if("Northness" %in% p_no == T){p_no <- p_no[!p_no == "Northness"]}
  
  p_wet <- dorf(td_wet, p_yes, mtry=rangermodel$wet_withtopo$mtry)
  out <- append(out, list(wet_withtopo = p_wet))
  p_dry <- dorf(td_dry, p_yes, mtry=rangermodel$dry_withtopo$mtry)
  out <- append(out, list(dry_withtopo = p_dry))
  p_wet <- dorf(td_wet, p_no, mtry=rangermodel$wet_notopo$mtry)
  out <- append(out, list(wet_notopo = p_wet))
  p_dry <- dorf(td_dry, p_no, mtry=rangermodel$dry_notopo$mtry)
  out <- append(out, list(dry_notopo = p_dry))
  return(out)
}

rangermodel <- readRDS("D:/landsat/train/cv_pred.rds")
rangermodel <- list(wet_withtopo=rangermodel[1:4], dry_withtopo=rangermodel[5:8],
                    wet_notopo=rangermodel[9:12], dry_notopo=rangermodel[13:16])
rf_model <- buildfinmodel(traindata, predictors, 
                          rangermodel=rangermodel)
saveRDS(rf_model, "D:/landsat/train/rf_model.rds")


