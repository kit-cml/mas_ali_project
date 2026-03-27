# File:         compute_TdP_error.R
# Author:       Kelly Chang
#               Zhihua Li
# Modifier:     ALI X AYD x LTF
# Date:         Nov 2017 (Nov 2023)
# Version:      1.0 (1.1)
# 
# Description:  R script to perform Torsade de Pointes (TdP) risk
#               classification using ordinal logistic regression and
#               leave-one-out cross validation (LOOCV).
#               For help, run this script with command line option "-h".
#

set.seed(100)

#--- load libraries
library(rms)
require(ROCR)
library(ggplot2)
library(gtools)
library(comprehenr)
library(MASS)
source("functions.R")

#--- specify command line arguments
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-d", "--drug_candidate_file"), default = "results_sensitivity/drug_candidates.csv", help = "Filepath to drug candidate list")
parser <- add_option(parser, c("-m", "--metrics_file"), default = "data/metrics.csv", help = "Filepath to metrics file")
parser <- add_option(parser, c("-z", "--drug_size"), default = 6, type="integer", help = "Number of drugs evaluated")
parser <- add_option(parser, c("-s", "--simdir"), default = "results_performance_evals", help = "Filepath to directory for storing results")
parser <- add_option(parser, c("-e", "--evaluation_mode"), default = 1, type="integer", help = "Evaluation mode printed in summary files for performance measures: 1 = worst cases; 2 = median cases")
parser <- add_option(parser, c("-t", "--sensitivity_file"), default = NULL, help = "Filepath to sensitivity file")
parser <- add_option(parser, c("-l", "--th1"), default = "NULL", help = "The referenced Threshold 1")
parser <- add_option(parser, c("-r", "--th2"), default = "NULL", help = "The referenced Threshold 2")
args <- parse_args(parser)

# Input information
# drug_candidate_file <- "results_sensitivity_nanion/tms_metrics_2_19.csv/nanion_drug_candidates.csv"
# metrics_file <- "results_tms_nanion/tms_metrics_2_19.csv"
# drug_size = 6 #Adjust combine_drug value based on number of drugs
# simdir <- "results_performance_evals_nanion"
# eval_mode <- 1
# sensitivity_file <- "results_sensitivity_nanion/tms_metrics_2_19.csv/nanion_28_sens.csv"
# th1 <- NULL
# th2 <- NULL
# print("check!")
drug_candidate_file <- args$drug_candidate_file
metrics_file <- args$metrics_file
drug_size <- as.integer(args$drug_size)
simdir <- args$simdir
eval_mode <- as.integer(args$evaluation_mode)
sensitivity_file <- args$sensitivity_file

# sensitivity_file <- as.character(args$sensitivity_file)

sprintf("drug_candidate_file = %s",drug_candidate_file)
sprintf("metrics_file = %s", metrics_file)
sprintf("drug_size = %s",drug_size)
sprintf("simdir = %s",simdir)
sprintf("eval_mode = %s",eval_mode)
if (is.null(sensitivity_file)) {
  sprintf("sensitivity_file = NULL")
  th1 <- as.numeric(args$th1)
  th2 <- as.numeric(args$th2)
  sprintf("Th1 = %s",as.character(th1))
  sprintf("Th2 = %s",as.character(th2))
} else {
  sprintf("sensitivity_file = %s", sensitivity_file)
}
check <- (eval_mode == 1 | eval_mode == 2)
if (!check) {
  stop("Wrong evaluation mode entered. The possible options are 1 = worst case evaluation and 2 = median case evaluation")
}

#--- get arguments
outputI <- TRUE
metric <- "tms"
metricv <- c("tms") 
scale <- "sara"
CL <- 2000
datasetvec <- c("chantest")  

# Read drug candidates data
drug_candidates <- as.data.frame(read.csv(drug_candidate_file)) 
if (is.na(drug_size)) {
  drug_size <- nrow(drug_candidates)
}

# setup output directory
input <- "uncertainty"
simdir <- paste0(simdir,"_",drug_size)

# Check if the folder exists
if (!dir.exists(simdir)) {
  # The folder does not exist, so create it
  dir.create(simdir)
  cat("Folder created:", simdir, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = simdir, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

#--- creating drug combination
drug_combination <- combine_drug(drug_size,drug_candidates)

# Obtain the Th1 and Th2 from drug candidates
if (!is.null(sensitivity_file)) {
  sensdf <- read.csv(sensitivity_file)
  Th1_ref <- unique(sensdf["V4"][sensdf["V5"]=="threshold1"])
  Th2_ref <- unique(sensdf["V4"][sensdf["V5"]=="threshold2"])
} else {
  Th1_ref <- th1
  Th2_ref <- th2
}
measurevec <- c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
                "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
                "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
                "Threshold1", "Threshold2", "Th1_ref", "Th2_ref", "Th1_changes", "Th2_changes", "Th_changes_overall", "Normalized_logLik")


columns <- c(to_vec(for(i in 1:drug_size) paste0("drug",i)), measurevec) 
result_df <- data.frame(matrix(nrow = 0, ncol = length(columns)))
colnames(result_df) <- columns
table7_all <- data.frame()

for(drug_index in seq_len(nrow(drug_combination))) {
# for(drug_index in 1:2) {
  print(paste0("-------------------------------------", drug_index)) 
  drug_pair <- drug_combination[drug_index,]
  drugtable <- data.frame(drug = c(drug_pair))
  for (row in 1:nrow(drugtable)) {
    drugtable$CiPA[row] = drug_candidates$risk[drug_candidates$drug==drugtable$drug[row]]
  }
  drugnames <- as.character(drugtable$drug)
  
  # read in dataset
  infile <- metrics_file
  df <- read.csv(infile)
  df <- df[df$drug != "control" & df$dose != 0, ]
  df <- df[!is.na(df$max_dv), ] # depolarization failures
  df <- df[(df$dose) %in% c("1", "2", "3", "4"), ]                  #very important cause some drugs may have additional doses in metrics.rds!!!!
  if (input == "uncertainty") {
    qNettable <- aggregate(df[, metric], list(df$drug, df$sample), mean)  #name is qNettable, but it could be any metric
    colnames(qNettable) <- c("drug", "sample", metric)
  } else {
    qNettable <- aggregate(df[, metric], list(df$drug), mean)
    colnames(qNettable) <- c("drug", "qNet")
    qNettable$sample <- NA
  }
  df <- merge(qNettable, drugtable[, c("drug", "CiPA")], by.x = "drug", by.y = "drug")

  #organize df
  df$drug <- factor(df$drug, levels = drugnames)
  df <- df[order(df$drug, df$sample), ]
  df$CiPA <- ordered(df$CiPA)
  
  # perform logistic regression and leave-one-out cross validation for each drug
  errdf <- data.frame()
  probdf <- data.frame()
  allprobdf <- data.frame()   #corresponding to cvdf but without LOOCV
  cverrdf <- data.frame()
  cvprobdf <- data.frame()
  cvdf <- data.frame()                              #used to store individual sample's probs
  is_training_all_error <- FALSE
  
  for (drug in c(NA, drugnames)) {
    
    # fit logistic regression model
    if (is.na(drug)) {
      print("training on all drugs")
      datadf <- df[, c("drug", "CiPA", "sample", metric)]
      testdf <- datadf
      traindf <- datadf
    } else {
      print(sprintf("cross validating with %s", drug))
      testdf <- datadf[datadf$drug == drug, ]
      traindf <- datadf[datadf$drug != drug, ]
      if (nrow(testdf) == 0) {
        print("no data to test! skipping...")
        next
      }
    }
    
    # lmod <- try_lrm_loocv(traindf)
    dd <- datadist(traindf)
    options(datadist = "dd")
    lmod <- try({ lrm(CiPA ~ get(metric), data = traindf, penalty = 0, x = TRUE, y = TRUE, tol = 1e-10, maxit = 1e6) })
    # print(lmod)
    
    if (inherits(lmod, "try-error")) {
      if (is.na(drug)) {
        print("fitting all drugs failed? skipping...")
        is_training_all_error <- TRUE
        break
      } else {
        print(sprintf("cross validating with %s failed! skipping...", drug))
        next
      }
    }
    
    # save coefficients
    print(sprintf("Convergence failure: %s", lmod$fail))
    cf <- coefficients(lmod)
    cfvec <- c()
    for (kint in 1:(length(cf) - 1))
      cfvec[[paste0("intercept", kint)]] <- cf[[kint]]
    cfvec[["slope"]] <- cf[[length(cf)]]
    
    #use formal math for threshold 1 and 2
    ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
    t1 <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
    t2 <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]
    
    # get training/prediction error
    y0 <- testdf$CiPA
    if (is.na(drug)) {
      probs <- predict(lmod, type = "fitted.ind")
    } else {
      probs <- matrix(predict(lmod, newdata = testdf, type = "fitted.ind"), ncol = length(levels(y0)))
    }
    
    #sometimes APD90 is NA, giving rise to NA in probs
    idx <- apply(probs, 1, function(x) any(is.na(x)))
    probs[idx, ] <- matrix(rep(c(0, 0, 1), sum(idx)), nrow = sum(idx), byrow = T) #this is to assume all NAs are due to repolarization failure, and thus high risk
    
    yPred <- apply(probs, 1, function(x) which.max(x))
    pred_err <- mean(abs(yPred - as.integer(y0)))
    
    if (is.na(drug)) {
      # store the value of thresholds
      threshold1 <- t1
      threshold2 <- t2
      
      # Loglikelihood value of the model
      logLikVal <- numeric(length(y0))
      
      for (i in 1:length(y0)) {
        logLikVal[i] <- log(probs[i,y0[i]])
      }
      
      # Normalized loglikelihood value
      nrows_trainingdf <- nrow(traindf)
      normalized_logLik <- logLikVal / nrows_trainingdf
    }
    
    # append to data frame
    if (is.na(drug)) {
      newrow <- data.frame(t(cfvec), error = pred_err)
      errdf <- rbind(errdf, newrow)
      
      # allprobdf <- cbind(traindf[, 1:3], probs)
      # colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob")
      allprobdf <- cbind(traindf[, 1:3], probs, abs(yPred - as.integer(y0)), normalized_logLik)
      colnames(allprobdf) <- c("drug", "CiPA", "sample", "low_prob", "inter_prob", "high_prob", "pred_err", "normalized_logLik")
      
    } else {
      newrow <- data.frame(drug = drug, t(cfvec), error = pred_err)
      cverrdf <- rbind(cverrdf, newrow)
      if (outputI) {                            #make detailed prob table here
        testdf <- cbind(testdf, high_prob = probs[, 3], inter_prob = probs[, 2], low_prob = probs[, 1], pred_err = abs(yPred - as.integer(y0)))
        cvdf <- rbind(cvdf, testdf)
      } #if outputI
    }
    
    # detailed results
    newcols <- sapply(levels(y0), function(s) paste0("predict_", s), USE.NAMES = F)
    predictions <- factor(newcols[yPred], levels = newcols)
    tmpdf <- as.data.frame(tapply(1:nrow(testdf), list(drug = as.character(testdf$drug), predict = predictions), FUN = length))
    tmpdf[is.na(tmpdf)] <- 0
    tmpdf$drug <- factor(rownames(tmpdf), levels = drugnames)
    
    tmpdf <- merge(drugtable[, c("drug", "CiPA")], tmpdf)
    tmpdf <- tmpdf[order(tmpdf$drug), ]
    rownames(tmpdf) <- NULL
    if (is.na(drug)) {
      probdf <- rbind(probdf, tmpdf)
      
    } else {
      cvprobdf <- rbind(cvprobdf, tmpdf)
    }
    
  } # for drug
  
  if (!is_training_all_error) {
    outfile <- paste0(simdir, "/",metric, "_", drug_index, "_training_probs.csv")   #here training means "trained on all drugs"
    # write.csv(probdf, outfile, row.names = F, quote = F)
    # print(head(probdf))
    
    outfile <- paste0(simdir, "/",metric,  "_", drug_index, "_training_errors.csv")
    list_columns <- sapply(errdf, is.list) # Identify columns in errdf that are of type list
    if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
      errdf[list_columns] <- lapply(errdf[list_columns], as.numeric)
    }
    # write.csv(errdf, outfile, row.names = F, quote = F)
    # print(head(errdf))
    
    outfile <- paste0(simdir, "/",metric,  "_", drug_index, "_LOOCV_probs.csv")
    # write.csv(cvprobdf, outfile, row.names = F, quote = F)
    # print(head(cvprobdf))
    
    outfile <- paste0(simdir, "/",metric,  "_", drug_index, "_LOOCV_errors.csv")
    list_columns <- sapply(cverrdf, is.list) # Identify columns in errdf that are of type list
    if (any(list_columns)) { # Convert list columns to character (or other suitable type) if any
      cverrdf[list_columns] <- lapply(cverrdf[list_columns], as.numeric)
    }
    # write.csv(cverrdf, outfile, row.names = F, quote = F)
    # print(head(cverrdf))
    
    if (outputI) {
      outfile <- paste0(simdir,"/", metric,  "_", drug_index, "_LOOCV_allprobes.csv")
      write.csv(cvdf, outfile, row.names = F, quote = F)
      # print(head(cvdf))
      outfile <- paste0(simdir,"/", metric,  "_", drug_index, "_training_allprobes.csv")
      write.csv(allprobdf, outfile, row.names = F, quote = F)
    }
  }

  if (!is_training_all_error) {
    
    # Initialize data frames
    # all28df <- data.frame(metric=c("AUC1","AUC2","Pairwise","LR1plus","LR1minus","LR2plus","LR2minus","Mean_error",
    #                                "AUC1_LOOCV","AUC2_LOOCV","Pairwise_LOOCV","LR1plus_LOOCV","LR1minus_LOOCV","LR2plus_LOOCV","LR2minus_LOOCV","Mean_error_LOOCV",
    #                                "AUC1_training","AUC2_training","Pairwise_training","LR1plus_training","LR1minus_training","LR2plus_training","LR2minus_training","Mean_error_training",
    #                                "Threshold1", "Threshold2", "Normalized_logLik"))
    all28df <- data.frame(metric = measurevec )
    colordf <- data.frame(
      drug = c("dofetilide","bepridil","quinidine","sotalol","ibutilide", "azimilide", "disopyramide", "vandetanib",
               "cisapride","terfenadine","ondansetron","chlorpromazine","clarithromycin","risperidone","domperidone","astemizole", "pimozide","droperidol","clozapine",
               "ranolazine","mexiletine","diltiazem","verapamil","metoprolol","loratadine","tamoxifen","nifedipine","nitrendipine"),
      classidx=c(rep(2,8),rep(1,11),rep(0,9)),
      coloridx=c(rep(1,8),rep(3,11),rep(2,9)),
      isTraining=c(1,1,1,1,0,0,0,0,
                   1,1,1,1,0,0,0,0,0,0,0,
                   1,1,1,1,0,0,0,0,0)
    )
    
    # Loop over datasets and metrics
    for(dataset in datasetvec){
      for(metric in metricv){
        altstr <- paste0(dataset,"_",metric)
        altvec <- 0
        valdf <- read.delim(paste0(simdir,"/",metric,"_", drug_index, "_LOOCV_allprobes.csv"),sep=",",as.is=T)
        trainingdf <- read.delim(paste0(simdir,"/",metric,"_", drug_index, "_training_allprobes.csv"),sep=",",as.is=T)
        file.remove(paste0(simdir,"/",metric,"_", drug_index, "_LOOCV_allprobes.csv"))
        file.remove(paste0(simdir,"/",metric,"_", drug_index, "_training_allprobes.csv"))
        
        # Classify data
        valdf$classidx <- valdf$CiPA
        valdf$class <- 2 - valdf$classidx 
        
        trainingdf$classidx <- trainingdf$CiPA
        trainingdf$class <- 2 - trainingdf$classidx 
        
        # Prepare data for ROC analysis
        forROC1 <- valdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1$inter_prob <- rowSums(valdf[,c("inter_prob","high_prob")])
        colnames(forROC1)[1] <- "positive_prob"
        forROC1$class <- forROC1$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        forROC1_training <- trainingdf[,c("inter_prob","classidx","sample", "drug")]
        forROC1_training$inter_prob <- rowSums(trainingdf[,c("inter_prob","high_prob")])
        colnames(forROC1_training)[1] <- "positive_prob"
        forROC1_training$class <- forROC1_training$classidx > 0  # now high and intermediate are both "TRUE" or "positive" 
        
        # Generate predictions and labels
        forROC1predictions <- do.call(cbind,by(forROC1[,"positive_prob"],forROC1$sample,function(x) x))
        forROC1labels <- do.call(cbind,by(forROC1$class,forROC1$sample,function(x) x))
        
        forROC1predictions_training <- do.call(cbind,by(forROC1_training[,"positive_prob"],forROC1_training$sample,function(x) x))
        forROC1labels_training <- do.call(cbind,by(forROC1_training$class,forROC1_training$sample,function(x) x))
        
        # Perform ROC analysis
        ROC1obj <- prediction(forROC1predictions, forROC1labels)        
        AUC1 <- performance(ROC1obj,"auc")  #the value is in y.values
        AUC1dist <- formatC(quantile(unlist(AUC1@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC1dist[2]," (",AUC1dist[1]," - ",AUC1dist[3],")"))
        
        ROC1obj_training <- prediction(forROC1predictions_training, forROC1labels_training)        
        AUC1_training <- performance(ROC1obj_training,"auc")  #the value is in y.values
        AUC1dist_training <- formatC(quantile(unlist(AUC1_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC1dist_training[2]," (",AUC1dist_training[1]," - ",AUC1dist_training[3],")"))
        
        # Prepare data for second ROC analysis
        forROC2 <- valdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2)[1] <- "positive_prob"
        forROC2$class <- forROC2$classidx > 1                         
        
        forROC2_training <- trainingdf[,c("high_prob","classidx","sample","drug")]
        colnames(forROC2_training)[1] <- "positive_prob"
        forROC2_training$class <- forROC2_training$classidx > 1
        
        # Generate predictions and labels
        forROC2predictions <- do.call(cbind,by(forROC2[,"positive_prob"],forROC2$sample,function(x) x))  
        forROC2labels <- do.call(cbind,by(forROC2$class,forROC2$sample,function(x) x))
        
        forROC2predictions_training <- do.call(cbind,by(forROC2_training[,"positive_prob"],forROC2_training$sample,function(x) x))  
        forROC2labels_training <- do.call(cbind,by(forROC2_training$class,forROC2_training$sample,function(x) x))
        
        # Perform second ROC analysis
        ROC2obj <- prediction(forROC2predictions, forROC2labels)        
        AUC2 <- performance(ROC2obj,"auc")  #the value is in y.values 
        AUC2dist <- formatC(quantile(unlist(AUC2@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(AUC2dist[2]," (",AUC2dist[1]," - ",AUC2dist[3],")"))
        
        ROC2obj_training <- prediction(forROC2predictions_training, forROC2labels_training)        
        AUC2_training <- performance(ROC2obj_training,"auc")  #the value is in y.values 
        AUC2dist_training <- formatC(quantile(unlist(AUC2_training@y.values),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(AUC2dist_training[2]," (",AUC2dist_training[1]," - ",AUC2dist_training[3],")"))
        
        # Pairwise ranking
        pairwisefun <- function(fulltable){
          cmb <- combn(seq_len(nrow(fulltable)), 2)
          mergedtable <- cbind(fulltable[cmb[1,],], fulltable[cmb[2,],])
          validpairidx <- (mergedtable[,7]!=mergedtable[,16])&(!mergedtable[,9]|!mergedtable[,18])
          correctidx1 <- ((mergedtable[,7]>mergedtable[,16])&(mergedtable[,6]<mergedtable[,15]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,6]>mergedtable[,15])) #when predicted class are different
          correctidx2 <- (mergedtable[,6]==1)&(mergedtable[,15]==1)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,3]>mergedtable[,12]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,3]<mergedtable[,12]))) #when predicted class are both high
          correctidx3 <- (mergedtable[,6]==3)&(mergedtable[,15]==3)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both low
          correctidx4 <- (mergedtable[,6]==2)&(mergedtable[,15]==2)&(((mergedtable[,7]>mergedtable[,16])&(mergedtable[,5]<mergedtable[,14]))|((mergedtable[,7]<mergedtable[,16])&(mergedtable[,5]>mergedtable[,14]))) #when predicted class are both intermediate
          correctidx <- correctidx1|correctidx2|correctidx3|correctidx4
          sum(validpairidx&correctidx)/sum(validpairidx)
        }
        tempdf <- valdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf$pred <- apply(valdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf <- merge(tempdf,colordf,by="drug")
        distpairwise <- by(tempdf, valdf$sample, pairwisefun)
        pairwisedist <- formatC(quantile(unlist(distpairwise),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(pairwisedist[2]," (",pairwisedist[1]," - ",pairwisedist[3],")"))
        
        tempdf_training <- trainingdf[,c("drug","sample","high_prob","inter_prob","low_prob")]; 
        tempdf_training$pred <- apply(trainingdf[,5:7],1,which.max) #note here 1 is high, 2 is int, 3 is low
        tempdf_training <- merge(tempdf_training,colordf,by="drug")
        distpairwise_training <- by(tempdf_training, trainingdf$sample, pairwisefun)
        pairwisedist_training <- formatC(quantile(unlist(distpairwise_training),prob=c(0.025,0.5,0.975),na.rm=TRUE),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(pairwisedist_training[2]," (",pairwisedist_training[1]," - ",pairwisedist_training[3],")"))
        
        # Calculate sensitivity and specificity
        sens1 <- with(valdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1 <-  with(valdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1 <- sens1/(1-spec1)
        LRminus1 <- (1-sens1)/spec1
        
        sens2 <- with(valdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2 <-  with(valdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2 <- sens2/(1-spec2)
        LRminus2 <- (1-sens2)/spec2
        
        sens1_training <- with(trainingdf, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0))
        spec1_training <-  with(trainingdf, sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0))
        LRplus1_training <- sens1_training/(1-spec1_training)
        LRminus1_training <- (1-sens1_training)/spec1_training
        
        sens2_training <- with(trainingdf, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2))
        spec2_training <-  with(trainingdf, sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2))
        LRplus2_training <- sens2_training/(1-spec2_training)
        LRminus2_training <- (1-sens2_training)/spec2_training
        
        # Calculate LR using random sampling
        sens1UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ <- (sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(1-spec1UQ+rnorm(length(sens1UQ),1e-6,1e-12))
        LRplus1dist <- formatC(quantile(unlist(LRplus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus1dist[2]," (",LRplus1dist[1]," - ",LRplus1dist[3],")"))
        LRminus1UQ <- (1-sens1UQ+rnorm(length(sens1UQ),1e-6,1e-12))/(spec1UQ+rnorm(length(sens1UQ),1e-12,1e-12))
        LRminus1dist <- formatC(quantile(unlist(LRminus1UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus1dist[2]," (",LRminus1dist[1]," - ",LRminus1dist[3],")"))
        
        sens2UQ <- by(valdf, valdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ <- by(valdf, valdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ <- (sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(1-spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRplus2dist <- formatC(quantile(unlist(LRplus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRplus2dist[2]," (",LRplus2dist[1]," - ",LRplus2dist[3],")"))
        LRminus2UQ <- (1-sens2UQ+rnorm(length(sens2UQ),1e-6,1e-12))/(spec2UQ+rnorm(length(sens2UQ),1e-6,1e-12))
        LRminus2dist <- formatC(quantile(unlist(LRminus2UQ),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(LRminus2dist[2]," (",LRminus2dist[1]," - ",LRminus2dist[3],")"))
        
        sens1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx>0&(low_prob<high_prob|low_prob<inter_prob))/sum(classidx>0)))
        spec1UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx==0&low_prob>high_prob&low_prob>inter_prob)/sum(classidx==0)))
        LRplus1UQ_training <- (sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(1-spec1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))
        LRplus1dist_training <- formatC(quantile(unlist(LRplus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus1dist_training[2]," (",LRplus1dist_training[1]," - ",LRplus1dist_training[3],")"))
        LRminus1UQ_training <- (1-sens1UQ_training+rnorm(length(sens1UQ_training),1e-6,1e-12))/(spec1UQ_training+rnorm(length(sens1UQ_training),1e-12,1e-12))
        LRminus1dist_training <- formatC(quantile(unlist(LRminus1UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus1dist_training[2]," (",LRminus1dist_training[1]," - ",LRminus1dist_training[3],")"))
        
        sens2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x, sum(classidx==2&low_prob<high_prob&inter_prob<high_prob)/sum(classidx==2)))
        spec2UQ_training <- by(trainingdf, trainingdf$sample, function(x) with(x,sum(classidx<2&(low_prob>high_prob|inter_prob>high_prob))/sum(classidx<2)))
        LRplus2UQ_training <- (sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(1-spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRplus2dist_training <- formatC(quantile(unlist(LRplus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRplus2dist_training[2]," (",LRplus2dist_training[1]," - ",LRplus2dist_training[3],")"))
        LRminus2UQ_training <- (1-sens2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))/(spec2UQ_training+rnorm(length(sens2UQ_training),1e-6,1e-12))
        LRminus2dist_training <- formatC(quantile(unlist(LRminus2UQ_training),prob=c(0.025,0.5,0.975)),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(LRminus2dist_training[2]," (",LRminus2dist_training[1]," - ",LRminus2dist_training[3],")"))
        
        # Calculate mean error
        meanerror <- mean(valdf$pred_err)
        lowererror <- meanerror - 1.96* sd(valdf$pred_err)/sqrt(dim(valdf)[1]); uppererror <- meanerror+1.96*sd(valdf$pred_err)/sqrt(dim(valdf)[1])
        errordist <- formatC(c(meanerror,lowererror,uppererror),digits=3,format="g",flag="-",drop0trailing=T)
        altvec <- c(altvec,paste0(errordist[1]," (",errordist[2]," - ",errordist[3],")"))
        
        meanerror_training <- mean(trainingdf$pred_err)
        lowererror_training <- meanerror_training - 1.96* sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1]); uppererror_training <- meanerror_training+1.96*sd(trainingdf$pred_err)/sqrt(dim(trainingdf)[1])
        errordist_training <- formatC(c(meanerror_training,lowererror_training,uppererror_training),digits=3,format="g",flag="-",drop0trailing=T)
        # altvec <- c(altvec,paste0(errordist_training[1]," (",errordist_training[2]," - ",errordist_training[3],")"))
        
        if (eval_mode == 1) { # Worst case evaluation
          # LOOCV results for ranking
          altvec <- c(altvec,paste0(AUC1dist[1]))
          altvec <- c(altvec,paste0(AUC2dist[1]))
          altvec <- c(altvec,paste0(pairwisedist[1]))
          altvec <- c(altvec,paste0(LRplus1dist[1]))
          altvec <- c(altvec,paste0(LRminus1dist[3]))
          altvec <- c(altvec,paste0(LRplus2dist[1]))
          altvec <- c(altvec,paste0(LRminus2dist[3]))
          altvec <- c(altvec,paste0(errordist[3]))
          
          # Training results for ranking
          altvec <- c(altvec,paste0( AUC1dist_training[1]))
          altvec <- c(altvec,paste0( AUC2dist_training[1]))
          altvec <- c(altvec,paste0(pairwisedist_training[1]))
          altvec <- c(altvec,paste0(LRplus1dist_training[1]))
          altvec <- c(altvec,paste0(LRminus1dist_training[3]))
          altvec <- c(altvec,paste0(LRplus2dist_training[1]))
          altvec <- c(altvec,paste0(LRminus2dist_training[3]))
          altvec <- c(altvec,paste0(errordist_training[3]))
        } else if (eval_mode == 2) { # Median evaluation
          # LOOCV results for ranking
          altvec <- c(altvec,paste0(AUC1dist[2]))
          altvec <- c(altvec,paste0(AUC2dist[2]))
          altvec <- c(altvec,paste0(pairwisedist[2]))
          altvec <- c(altvec,paste0(LRplus1dist[2]))
          altvec <- c(altvec,paste0(LRminus1dist[2]))
          altvec <- c(altvec,paste0(LRplus2dist[2]))
          altvec <- c(altvec,paste0(LRminus2dist[2]))
          altvec <- c(altvec,paste0(errordist[1]))
          
          # Training results for ranking
          altvec <- c(altvec,paste0( AUC1dist_training[2]))
          altvec <- c(altvec,paste0( AUC2dist_training[2]))
          altvec <- c(altvec,paste0(pairwisedist_training[2]))
          altvec <- c(altvec,paste0(LRplus1dist_training[2]))
          altvec <- c(altvec,paste0(LRminus1dist_training[2]))
          altvec <- c(altvec,paste0(LRplus2dist_training[2]))
          altvec <- c(altvec,paste0(LRminus2dist_training[2]))
          altvec <- c(altvec,paste0(errordist_training[1]))
        }
        
        # Thresholds
        altvec <- c(altvec,paste0(threshold1))
        altvec <- c(altvec,paste0(threshold2))
        
        # if (!is.null(sensitivity_file)) {
          # Referenced thresholds
          altvec <- c(altvec,paste0(Th1_ref))
          altvec <- c(altvec,paste0(Th2_ref))
          
          # Threshold 1 & 2 changes
          altvec <- c(altvec,paste0(abs((threshold1-Th1_ref)/Th1_ref)))
          altvec <- c(altvec,paste0(abs((threshold2-Th2_ref)/Th2_ref)))
          
          # Threshold changes overall
          altvec <- c(altvec,paste0(abs((threshold1-Th1_ref)/Th1_ref) * abs((threshold2-Th2_ref)/Th2_ref)))
        # }
        
        # The normalized loglikelihood value
        altvec <- c(altvec,paste0(mean(normalized_logLik)))
      
        # Update altvec and all28df
        altvec <- altvec[-1]
        all28df[,altstr] <- altvec
      }
    }
    
    # Write results to file
    
    table7 <- data.frame(metric=rep(measurevec,each=length(datasetvec)), 
                         dataset=rep(datasetvec, length(measurevec)),
                         qNet=NA)
    table7[,-1:-2] <- matrix(unlist(t(all28df[,-1])),nrow= length(measurevec)*length(datasetvec),byrow=T)
    table7$drug_index <- drug_index
    # outfile <- paste0(simdir,"/",drug_index,"_Table7.txt")
    # write.table(table7, outfile, col.names=T,row.names=F,sep="\t",quote=F)
    table7_all <- rbind(table7_all, table7)
    write.table(table7_all, file = paste0(simdir, "/All_Table7.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    table7$drug_index <- NULL
    table7 <- subset(table7, select=-c(dataset))
    #table7 <- subset(table7, select=-c(1))
    names(table7) <- NULL
    table7 <- setNames(data.frame(t(table7[,-1])), table7[,1])
    drugnames <- t(drugnames)
    colnames(drugnames) <- c(to_vec(for(i in 1:drug_size) paste0("drug_",i)))
    drug_pair_df <- data.frame(
      drug_pair = paste(drug_pair,collapse = "-")
    )
    final_result <- cbind(drug_pair_df, table7[1,])
    
    # Add number of drugs and Rank score LOOCV
    final_result["number_of_drugs"] <- drug_size
    final_result$Rank_score_LOOCV <- rankscorefun_loocv(final_result)
    result_df <- rbind(result_df, final_result)
    outfile <- paste0(simdir, "/summary_",metric,"_",drug_size, "_metrics.csv")
    write.csv(result_df, outfile ,row.names=F)
  } #if is_training_all_error
}

# # Add number of drugs and Rank score LOOCV
# result_df["number_of_drugs"] <- drug_size
# result_df$Rank_score_LOOCV <- NA
# for (j in 1:nrow(result_df)) {
#   result_df[j,]$Rank_score_LOOCV <- rankscorefun_loocv(result_df[j,])
# }
# write.csv(result_df, outfile ,row.names=F)