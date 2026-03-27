# Load required library
library(readr)
source("functions.R")

# Input information
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-m", "--original_metrics"), default = "data/metrics_manual_scaled_rev.csv", help = "Filepath to metrics file")
parser <- add_option(parser, c("-t", "--original_training"), default = "data/manual_training_scaled_rev.csv", help = "Filepath to training data")
parser <- add_option(parser, c("-s", "--original_testing"), default = "data/manual_testing_scaled_rev.csv", help = "Filepath to testing data")
parser <- add_option(parser, c("-a", "--accepted_models_dir"), default = "accepted_models_manual_latest", help = "Filepath to accepted models folder")
parser <- add_option(parser, c("-r", "--result_folder"), default = "results_tms_manual", help = "Filepath to results folder")
parser <- add_option(parser, c("-x", "--max_dimension"), default = 1, type="integer", help = "Maximum number of axis to evaluate")
args <- parse_args(parser)

# original_metrics <- "data/metrics_manual_scaled_rev.csv"
# original_training <- "data/manual_training_scaled_rev.csv"
# original_testing <- "data/manual_testing_scaled_rev.csv"
# accepted_models_dir <- "accepted_models_manual_latest"
# result_folder <- "results_tms_manual"
# max_dimension <- 10

original_metrics <- args$original_metrics
original_training <- args$original_training
original_testing <- args$original_testing
accepted_models_dir <- args$accepted_models_dir
result_folder <- args$result_folder
max_dimension <- args$max_dimension

# Check if the folder exists
if (!dir.exists(result_folder)) {
  # The folder does not exist, so create it
  dir.create(result_folder)
  cat("Folder created:", result_folder, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = result_folder, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

metricsdf <- read.csv(original_metrics)
trainingdf <- read.csv(original_training)
testingdf <- read.csv(original_testing)

for (dimension in 1:max_dimension) {
  summaryfile <- paste0(accepted_models_dir,"/summary_",dimension,".csv")
  summarydf <- read.csv(summaryfile)
  for (row in 1:nrow(summarydf)) {
    feature_pair <- summarydf$Feature_Pair[row]
    features <- unlist(strsplit(feature_pair, "-"))
    alphas <- c(summarydf[row,"Alpha_1"],summarydf[row,"Alpha_2"])
    betas <- as.numeric(summarydf[row, paste0("Beta_", 1:dimension)])
    tms_training <- (alphas[1] + alphas[2]) / 2.0
    tms_testing <- (alphas[1] + alphas[2]) / 2.0
    tms_metrics <- (alphas[1] + alphas[2]) / 2.0
    for (axis in 1:dimension) { 
      tms_training <- tms_training - betas[axis] * as.numeric(trainingdf[,features[axis]])
      tms_testing <- tms_testing - betas[axis] * as.numeric(testingdf[,features[axis]])
      tms_metrics <- tms_metrics - betas[axis] * as.numeric(metricsdf[,features[axis]])
    }
    tms_trainingdf <- trainingdf
    tms_trainingdf$tms <- tms_training
    tms_trainingdf <- tms_trainingdf[, names(tms_trainingdf) %in% c("tms","drug_name","Sample_ID","label")]
    # write.csv(tms_trainingdf, paste0(result_folder,"/tms_training_",dimension,"_",row,".csv") ,row.names=F)
        
    tms_testingdf <- testingdf
    tms_testingdf$tms <- tms_testing
    tms_testingdf <- tms_testingdf[, names(tms_testingdf) %in% c("tms","drug_name","Sample_ID","label")]
    # write.csv(tms_testingdf, paste0(result_folder,"/tms_testing_",dimension,"_",row,".csv") ,row.names=F)
    
    tms_df <- rbind(tms_trainingdf,tms_testingdf)
    write.csv(tms_df, paste0(result_folder,"/tms_",dimension,"_",row,".csv") ,row.names=F)
    
    tms_metricsdf <- metricsdf
    tms_metricsdf$tms <- tms_metrics
    tms_metricsdf <- tms_metricsdf[, names(tms_metricsdf) %in% c("dose","beat","tms","max_dv","EADs","noDepols","noRepols","drug","sample")]
    write.csv(tms_metricsdf, paste0(result_folder,"/tms_metrics_",dimension,"_",row,".csv") ,row.names=F)
  }
}
  
  
  