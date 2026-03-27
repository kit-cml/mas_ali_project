library(ggplot2)
library(readr)
library(rms)
require(ROCR)
library(ggplot2)
library(gtools)
library(comprehenr)
library(MASS)
source("functions.R")
# Input information======================================================
#--- specify command line arguments
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-a", "--accepted_model_folder"), default = "accepted_models_manual_latest", help = "Filepath to drug candidate list")
parser <- add_option(parser, c("-d", "--tms_dimension"), default = 1, type = "integer", help = "TMS dimension")
parser <- add_option(parser, c("-i", "--tms_pairid"), default = 1, type="integer", help = "Biomarker Pair ID")
parser <- add_option(parser, c("-t", "--mdl_tms_file"), default = "results_tms_manual/tms_1_1.csv", help = "Filepath to TMS file from model development lab")
parser <- add_option(parser, c("-s", "--mdl_sens_file"), default = "results_sensitivity_manual/tms_metrics_1_1.csv/manual_28_sens.csv", help = "Filepath to metrics file of the TMS from model development lab")
parser <- add_option(parser, c("-c", "--calibration_drugs_file"), default = "results_get_calibration_drugs/tms_metrics_1_1.csv/calibration_drugs.csv", help = "Filepath calibration drug file")
parser <- add_option(parser, c("-r", "--mil_training_file"), default = "data/nanion_training_scaled.csv", help = "Filepath to the training data of model implementation lab")
parser <- add_option(parser, c("-v", "--mil_testing_file"), default = "data/nanion_testing_scaled.csv", help = "Filepath to the testing data of model implementation lab")
parser <- add_option(parser, c("-m", "--mil_metrics_file"), default = "data/metrics_nanion_scaled.csv", help = "Filepath to the metrics data of model implementation lab")
parser <- add_option(parser, c("-u", "--results_folder"), default = "results_validate_calibration_drugs", help = "Filepath to the results folder")
args <- parse_args(parser)

accepted_model_folder <- args$accepted_model_folder
tms_dimension <- args$tms_dimension
tms_pairid <- args$tms_pairid
mdl_tms_file <- args$mdl_tms_file
mdl_sens_file <- args$mdl_sens_file
calibration_drugs_file <- args$calibration_drugs_file
mil_training_file <- args$mil_training_file
mil_testing_file <- args$mil_testing_file
mil_metrics_file <- args$mil_metrics_file
results_folder <- args$results_folder

# accepted_model_folder <- "accepted_models_manual_latest"
# tms_dimension <- 1
# tms_pairid <- 2
# mdl_tms_file <- "results_tms_manual/tms_1_2.csv"
# mdl_sens_file <- "results_sensitivity_manual/tms_metrics_1_2.csv/manual_28_sens.csv"
# calibration_drugs_file <- "results_get_calibration_drugs/tms_metrics_1_2.csv/calibration_drugs.csv"
# mil_training_file <- "data/nanion_training_scaled.csv"
# mil_testing_file <- "data/nanion_testing_scaled.csv"
# mil_metrics_file <- "data/metrics_nanion_scaled.csv"
# results_folder <- "results_validate_calibration_drugs"

figure_title <- "28 drugs"
tms_name <- "tms"
tms_unit <- ""
# Check if the folder exists
if (!dir.exists(results_folder)) {
  # The folder does not exist, so create it
  dir.create(results_folder)
  cat("Folder created:", results_folder, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = results_folder, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}
# Get the Threshold 1 and 2 for 28 drugs from model development lab
tms_df <- read.csv(mdl_tms_file)
print(dim(tms_df))  # Check if data is loaded correctly
print(head(tms_df)) # Print first few rows for debugging
tms_df$risk <- tms_df$label
levels <- c("low", "intermediate", "high")
values <- c(1, 2, 3)
tms_df$label <- factor(tms_df$label, levels = levels, labels = values, ordered = TRUE)
sensdf <- read.csv(mdl_sens_file)
th1 <- unique(sensdf["V4"][sensdf["V5"]=="threshold1"])
th2 <- unique(sensdf["V4"][sensdf["V5"]=="threshold2"])
print("======================================================")
print("TMS Thresholds for 28 drugs from model development lab")
print("======================================================")
print(paste0("Threshold 1: ",th1))
print(paste0("Threshold 2: ",th2))

# Get the Threshold 1 and 2 for calibration drugs from model development lab
calibrationdf <- read.csv(calibration_drugs_file)
calibration_drugs <- unlist(strsplit(calibrationdf$drug_pair, "-"))
th1_new <- calibrationdf$Threshold1
th2_new <- calibrationdf$Threshold2
print("===============================================================")
print("TMS Thresholds for calibration drugs from model development lab")
print("===============================================================")
print(paste0("Threshold 1: ",th1_new))
print(paste0("Threshold 2: ",th2_new))

figure_title <- "Development 28 drugs"
# Plot the TMS distribution of 28 drugs from model development lab
drug_colors <- c("azimilide" = "red", "bepridil" = "red", "disopyramide" = "red", "dofetilide" = "red",
                 "ibutilide" = "red", "quinidine" = "red", "sotalol" = "red", "vandetanib" = "red",
                 "astemizole" = "blue", "chlorpromazine" = "blue", "cisapride" = "blue", "clarithromycin" = "blue",
                 "clozapine" = "blue", "domperidone" = "blue", "droperidol" = "blue", "ondansetron" = "blue",
                 "pimozide" = "blue", "risperidone" = "blue", "terfenadine" = "blue", "diltiazem" = "green",
                 "loratadine" = "green", "metoprolol" = "green", "mexiletine" = "green", "nifedipine" = "green",
                 "nitrendipine" = "green", "ranolazine" = "green", "tamoxifen" = "green", "verapamil" = "green")
# drug_colors_calibration <- drug_colors[calibration_drugs]
data <- tms_df
data_wo_outliers <- remove_outliers_grouped(data, group_col = "drug_name", column_name = tms_name)
tms_range <- range(data_wo_outliers[[tms_name]], na.rm = TRUE)
if (is.na(th1_new) & is.na(th2_new)) {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/28.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
} else {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              th1_new = th1_new,
              th2_new = th2_new,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/28.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
}
th1 <- th1_new
th2 <- th2_new
th1_new <- NA
th2_new <- NA
tms_calibration_df <- tms_df[tms_df$drug_name %in% calibration_drugs,]
data <- tms_calibration_df
data_wo_outliers <- remove_outliers_grouped(data, group_col = "drug_name", column_name = tms_name)
tms_range <- range(data_wo_outliers[[tms_name]], na.rm = TRUE)
figure_title <- "Development calibration drugs"
if (is.na(th1_new) & is.na(th2_new)) {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/development.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
} else {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              th1_new = th1_new,
              th2_new = th2_new,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/development.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
}
# Get the Threshold 1 and 2 for calibration drugs from model implementation lab
mil_trainingdf <- read.csv(mil_training_file)
mil_testingdf <- read.csv(mil_testing_file)
mildf28 <- rbind(mil_trainingdf,mil_testingdf)
mil_metricsdf <- read.csv(mil_metrics_file)
accepteddf <- read.csv(paste0(accepted_model_folder,"/summary_",tms_dimension,".csv"))
# accepteddf <- accepteddf[tms_pairid,]
feature_pair <- accepteddf$Feature_Pair[tms_pairid]
features <- unlist(strsplit(feature_pair, "-"))
alphas <- c(accepteddf[tms_pairid,"Alpha_1"],accepteddf[tms_pairid,"Alpha_2"])
betas <- as.numeric(accepteddf[tms_pairid, paste0("Beta_", 1:tms_dimension)])
tms_mil <- (alphas[1] + alphas[2]) / 2.0
tms_metrics_mil <- (alphas[1] + alphas[2]) / 2.0
for (axis in 1:tms_dimension) { 
  tms_mil <- tms_mil - betas[axis] * as.numeric(mildf28[,features[axis]])
  tms_metrics_mil <- tms_metrics_mil - betas[axis] * as.numeric(mil_metricsdf[,features[axis]])
}
mildf28$tms <- tms_mil
mildf28$risk <- mildf28$label
mildf28 <- mildf28[, names(mildf28) %in% c("tms","drug_name","Sample_ID","label","risk")]
mildf28$label <- factor(mildf28$label, levels = levels, labels = values, ordered = TRUE)
mildf <- mildf28[mildf28$drug_name %in% calibration_drugs,]
mil_metricsdf$tms <- tms_metrics_mil
mil_metricsdf <- mil_metricsdf[, names(mil_metricsdf) %in% c("dose","beat","tms","max_dv","EADs","noDepols","noRepols","drug","sample")]
dd <- datadist(mildf)
options(datadist = "dd")
lmod <- try({ lrm(label ~ tms, data = mildf, penalty = 0, x = TRUE, y = TRUE, tol = 1e-10, maxit = 1e6) })
cf <- coefficients(lmod)
cfvec <- c()
for (kint in 1:(length(cf) - 1))
  cfvec[[paste0("intercept", kint)]] <- cf[[kint]]
cfvec[["slope"]] <- cf[[length(cf)]]
ebeta1<-exp(-cfvec[["intercept1"]]); ebeta2<- exp(-cfvec[["intercept2"]])
th1_mil <- log(ebeta1*ebeta2/-(2*ebeta1-ebeta2))/cfvec[["slope"]]
th2_mil <- log(exp(-cfvec[["intercept2"]])-2*exp(-cfvec[["intercept1"]]))/cfvec[["slope"]]
print("==================================================================")
print("TMS Thresholds for calibration drugs from model implementation lab")
print("==================================================================")
print(paste0("Threshold 1: ",th1_mil))
print(paste0("Threshold 2: ",th2_mil))
th1 <- th1_mil
th2 <- th2_mil
th1_new <- NA
th2_new <- NA
data <- mildf
data_wo_outliers <- remove_outliers_grouped(data, group_col = "drug_name", column_name = tms_name)
tms_range <- range(data_wo_outliers[[tms_name]], na.rm = TRUE)
figure_title <- "Implementation calibration drugs"
if (is.na(th1_new) & is.na(th2_new)) {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/implementation.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
} else {
  tmsplotfun2(data = data,
              th1 = th1,
              th2 = th2,
              th1_new = th1_new,
              th2_new = th2_new,
              drug_colors = drug_colors,
              title = figure_title,
              file_name = paste0(results_folder, "/implementation.jpg"),
              tms_name = tms_name,
              tms_range = tms_range)
}
output1df <- data.frame(
  drug = calibration_drugs
)
for (row in 1:length(calibration_drugs)) {
  output1df$risk[row] <- as.integer(unique(mildf$label[mildf$drug_name==calibration_drugs[row]]))-1
}
write.csv(output1df,paste0(results_folder,"/calibration_drugs.csv"),row.names = FALSE)
write.csv(mil_metricsdf,paste0(results_folder,"/mil_tms_metrics.csv"),row.names = FALSE)
thmildf <- data.frame(
  th1 = th1_mil,
  th2 = th2_mil
)
write.csv(thmildf,paste0(results_folder,"/mil_thresholds.csv"),row.names = FALSE)