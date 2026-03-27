library(ggplot2)
library(readr)
library(rms)
require(ROCR)
library(ggplot2)
library(gtools)
library(comprehenr)
library(MASS)

library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input_dir"), default = "results_performance_evals", help = "Filepath to results from performance evaluation stage")
parser <- add_option(parser, c("-o", "--output_dir"), default = "results_calibration_drugs", help = "Filepath to output folder")
args <- parse_args(parser)

# input_dir <- "results_performance_evals_manual/tms_metrics_1_1.csv/tms_metrics_1_1.csv"
# output_dir <- "results_get_calibration_drugs/tms_metrics_1_1.csv"

input_dir <- args$input_dir
output_dir <- args$output_dir

sprintf("input_dir = %s",input_dir)
sprintf("output_dir = %s", output_dir)

# Check if the output_dir folder exists
if (!dir.exists(output_dir)) {
  # The folder does not exist, so create it
  dir.create(output_dir)
  cat("Folder created:", output_dir, "\n")
} else {
  # The folder already exists
  # List all files in the folder
  files <- list.files(path = output_dir, full.names = TRUE)
  # Remove all files in the folder
  if (length(files) > 0) {
    file.remove(files)
    cat("All files in the folder have been removed.\n")
  } else {
    cat("The folder is already empty.\n")
  }
}

results <- data.frame()
for (i in 6:12) {
  filename <- paste0(input_dir,"_",i,"/summary_tms_",i,"_metrics.csv")
  if (file.exists(filename)) {
    print("File exists")
    data <- read.csv(filename)
    data <- na.omit(data)
    results <- rbind(results,data)
  } else {
    print("File doesnot exist")
    next
  }
}

# Minmax normalization 
# For threshold changes
results$Th1_changes_new <- (results$Th1_changes - min(results$Th1_changes))/(max(results$Th1_changes) - min(results$Th1_changes))
results$Th2_changes_new <- (results$Th2_changes - min(results$Th2_changes))/(max(results$Th2_changes) - min(results$Th2_changes))
# For loglik values
results$loglike_rank <- rank(-results$Normalized_logLik, ties.method = "first")
results$loglike_rank_new <- (results$loglike_rank - min(results$loglike_rank))/(max(results$loglike_rank) - min(results$loglike_rank))
results$distance <- sqrt(results$Th1_changes_new^2 + results$Th2_changes_new^2 + results$loglike_rank_new^2)

# Final calibration drugs with lowest distance
final_results <- subset(results,results$distance == min(results$distance))
write.csv(final_results, paste0(output_dir,"/calibration_drugs.csv"), row.names = F)

