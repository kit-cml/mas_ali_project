library(readr)

#--- specify command line arguments
library(optparse)
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input_folder"), default = "results_testing_manual", help = "Filepath to the folder of validation results")
parser <- add_option(parser, c("-r", "--results_filename"), default = "selected_model_manual.csv", help = "Filepath to result file")
parser <- add_option(parser, c("-x", "--max_input"), default = 2, type="integer", help = "Maximum number of dimension to evaluate")
args <- parse_args(parser)

input_folder <- args$input_folder
results_filename <- args$results_filename
max_input <- args$max_input

# input_folder <- "results_testing_nanion"
# results_filename <- "selected_model_nanion.csv"
# max_input <- 11

results <- data.frame()
for (i in 1:max_input) {
  filepath <- paste0(input_folder,"/summary_",i,".csv")
  data <- read_csv(filepath, show_col_types = FALSE)
  data$input <- i
  data <- na.omit(data)
  for (j in 1:i) {
    colname <- paste0("Beta_",j)
    data[colname] <- NULL
    colname <- paste0("Mean_",j)
    data[colname] <- NULL
    colname <- paste0("SD_",j)
    data[colname] <- NULL
  }
  results <- rbind(results,data)
}
round(results$Rank_score, 3)
max_validation_score <- max(results$Rank_score)
results <- subset(results, results$Rank_score == max_validation_score)
min_input <- min(results$input)
results <- subset(results, results$input == min_input)
max_logLik <- max(results$logLik)
results <- subset(results, results$logLik == max_logLik)
write.csv(results,results_filename,row.names = FALSE)
