library(readr)
library(MASS)
library(pROC)
library(ggplot2)
library(foreach)
library(doParallel)
source("functions.R")

# Declare features and units
features <- "tms"
# features <- c("qNet",
#               "dvdtmax",
#               "vmax",
#               "vrest",
#               "APD50",
#               "APD90",
#               "max_dv",
#               "camax",
#               "carest",
#               "CaTD50",
#               "CaTD90")

units <- ""
# units <- c("(nC/uF)",
#            "(mV/ms)",
#            "(mV)",
#           "(mV)",
#           "(ms)",
#           "(ms)",
#           "(mV/ms)",
#           "(mM)",
#           "(mM)",
#           "(ms)",
#            "(ms)")

# Declare the file paths
filepath_training <- "data/chantest_training_scaled_w_tms_vrest-APD90-CaTD90.csv"
filepath_testing <- "data/chantest_testing_scaled_w_tms_vrest-APD90-CaTD90.csv"
# filepath_training <- "data/manual_training_scaled.csv"
# filepath_testing <- "data/manual_testing_scaled.csv"

# Set feature dimension
dimension <- 1

# Create pairsdf with all unique combinations
pairsdf <- pairsdfinitfun(features = features, units = units, dimension = dimension)

# Create the results folder
results_folder <- "results_chantest_w_tms_vrest-APD90-CaTD90_latest"

# Choose whether data needs to be normalized
is_normalized <- TRUE

# Choose how many attemps required before skipping the fitting process (when divergence is found)
max_attempts <- 10000

# Register parallel backend
numCores <- 1

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

# Execute the tasks in parallel
cl <- makeCluster(numCores)
registerDoParallel(cl)

is_loocv <- FALSE
summarydf <- foreach( pair_id = 1:nrow(pairsdf),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
                       # Select the row corresponding to pair_id
                       pair_row <- pairsdf[pair_id, ]
                       
                       # Extract vectors of features and units
                       features_vector <- pair_row[grep("feature_", names(pair_row))]
                       units_vector <- pair_row[grep("unit_", names(pair_row))]
                       
                       result <- solvepair(results_folder = results_folder,
                                           pair_id = pair_id,
                                           filepath_training = filepath_training,
                                           filepath_testing = filepath_testing,
                                           features_vector = features_vector,
                                           units_vector = units_vector,
                                           is_normalized = is_normalized,
                                           max_attempts = max_attempts,
                                           is_loocv = is_loocv)
                       
                       # Return the result
                       result
                     }
write.csv(summarydf, paste(results_folder,"/","summary_training.csv", sep = ""), row.names = FALSE)

is_loocv <- TRUE
# Remove NA from training results
summarydf <- na.omit(summarydf)
summaryloocvdf <- foreach(pair_id = 1:nrow(pairsdf),
                     .combine = 'rbind',
                     .packages = c("readr","MASS","pROC","ggplot2")) %dopar% {
                       # Select the row corresponding to pair_id
                       pair_row_loocv <- pairsdf[pair_id, ]

                       # Extract vectors of features and units
                       features_vector_loocv <- pair_row_loocv[grep("feature_", names(pair_row_loocv))]
                       units_vector_loocv <- pair_row_loocv[grep("unit_", names(pair_row_loocv))]

                       # Check whether the the feature pair is available or not
                       feature_pair <- if (dimension != 1) {
                         paste(features_vector_loocv, collapse="-")
                       } else {
                         features_vector_loocv
                       }
                       is_available <- sum(summarydf[1] == as.character(feature_pair))

                       if (is_available != 0) {
                         resultloocv <- solvepair(results_folder = results_folder,
                                                  pair_id = pair_id,
                                                  filepath_training = filepath_training,
                                                  filepath_testing = filepath_testing,
                                                  features_vector = features_vector_loocv,
                                                  units_vector = units_vector_loocv,
                                                  is_normalized = is_normalized,
                                                  max_attempts <- max_attempts,
                                                  is_loocv = is_loocv)

                         # Return the result
                         resultloocv
                       }
                     }

# Stop the parallel cluster
stopCluster(cl)
# Save the summarydf dataframe
write.csv(summaryloocvdf, paste(results_folder,"/","summary_loocv.csv", sep = ""), row.names = FALSE)