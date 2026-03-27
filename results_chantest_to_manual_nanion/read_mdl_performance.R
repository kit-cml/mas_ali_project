# Set the parent directory
parent_dir <- "results_get_calibration_drugs"

# Get all subdirectories
all_dirs <- list.dirs(parent_dir, recursive = TRUE, full.names = TRUE)

# Find folders that contain 'calibration_drugs.csv'
folders_with_calibration <- all_dirs[sapply(all_dirs, function(dir) {
  file.exists(file.path(dir, "calibration_drugs.csv"))
})]

# Initialize list to collect data
all_calibration_data <- list()

# Read each calibration_drugs.csv file
for (folder in folders_with_calibration) {
  file_path <- file.path(folder, "calibration_drugs.csv")
  
  calibration_data <- tryCatch({
    read.csv(file_path, stringsAsFactors = FALSE)
  }, error = function(e) {
    message("Error reading file: ", file_path)
    return(NULL)
  })
  
  # Skip if NULL or has no rows
  if (!is.null(calibration_data) && nrow(calibration_data) > 0) {
    calibration_data$source_folder <- folder
    calibration_data$calibration_file <- basename(file_path)
    
    all_calibration_data[[length(all_calibration_data) + 1]] <- calibration_data
  } else {
    message("Skipping empty file: ", file_path)
  }
}

# Combine all data frames
if (length(all_calibration_data) > 0) {
  final_calibration_df <- do.call(rbind, all_calibration_data)
  
  # Save to CSV
  write.csv(final_calibration_df, file = "combined_calibration_drugs_mdl.csv", row.names = FALSE)
  cat("Combined calibration data saved to 'combined_calibration_drugs_mdl.csv'.\n")
} else {
  cat("No valid calibration_drugs.csv files with data rows found.\n")
}
