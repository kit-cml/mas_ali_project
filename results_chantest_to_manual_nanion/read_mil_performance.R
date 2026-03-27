# Set the parent directory
parent_dir <- "results_validate_calibration_drugs_chantest_to_nanion"

# Get all subdirectories
all_dirs <- list.dirs(parent_dir, recursive = TRUE, full.names = TRUE)

# Find folders that contain 'calibration_drugs.csv'
folders_with_file <- all_dirs[sapply(all_dirs, function(dir) {
  file.exists(file.path(dir, "calibration_drugs.csv"))
})]

# Initialize an empty list to store data frames
all_metrics_list <- list()

# Loop through each folder
for (folder in folders_with_file) {
  # Find subfolders whose name starts with 'mil_performance_evals'
  eval_dirs <- list.dirs(folder, full.names = TRUE, recursive = FALSE)
  eval_dirs <- eval_dirs[grepl("^.*/mil_performance_evals_[0-9]*$", eval_dirs)]
  
  # If found any matching 'mil_performance_evals' subfolder
  for (eval_dir in eval_dirs) {
    # List all files ending with 'metrics.csv' inside eval_dir
    metrics_files <- list.files(eval_dir, pattern = "metrics\\.csv$", full.names = TRUE)
    
    if (length(metrics_files) > 0) {
      # Loop through each metrics file
      for (metrics_file in metrics_files) {
        # Read the metrics file
        metrics_data <- tryCatch({
          read.csv(metrics_file, stringsAsFactors = FALSE)
        }, error = function(e) {
          message("Error reading: ", metrics_file)
          NULL
        })
        
        # If successfully read, add folder info
        if (!is.null(metrics_data)) {
          metrics_name <- basename(metrics_file)  # <-- Only filename
          metrics_data$tms_name <- basename(folder)    # Folder where calibration_drugs.csv is found
          metrics_data$metrics_name <- metrics_name  # <-- Only the file name
          all_metrics_list[[length(all_metrics_list) + 1]] <- metrics_data
        }
      }
    }
  }
}

# Combine all data frames into one
if (length(all_metrics_list) > 0) {
  final_df <- do.call(rbind, all_metrics_list)
  
  # Save to CSV
  write.csv(final_df, file = "chantest_to_nanion.csv", row.names = FALSE)
  
  cat("Combined metrics file created successfully.\n")
} else {
  cat("No metrics files found.\n")
}
