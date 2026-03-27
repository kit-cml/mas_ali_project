# Set the parent directory where to search
parent_dir <- "results_validate_calibration_drugs_chantest_to_nanion"

# Get all subdirectories
all_dirs <- list.dirs(parent_dir, recursive = TRUE, full.names = TRUE)

# Find folders that contain 'calibration_drugs.csv'
folders_with_file <- all_dirs[sapply(all_dirs, function(dir) {
  file.exists(file.path(dir, "calibration_drugs.csv"))
})]

# Count how many
count <- length(folders_with_file)
cat("Total number of folders containing 'calibration_drugs.csv':", count, "\n")

# Create a dataframe to store the folder paths
result_df <- data.frame(folder_path = folders_with_file, stringsAsFactors = FALSE)

# Save to CSV
write.csv(result_df, file = "chantest_to_nanion.csv", row.names = FALSE)
