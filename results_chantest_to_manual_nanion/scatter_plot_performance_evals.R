library(ggplot2)
library(readr)
results <- data.frame()
for (i in 6:12) {
  filename <- paste0("results_performance_evals_check_",i,"/summary_tms_",i,"_metrics.csv")
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

# Scatter plot
ggplot(results, aes(x = Th1_changes, y = Th2_changes, color = factor(number_of_drugs))) + # Use `factor()` if the column is categorical
  geom_point(size = 3) +                               # Set point size
  labs(title = "Scatter Plot", 
       x = "Th1_changes", 
       y = "Th2_changes", 
       color = "Number of drugs") +
  theme_minimal()  
ggplot(results, aes(x = Th1_changes, y = Th2_changes, color = Normalized_logLik)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") + # Color gradient
  labs(title = "Scatter Plot with Color Bar",
       x = "Th1_changes",
       y = "Th2_changes",
       color = "Normalized_logLik") + # Label for the color bar
  theme_minimal()


# Revised scatter plot
ggplot(results, aes(x = Th1_changes_new, y = Th2_changes_new, color = factor(number_of_drugs))) + # Use `factor()` if the column is categorical
  geom_point(size = 3) +                               # Set point size
  labs(title = "Scatter Plot", 
       x = "Th1_changes_new", 
       y = "Th2_changes_new", 
       color = "Number of drugs") +
  theme_minimal()                                      # Use a clean theme

# Scatter plot with color bar
ggplot(results, aes(x = Th1_changes_new, y = Th2_changes_new, color = loglike_rank_new)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") + # Color gradient
  labs(title = "Scatter Plot with Color Bar",
       x = "Th1_changes_new",
       y = "Th2_changes_new",
       color = "loglike_rank_new") + # Label for the color bar
  theme_minimal()

# Scatter plot with color bar
ggplot(results, aes(x = Th1_changes_new, y = loglike_rank_new, color = Th2_changes_new)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") + # Color gradient
  labs(title = "Scatter Plot with Color Bar",
       x = "Th1_changes_new",
       y = "loglike_rank_new",
       color = "Th2_changes") + # Label for the color bar
  theme_minimal()

# Scatter plot with color bar
ggplot(results, aes(x = Th2_changes_new, y = loglike_rank_new, color = Th1_changes_new)) +
  geom_point(size = 3) +
  scale_color_gradient(low = "blue", high = "red") + # Color gradient
  labs(title = "Scatter Plot with Color Bar",
       x = "Th2_changes_new",
       y = "loglike_rank_new",
       color = "Th1_changes") + # Label for the color bar
  theme_minimal()