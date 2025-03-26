library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(gridExtra)
library(grid)


# Rolling-Block 4-Year Cross-Validation Function
cross_validation_ts <- function(tsdata, sem, fit_function, loo_function, chunk_size = 4, seed = 42, output_df_name = "loodf_cv") {
  
  # Set the seed for reproducibility
  set.seed(seed)
  
  # Extract time index from the ts object
  time_index <- as.numeric(time(tsdata))  # Converts ts time attribute into a numeric vector
  
  # Convert tsdata into a dataframe (keeping time index as first column)
  tsdata_df <- as.data.frame(tsdata) %>%
    mutate(Time = time_index)  # Add extracted time column
  
  # Identify the start and end year
  min_year <- min(tsdata_df$Time)
  max_year <- max(tsdata_df$Time)
  
  # Define chunk start years (we move in 4-year blocks)
  chunk_starts <- seq(min_year, max_year, by = chunk_size)
  
  # Initialize dataframe to store final test-set predictions
  final_loodf_cv <- data.frame()
  
  # Loop over each 4-year block
  for (start_year in chunk_starts) {
    
    # Define the 4-year chunk to remove
    test_years <- seq(start_year, start_year + chunk_size - 1)
    
    # Ensure we don’t go beyond the dataset range
    test_data <- tsdata_df %>% filter(Time %in% test_years)
    train_data <- tsdata_df %>% filter(!Time %in% test_years)
    
    # Convert train data back to time series (ts object)
    train_data_ts <- ts(train_data %>% select(-Time), start = start(tsdata), frequency = frequency(tsdata))
    
    # Fit the model using the training data (now as ts object)
    fit_model <- fit_function(sem = sem,
                              tsdata = train_data_ts,
                              control = dsem_control(quiet = TRUE))
    
    # Get model predictions using loo_function (e.g., loo_residuals)
    loodf_cv <- loo_function(fit_model, what = "loo", track_progress = FALSE)
    
    # Keep only test set predictions (where those years were missing in training)
    test_set_predictions <- loodf_cv %>%
      filter(Var1 %in% test_years)  # Only keep years that were omitted
    
    # Store these test set predictions
    final_loodf_cv <- bind_rows(final_loodf_cv, test_set_predictions)
  }
  
  # 🛠️ Fix: Remove NA values before calculating RMSE and correlation
  final_loodf_cv <- final_loodf_cv %>%
    filter(!is.na(obs) & !is.na(est))  # Ensure no NA values
  
  # 🛠️ Fix: Calculate RMSE and correlation **per variable**
  final_metrics <- final_loodf_cv %>%
    group_by(Var2) %>%
    summarise(
      correlation = cor(obs, est, use = "complete.obs"),
      RMSE = sqrt(mean((obs - est)^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # ✅ Assign final_loodf_cv (only test set predictions) to the global environment with the specified name
  assign(output_df_name, final_loodf_cv, envir = .GlobalEnv)
  
  # Return the final RMSE and correlation per variable
  return(final_metrics)
}


# Function to plot observed vs predicted values for a specified variable
plot_observed_vs_predicted <- function(df, response_var, yaxis_title, plot_title, predictions_color) {
  
  # Filter the data for the specified response variable
  plot_data <- df %>% 
    filter(Var2 == response_var) %>%
    select(Var1, obs, est, se)
  
  # Ensure Var1 (Year) is numeric
  plot_data <- plot_data %>% mutate(Var1 = as.numeric(Var1))
  
  # Calculate confidence intervals (95% CI) for the predictions
  plot_data <- plot_data %>%
    mutate(
      lower_ci = est - 1.96 * se,
      upper_ci = est + 1.96 * se
    )
  
  # Create the base plot
  p <- ggplot(plot_data, aes(x = Var1)) +
    # Observed values (dots + line)
    geom_point(aes(y = obs), color = "black") +
    geom_line(aes(y = obs), color = "black", linewidth = 1) +
    # Predicted values (dots + line)
    geom_point(aes(y = est), color = predictions_color, size=4) +
    geom_line(aes(y = est), color = predictions_color, linewidth = 2) +
    # Confidence interval shading
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = predictions_color, alpha = 0.25) +
    # Formatting
    scale_x_continuous(breaks = seq(min(plot_data$Var1), max(plot_data$Var1), by = 2)) +  # Show every other year
    labs(
      title = plot_title,
      y = ifelse(is.na(yaxis_title), "", yaxis_title),  # Set y-axis title to "" if NA
      x = "Year"
    ) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),  # Rotate x-axis labels
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 16, hjust = 0.5)
    )
  
  # If yaxis_title is NA, remove the y-axis title completely
  if (is.na(yaxis_title)) {
    p <- p + theme(axis.title.y = element_blank())
  }
  
  # If plot title is NA, remove it completely
  if (is.na(plot_title)) {
    p <- p + theme(plot.title = element_blank())
  }
  
  # Return the plot
  return(p)
}


# Function to remove x-axis text while keeping the optional y-axis title removal
plot_observed_vs_predicted_no_xtext <- function(df, response_var, yaxis_title, plot_title, predictions_color) {
  p <- plot_observed_vs_predicted(df, response_var, yaxis_title, plot_title, predictions_color) +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis text
      axis.title.x = element_blank()  # Remove x-axis title
    )
  return(p)
}



###################


#### Testing the prediction skill function on the full bottom-up model for Aransas Bay for the Keystone Predator system. 

#AB-keystone-full bottom up
results_semAB_Pred <- cross_validation_ts(
  tsdata = AB_Pred_TS,
  sem = semAB_Pred_fullbottomup,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Pred" 
)
print(results_semAB_Pred)
desktop_path <- file.path(Sys.getenv("HOME"), "Desktop", "results_semAB_Pred.xlsx")
write_xlsx(results_semAB_Pred__fullbottomup, path = desktop_path)

