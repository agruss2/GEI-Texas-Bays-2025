library(dplyr)
library(tidyr)
library(ggplot2)
library(writexl)
library(gridExtra)
library(grid)
library(dsem)
library(ggpubr)
library(ggraph)
library(phylopath)
library(ggdag)
library(readxl)
library(qgraph)
library(DHARMa)



########## Get predictions for each year, variable, minor bay, trophic system and major bay


# Rolling-Block 4-Year Cross-Validation Function with 3-Year Last Chunk
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
  total_obs <- nrow(tsdata_df)
  
  # Define chunk start years (adjusting the last chunk if needed)
  chunk_starts <- seq(min_year, max_year, by = chunk_size)
  
  # Adjust for last chunk being 3 years if total observations is not divisible by 4
  if ((total_obs %% chunk_size) != 0) {
    chunk_starts[length(chunk_starts)] <- max_year - 2  # Last chunk starts 3 years before the last year
  }
  
  # Initialize dataframe to store final test-set predictions
  final_loodf_cv <- data.frame()
  
  # Loop over each chunk
  for (start_year in chunk_starts) {
    
    # Determine the chunk size dynamically
    current_chunk_size <- ifelse(start_year == max_year - 2, 3, chunk_size)
    
    # Define the chunk to remove
    test_years <- seq(start_year, start_year + current_chunk_size - 1)
    
    # Ensure we don’t go beyond the dataset range
    test_data <- tsdata_df %>% filter(Time %in% test_years)
    train_data <- tsdata_df %>% filter(!Time %in% test_years)
    
    # Convert train data back to time series (ts object)
    train_data_ts <- ts(train_data %>% select(-Time), start = start(tsdata), frequency = frequency(tsdata))
    
    # Fit the model using the training data (now as ts object)
    fit_model <- fit_function(sem = sem,
                              tsdata = train_data_ts,
                              control = dsem_control(quiet = TRUE))
    
    # Get model predictions using loo_function
    loodf_cv <- loo_function(fit_model, what = "loo", track_progress = FALSE)
    
    # Keep only test set predictions (where those years were missing in training)
    test_set_predictions <- loodf_cv %>%
      filter(Var1 %in% test_years)  # Only keep years that were omitted
    
    # Store these test set predictions
    final_loodf_cv <- bind_rows(final_loodf_cv, test_set_predictions)
  }
  
  # Remove NA values before calculating metrics
  final_loodf_cv <- final_loodf_cv %>%
    filter(!is.na(obs) & !is.na(est))
  
  # Calculate RMSE, correlation, and R-squared per variable
  final_metrics <- final_loodf_cv %>%
    group_by(Var2) %>%
    summarise(
      correlation = cor(obs, est, use = "complete.obs"),
      RMSE = sqrt(mean((obs - est)^2, na.rm = TRUE)),
      R_squared = 1 - (sum((obs - est)^2, na.rm = TRUE) / sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Assign final predictions dataframe to global environment
  assign(output_df_name, final_loodf_cv, envir = .GlobalEnv)
  
  # Return final metrics
  return(final_metrics)
}

results_semGB_Pred <- cross_validation_ts(
  tsdata = GB_Pred_TS,
  sem = semGB_Pred_fulltopdown,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Pred" 
)
print(results_semGB_Pred)

results_semGB_Sciaenid <- cross_validation_ts(
  tsdata = GB_Sciaenid_TS,
  sem = semGB_Sciaenid_fullbottomup,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Sciaenid" 
)
print(results_semGB_Sciaenid)

results_semAB_Sciaenid <- cross_validation_ts(
  tsdata = AB_Sciaenid_TS,
  sem = semAB_Sciaenid_notrophics,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Sciaenid" 
)
print(results_semAB_Sciaenid)

results_semAB_Pred <- cross_validation_ts(
  tsdata = AB_Pred_TS,
  sem = semAB_Pred_fulltopdown,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Pred" 
)
print(results_semAB_Pred)

########## Time Series Plots



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
      axis.title.y = element_text(size = 20),
      axis.title.x = element_text(size = 14),
      plot.title = element_text(size = 20, hjust = 0.5)
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
      axis.text.x = element_blank(),  
    )
  return(p)
}


# Make time series plots


# Function to extract R² value and add to plots
add_r_squared <- function(plot, species_bay) {
  r2_value <- results_semAB_Pred %>% filter(Var2 == species_bay) %>% pull(R_squared)
  plot + annotate("text", x = 2007, y = Inf, label = paste0("R² = ", round(r2_value, 2)),
                  size = 8, hjust = 1.2, vjust = 1.2, fontface = "bold")
}

# AB key stone predators first (I have yet to do the 3 other trophic models/major bays)

AlligatorGar_CopanoBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AlligatorGar_CopanoBay",
  predictions_color = "#f28482",
  yaxis_title = NA,
  plot_title = NA
), "AlligatorGar_CopanoBay")

AlligatorGar_AransasBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AlligatorGar_AransasBay",
  predictions_color = "#f28482",
  yaxis_title = "Alligator Gar CPUE",
  plot_title = NA
), "AlligatorGar_AransasBay")

AlligatorGar_MesquiteBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AlligatorGar_MesquiteBay",
  predictions_color = "#f28482",
  yaxis_title = NA,
  plot_title = NA
), "AlligatorGar_MesquiteBay")

BullShark_CopanoBay <- add_r_squared(plot_observed_vs_predicted(
  df = loodf_AB_Pred,  
  response_var = "BullShark_CopanoBay",
  predictions_color = "#f28482",
  yaxis_title = NA,
  plot_title = NA
), "BullShark_CopanoBay")

BullShark_AransasBay <- add_r_squared(plot_observed_vs_predicted(
  df = loodf_AB_Pred,  
  response_var = "BullShark_AransasBay",
  predictions_color = "#f28482",
  yaxis_title = "Bull Shark CPUE",
  plot_title = NA
), "BullShark_AransasBay")

BullShark_MesquiteBay <- add_r_squared(plot_observed_vs_predicted(
  df = loodf_AB_Pred,  
  response_var = "BullShark_MesquiteBay",
  predictions_color = "#f28482",
  yaxis_title = NA,
  plot_title = NA
), "BullShark_MesquiteBay")

Mullet_CopanoBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMullet_CopanoBay",
  predictions_color = "#2a9d8f",
  yaxis_title = NA,
  plot_title = "Copano Bay"
), "AllMullet_CopanoBay")

Mullet_AransasBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMullet_AransasBay",
  predictions_color = "#2a9d8f",
  yaxis_title = "Mullet CPUE",
  plot_title = "Aransas Bay"
), "AllMullet_AransasBay")

Mullet_MesquiteBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMullet_MesquiteBay",
  predictions_color = "#2a9d8f",
  yaxis_title = NA,
  plot_title = "Mesquite Bay"
), "AllMullet_MesquiteBay")

Menhaden_CopanoBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMenhaden_CopanoBay",
  predictions_color = "#2a9d8f",
  yaxis_title = NA,
  plot_title = NA
), "AllMenhaden_CopanoBay")

Menhaden_AransasBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMenhaden_AransasBay",
  predictions_color = "#2a9d8f",
  yaxis_title = "Menhaden CPUE",
  plot_title = NA
), "AllMenhaden_AransasBay")

Menhaden_MesquiteBay <- add_r_squared(plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred,  
  response_var = "AllMenhaden_MesquiteBay",
  predictions_color = "#2a9d8f",
  yaxis_title = NA,
  plot_title = NA
), "AllMenhaden_MesquiteBay")

ABpredictions_predators <- grid.arrange(
  Mullet_AransasBay, Mullet_CopanoBay, Mullet_MesquiteBay,
  Menhaden_AransasBay, Menhaden_CopanoBay, Menhaden_MesquiteBay,
  AlligatorGar_AransasBay, AlligatorGar_CopanoBay, AlligatorGar_MesquiteBay,
  BullShark_AransasBay, BullShark_CopanoBay, BullShark_MesquiteBay,
  ncol = 3, nrow = 4
)

ggsave("ABpredictions_predators.png", ABpredictions_predators, dpi = 150, bg = "white",
       width = 3200,
       height = 2000,
       units = "px") 



