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


cross_validation_ts <- function(tsdata, sem, fit_function, loo_function, chunk_size = 4, seed = 42, output_df_name = "loodf_cv", results_df_name = "results_semAB_Pred") {
  
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
  chunk_starts <- seq(min_year, max_year - 4, by = chunk_size)  # Skip last 4 years for chunking
  
  # Initialize dataframe to store final test-set predictions
  final_loodf_cv <- data.frame()
  
  # Loop over each chunk (1983 to 2018 or as per your data)
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
  
  # For the last 4 years (2019 to 2022), modify the approach: Predict for each year by excluding that year alone
  last_years <- 2019:2022
  for (year in last_years) {
    # Exclude the current year and create the training set
    train_data_last4 <- tsdata_df %>% filter(Time != year)
    
    # Convert to time series
    train_data_last4_ts <- ts(train_data_last4 %>% select(-Time), start = start(tsdata), frequency = frequency(tsdata))
    
    # Fit the model
    fit_model_last4 <- fit_function(sem = sem,
                                    tsdata = train_data_last4_ts,
                                    control = dsem_control(quiet = TRUE))
    
    # Get model predictions for that year (exclude the year we are predicting)
    loodf_last4 <- loo_function(fit_model_last4, what = "loo", track_progress = FALSE)
    
    # Filter predictions for the current year
    test_set_last4 <- loodf_last4 %>%
      filter(Var1 == year)
    
    # Add to the final predictions
    final_loodf_cv <- bind_rows(final_loodf_cv, test_set_last4)
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
  
  # Assign the predictions and metrics to the global environment with the user-defined names
  assign(output_df_name, final_loodf_cv, envir = .GlobalEnv)  # Assign predictions df
  assign(results_df_name, final_metrics, envir = .GlobalEnv)  # Assign metrics df
  
  # Return both the predictions and the metrics as a list
  return(list(predictions = final_loodf_cv, metrics = final_metrics))
}


# Get predictions for each trophic system and major bay
results <- cross_validation_ts(
  tsdata = AB_Pred_TS,
  sem = semAB_Pred_fulltopdown,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Pred",  
  results_df_name = "results_semAB_Pred"  
)
print(results_semAB_Pred) 

results <- cross_validation_ts(
  tsdata = GB_Pred_TS,
  sem = semGB_Pred_fullbottomup,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Pred",  
  results_df_name = "results_semGB_Pred"  
)
print(results_semGB_Pred)

results <- cross_validation_ts(
  tsdata = GB_Sciaenid_TS,
  sem = semGB_Sciaenid_fullbottomup,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_GB_Sciaenid",
  results_df_name = "results_semGB_Sciaenid"  
)
print(results_semGB_Sciaenid)

results <- cross_validation_ts(
  tsdata = AB_Sciaenid_TS,
  sem = semAB_Sciaenid_notrophics,
  fit_function = dsem,
  loo_function = loo_residuals,
  chunk_size = 4, 
  seed = 42,
  output_df_name = "loodf_AB_Sciaenid",
  results_df_name = "results_semAB_Sciaenid"  
)
print(results_semAB_Sciaenid)

# inspect residuals
samples <- loo_residuals(fit_semAB_Sciaenid_notrophics,  what="samples", track_progress=FALSE)
which_use = which(!is.na(AB_Sciaenid_TS))
fitResp = loo_residuals(fit_semAB_Sciaenid_notrophics, what="loo", track_progress=FALSE)[,'est']
simResp = apply(samples, MARGIN=3, FUN=as.vector)[which_use,]
res = DHARMa::createDHARMa(
  simulatedResponse = simResp,
  observedResponse = unlist(AB_Sciaenid_TS)[which_use],
  fittedPredictedResponse = fitResp )
plot(res)


# calculcate adjusted R2 values
results_semAB_Pred <- as.data.frame(results_semAB_Pred)
results_semGB_Pred <- as.data.frame(results_semGB_Pred)
results_semGB_Sciaenid <- as.data.frame(results_semGB_Sciaenid)
results_semAB_Sciaenid <- as.data.frame(results_semAB_Sciaenid)

results_semAB_Pred <- results_semAB_Pred %>% select(Var2, correlation, RMSE, R_squared)
results_semGB_Pred <- results_semGB_Pred %>% select(Var2, correlation, RMSE, R_squared)
results_semGB_Sciaenid <- results_semGB_Sciaenid %>% select(Var2, correlation, RMSE, R_squared)
results_semAB_Sciaenid <- results_semAB_Sciaenid %>% select(Var2, correlation, RMSE, R_squared)

results_semAB_Pred$adjusted_R2 <- NA
results_semGB_Pred$adjusted_R2 <- NA
results_semGB_Sciaenid$adjusted_R2 <- NA
results_semAB_Sciaenid$adjusted_R2 <- NA


calculate_adjusted_r2 <- function(df, predictors_count, observations_count, specific_var2) {
  df_filtered <- df %>%
    filter(Var2 == specific_var2) %>%
    mutate(
      adjusted_R2 = 1 - ((1 - R_squared) * (observations_count - 1)) / (observations_count - predictors_count - 1)
    )
  df$adjusted_R2[df$Var2 == specific_var2] <- df_filtered$adjusted_R2
  return(df)
}

results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMenhaden_MesquiteBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_AransasBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_CopanoBay")
results_semAB_Pred <- calculate_adjusted_r2(results_semAB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "AllMullet_MesquiteBay")

results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "BullShark_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 5, observations_count = 39, specific_var2 = "AlligatorGar_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Menhaden_EastBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_GalvestonBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_TrinityBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_WestBay")
results_semGB_Pred <- calculate_adjusted_r2(results_semGB_Pred, predictors_count = 7, observations_count = 39, specific_var2 = "Mullet_EastBay")

results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "SpottedSeatrout_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 7, observations_count = 39, specific_var2 = "RedDrum_EastBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_GalvestonBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_TrinityBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_WestBay")
results_semGB_Sciaenid <- calculate_adjusted_r2(results_semGB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_EastBay")


results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "SpottedSeatrout_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "Atlanticcroaker_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "RedDrum_MesquiteBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_AransasBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_CopanoBay")
results_semAB_Sciaenid <- calculate_adjusted_r2(results_semAB_Sciaenid, predictors_count = 4, observations_count = 39, specific_var2 = "BlueCrabSmall_MesquiteBay")

# exporting all 8 of the prediction DFs for easy future access because the prediction models take a while to run.
library(openxlsx)
desktop_path <- "~/Desktop"
write.xlsx(loodf_GB_Pred, file.path(desktop_path, "loodf_GB_Pred.xlsx"))
write.xlsx(results_semGB_Pred, file.path(desktop_path, "results_semGB_Pred.xlsx"))
write.xlsx(loodf_GB_Sciaenid, file.path(desktop_path, "loodf_GB_Sciaenid.xlsx"))
write.xlsx(results_semGB_Sciaenid, file.path(desktop_path, "results_semGB_Sciaenid.xlsx"))

write.xlsx(loodf_AB_Pred, file.path(desktop_path, "loodf_AB_Pred.xlsx"))
write.xlsx(results_semAB_Pred, file.path(desktop_path, "results_semAB_Pred.xlsx"))
write.xlsx(loodf_AB_Sciaenid, file.path(desktop_path, "loodf_AB_Sciaenid.xlsx"))
write.xlsx(results_semAB_Sciaenid, file.path(desktop_path, "results_semAB_Sciaenid.xlsx"))


########## Time Series Plots


plot_observed_vs_predicted <- function(df, response_var, yaxis_title, plot_title, predictions_color, results_df) {
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
  # Extract R² value using base R indexing (avoiding `pull()`)
  r2_value <- results_df$adjusted_R2[results_df$Var2 == response_var]
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
    ) +
    # Add R² annotation in the top right corner
    annotate("text", 
             x = 2001,  # Adjust for positioning near the right edge
             y = max(plot_data$obs, na.rm = TRUE) * 1.1,  # Position inside plot near the top
             label = paste0("R² = ", round(r2_value, 2)), 
             size = 6, 
             fontface = "bold", 
             color = "black")
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
plot_observed_vs_predicted_no_xtext <- function(df, response_var, yaxis_title, plot_title, predictions_color, results_df) {
  p <- plot_observed_vs_predicted(df, response_var, yaxis_title, plot_title, predictions_color, results_df) +
    theme(
      axis.text.x = element_blank(),  
    )
  return(p)
}


# Make plots


#AB PREDATORS
AlligatorGar_CopanoBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "AlligatorGar_CopanoBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#f28482", results_df = results_semAB_Pred)

AlligatorGar_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AlligatorGar_AransasBay",
  predictions_color = "#f28482",yaxis_title = "Alligator Gar CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

AlligatorGar_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AlligatorGar_MesquiteBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

BullShark_CopanoBay <-plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "BullShark_CopanoBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

BullShark_AransasBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "BullShark_AransasBay",
  predictions_color = "#f28482", yaxis_title = "Bull Shark CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

BullShark_MesquiteBay <- plot_observed_vs_predicted(
  df = loodf_AB_Pred, response_var = "BullShark_MesquiteBay",
  predictions_color = "#f28482",yaxis_title = NA,
  plot_title = NA,results_df = results_semAB_Pred)

Mullet_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Copano Bay",results_df = results_semAB_Pred)

Mullet_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_AransasBay",
  predictions_color = "#2a9d8f",yaxis_title = "Mullet CPUE",
  plot_title = "Aransas Bay",results_df = results_semAB_Pred)

Mullet_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMullet_MesquiteBay",
  predictions_color = "#2a9d8f",yaxis_title = NA,
  plot_title = "Mesquite Bay",results_df = results_semAB_Pred)

Menhaden_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMenhaden_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

Menhaden_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMenhaden_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Menhaden CPUE",
  plot_title = NA, results_df = results_semAB_Pred)

Menhaden_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Pred, response_var = "AllMenhaden_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Pred)

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

#GB PREDATORS
AlligatorGar_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "AlligatorGar_GalvestonBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#f28482", results_df = results_semGB_Pred)

AlligatorGar_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_TrinityBay",
  predictions_color = "#f28482", yaxis_title = "Alligator Gar CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

AlligatorGar_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_WestBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

AlligatorGar_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "AlligatorGar_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

BullShark_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "BullShark_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

BullShark_TrinityBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "BullShark_TrinityBay",
  predictions_color = "#f28482", yaxis_title = "Bull Shark CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

BullShark_WestBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "BullShark_WestBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

BullShark_EastBay <- plot_observed_vs_predicted(
  df = loodf_GB_Pred, response_var = "BullShark_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Mullet_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_GalvestonBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Galveston Bay", results_df = results_semGB_Pred)

Mullet_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_TrinityBay",
  predictions_color = "#2a9d8f", yaxis_title = "Mullet CPUE",
  plot_title = "Trinity Bay", results_df = results_semGB_Pred)

Mullet_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_WestBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "West Bay", results_df = results_semGB_Pred)

Mullet_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Mullet_EastBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "East Bay", results_df = results_semGB_Pred)

Menhaden_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Menhaden_GalvestonBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Menhaden_TrinityBay",
  predictions_color = "#2a9d8f", yaxis_title = "Menhaden CPUE",
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Menhaden_WestBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

Menhaden_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Pred, response_var = "Menhaden_EastBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Pred)

GBpredictions_predators <- grid.arrange(Mullet_TrinityBay,
                                        Mullet_GalvestonBay, Mullet_WestBay, Mullet_EastBay ,Menhaden_TrinityBay,
                                        Menhaden_GalvestonBay, Menhaden_WestBay, Menhaden_EastBay, AlligatorGar_TrinityBay,
                                        AlligatorGar_GalvestonBay, AlligatorGar_WestBay, AlligatorGar_EastBay,BullShark_TrinityBay,
                                        BullShark_GalvestonBay, BullShark_WestBay, BullShark_EastBay,ncol = 4, nrow = 4)

ggsave("GBpredictions_predators.png", GBpredictions_predators, dpi = 150, bg = "white",
       width = 3200,
       height = 2000,
       units = "px") 

##GB SCIAENID
RedDrum_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_GalvestonBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#f28482", results_df = results_semGB_Sciaenid)

RedDrum_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_TrinityBay",
  predictions_color = "#f28482", yaxis_title = "Red Drum CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

RedDrum_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_WestBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

RedDrum_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "RedDrum_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_GalvestonBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_GalvestonBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_TrinityBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_TrinityBay",
  predictions_color = "#f28482", yaxis_title = "Spotted Seatrout CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_WestBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_WestBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

SpottedSeatrout_EastBay <- plot_observed_vs_predicted(
  df = loodf_GB_Sciaenid, response_var = "SpottedSeatrout_EastBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

BlueCrabSmall_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_GalvestonBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Galveston Bay", results_df = results_semGB_Sciaenid)

BlueCrabSmall_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_TrinityBay",
  predictions_color = "#2a9d8f", yaxis_title = "Blue Crab CPUE",
  plot_title = "Trinity Bay", results_df = results_semGB_Sciaenid)

BlueCrabSmall_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_WestBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "West Bay", results_df = results_semGB_Sciaenid)

BlueCrabSmall_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "BlueCrabSmall_EastBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "East Bay", results_df = results_semGB_Sciaenid)

Atlanticcroaker_GalvestonBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_GalvestonBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_TrinityBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_TrinityBay",
  predictions_color = "#2a9d8f", yaxis_title = "Atlantic Croaker CPUE",
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_WestBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_WestBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

Atlanticcroaker_EastBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_GB_Sciaenid, response_var = "Atlanticcroaker_EastBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semGB_Sciaenid)

GBpredictions_Sciaenid <- grid.arrange(BlueCrabSmall_TrinityBay,
                                       BlueCrabSmall_GalvestonBay, BlueCrabSmall_WestBay, BlueCrabSmall_EastBay,Atlanticcroaker_TrinityBay,
                                       Atlanticcroaker_GalvestonBay, Atlanticcroaker_WestBay, Atlanticcroaker_EastBay,RedDrum_TrinityBay, 
                                       RedDrum_GalvestonBay,RedDrum_WestBay, RedDrum_EastBay,SpottedSeatrout_TrinityBay,
                                       SpottedSeatrout_GalvestonBay, SpottedSeatrout_WestBay, SpottedSeatrout_EastBay,
                                       ncol = 4, nrow = 4)

ggsave("GBpredictions_Sciaenid.png", GBpredictions_Sciaenid, dpi = 150, bg = "white",
       width = 3200,
       height = 2000,
       units = "px")

##AB SCIAENID
RedDrum_CopanoBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_CopanoBay",
  yaxis_title = NA, plot_title = NA,
  predictions_color = "#f28482", results_df = results_semAB_Sciaenid)

RedDrum_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_AransasBay",
  predictions_color = "#f28482", yaxis_title = "Red Drum CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

RedDrum_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "RedDrum_MesquiteBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

SpottedSeatrout_CopanoBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_CopanoBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

SpottedSeatrout_AransasBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_AransasBay",
  predictions_color = "#f28482", yaxis_title = "Spotted Seatrout CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

SpottedSeatrout_MesquiteBay <- plot_observed_vs_predicted(
  df = loodf_AB_Sciaenid, response_var = "SpottedSeatrout_MesquiteBay",
  predictions_color = "#f28482", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

BlueCrabSmall_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Copano Bay", results_df = results_semAB_Sciaenid)

BlueCrabSmall_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Blue Crab CPUE",
  plot_title = "Aransas Bay", results_df = results_semAB_Sciaenid)

BlueCrabSmall_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "BlueCrabSmall_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = "Mesquite Bay", results_df = results_semAB_Sciaenid)

Atlanticcroaker_CopanoBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_CopanoBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

Atlanticcroaker_AransasBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_AransasBay",
  predictions_color = "#2a9d8f", yaxis_title = "Atlantic Croaker CPUE",
  plot_title = NA, results_df = results_semAB_Sciaenid)

Atlanticcroaker_MesquiteBay <- plot_observed_vs_predicted_no_xtext(
  df = loodf_AB_Sciaenid, response_var = "Atlanticcroaker_MesquiteBay",
  predictions_color = "#2a9d8f", yaxis_title = NA,
  plot_title = NA, results_df = results_semAB_Sciaenid)

ABpredictions_Sciaenid <- grid.arrange(
  BlueCrabSmall_AransasBay, BlueCrabSmall_CopanoBay, BlueCrabSmall_MesquiteBay,
  Atlanticcroaker_AransasBay, Atlanticcroaker_CopanoBay, Atlanticcroaker_MesquiteBay,
  RedDrum_AransasBay, RedDrum_CopanoBay, RedDrum_MesquiteBay,
  SpottedSeatrout_AransasBay, SpottedSeatrout_CopanoBay, SpottedSeatrout_MesquiteBay,
  ncol = 3, nrow = 4
)

ggsave("ABpredictions_Sciaenid.png", ABpredictions_Sciaenid, dpi = 150, bg = "white",
       width = 3200,
       height = 2000,
       units = "px")


#summary plot for adjusted R2 values 

#merge the data frames into one
squared_df <- bind_rows(results_semAB_Pred, results_semGB_Pred, results_semAB_Sciaenid, results_semGB_Sciaenid)
squared_df <- squared_df %>%
  separate(Var2, into = c("species", "bay"), sep = "_", remove = FALSE) %>%
  mutate(bay = gsub("([a-z])([A-Z])", "\\1 \\2", bay)) %>%  
  mutate(species = case_when(
    species == "AllMullet" ~ "Mullet",
    species == "AllMenhaden" ~ "Menhaden",
    species == "BlueCrabSmall" ~ "BlueCrab",
    TRUE ~ species  
  )) %>%
  filter(!species %in% c("PDSI", "Salinity")) %>%  
  mutate(major_bay = case_when(
    bay %in% c("Aransas Bay", "Copano Bay", "Mesquite Bay") ~ "Aransas Bay",
    bay %in% c("Galveston Bay", "West Bay", "East Bay", "Trinity Bay") ~ "Galveston Bay",
    TRUE ~ NA_character_  
  )) %>%
  mutate(trophic_system = case_when(
    species %in% c("BullShark", "AlligatorGar", "Mullet", "Menhaden") ~ "Keystone Predator",
    species %in% c("RedDrum", "SpottedSeatrout", "BlueCrab", "Atlanticcroaker") ~ "Sciaenid",
    TRUE ~ NA_character_  
  ))

bay_order <- c("Aransas Bay", "Copano Bay", "Mesquite Bay", "Galveston Bay", "Trinity Bay", "East Bay", "West Bay")
squared_df$bay <- factor(squared_df$bay, levels = bay_order)
squared_df$major_bay <- factor(squared_df$major_bay, levels = c("Aransas Bay", "Galveston Bay"))

species_colors <- c(
  "BullShark" = "#f28482",
  "AlligatorGar" = "indianred3",
  "RedDrum" = "#f28482",
  "SpottedSeatrout" = "indianred3",
  "BlueCrab" = "#2a9d8f",
  "Atlanticcroaker" = "cadetblue3",
  "Mullet" = "#2a9d8f",
  "Menhaden" = "cadetblue3"
)

plot_keystone <- ggplot(squared_df %>% filter(trophic_system == "Keystone Predator"), 
                        aes(x = bay, y = adjusted_R2, fill = species)) +
  geom_point(size = 9, alpha = 0.7, shape = 21, color = "black", stroke = 1.2) +
  scale_fill_manual(values = species_colors) +  
  scale_y_continuous(limits = c(-0.35, 0.75), breaks = seq(-0.35, 0.75, by = 0.1)) +  
  theme_bw() +      theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                          axis.title.x = element_text(size = 12),  
                          axis.title.y = element_text(size = 12),  
                          legend.position = "top", 
                          legend.direction = "horizontal",  
                          legend.title = element_blank(),  # Removes the legend title
                          legend.text = element_text(size = 12),  
                          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Centers title above legend
                          axis.text.y = element_text(size = 12)) + 
  labs(title = "Keystone Predator Trophic System",  # Title above the legend
       x = "Minor Bay", 
       y = "Adjusted R²")

plot_sciaenid <- ggplot(squared_df %>% filter(trophic_system == "Sciaenid"), 
                        aes(x = bay, y = adjusted_R2, fill = species)) +
  geom_point(size = 9, alpha = 0.7, shape = 21, color = "black", stroke = 1.2) +
  scale_fill_manual(values = species_colors) +  
  scale_y_continuous(limits = c(-0.35, 0.75), breaks = seq(-0.35, 0.75, by = 0.1)) +  
  theme_bw() +    theme(axis.text.x = element_text(angle = 45, hjust = 1),  
                        axis.title.x = element_text(size = 12),  
                        axis.title.y = element_text(size = 12),  
                        legend.position = "top", 
                        legend.direction = "horizontal",  
                        legend.title = element_blank(),  # Removes the legend title
                        legend.text = element_text(size = 12),  
                        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Centers title above legend
                        axis.text.y = element_text(size = 12)) + 
  labs(title = "Sciaenid Trophic System",  # Title above the legend
       x = "Minor Bay", 
       y = "Adjusted R²")

r2plot<-grid.arrange(plot_keystone, plot_sciaenid, ncol = 2)

ggsave("r2plot.png", r2plot, dpi = 150, bg = "white",
       width = 2000,
       height = 1000,
       units = "px")
