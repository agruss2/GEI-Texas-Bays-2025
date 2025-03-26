process_PDSI_data <- function(df, export_name) {
  library(dplyr)
  
  # Step 1: Filter the dataframe to keep only rows where Var2 is 'PDSI'
  df_PDSI <- df %>%
    filter(Var2 == "PDSI") %>%
    arrange(desc(obs))  # Sort by observed PDSI values in descending order
  
  # Step 2: Identify the 3 years with the highest and lowest observed PDSI values
  highest_years <- head(df_PDSI, 3)  # Top 3 highest PDSI years
  lowest_years <- tail(df_PDSI, 3)   # Bottom 3 lowest PDSI years
  
  # Get the years for the highest and lowest PDSI
  highest_years_list <- unique(highest_years$Var1)
  lowest_years_list <- unique(lowest_years$Var1)
  
  # Combine the highest and lowest years
  years_to_keep <- c(highest_years_list, lowest_years_list)
  
  # Step 3: Filter the original dataframe to keep rows for these years across all variables
  df_PDSI_filtered <- df %>%
    filter(Var1 %in% years_to_keep)
  
  top_bottom_years <- df_PDSI_filtered %>%
    filter(Var2 == "PDSI") %>%
    select(Var1, obs) %>%
    arrange(desc(obs)) %>%
    slice(c(1:3, (n() - 2):n()))  # Get the 3 highest and 3 lowest PDSI years
  
  # Step 4: Create the 'weather' column based on the PDSI values of those years
  df_PDSI_filtered <- df_PDSI_filtered %>%
    mutate(weather = ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs < 0], "Severely Dry", 
                            ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs > 0], "Severely Wet", NA)))
  
  # Step 5: Calculate mean and standard error of 'est' for each species/minor bay and weather condition
  df_summary <- df_PDSI_filtered %>%
    group_by(Var2, weather) %>%
    summarise(
      mean_est = mean(est, na.rm = TRUE),  # Calculate mean of 'est'
      se_est = sd(est, na.rm = TRUE) / sqrt(n())  # Calculate standard error of 'est'
    ) %>%
    ungroup()
  
  # Step 6: Identify all species and minor bays dynamically
  excluded_vars <- c("PDSI")  # Exclude PDSI from species list
  all_species <- unique(df$Var2[!df$Var2 %in% excluded_vars])  # Identify species
  minor_bays <- unique(sub(".*_", "", all_species))  # Extract minor bay names
  
  # Step 7: Get mean PDSI values for Severely Dry and Severely Wet years
  mean_PDSI_dry <- mean(df_PDSI_filtered$obs[df_PDSI_filtered$Var2 == "PDSI" & df_PDSI_filtered$weather == "Severely Dry"], na.rm = TRUE)
  mean_PDSI_wet <- mean(df_PDSI_filtered$obs[df_PDSI_filtered$Var2 == "PDSI" & df_PDSI_filtered$weather == "Severely Wet"], na.rm = TRUE)
  
  # Step 8: Assign the summary tibble to the user-specified name
  assign(export_name, df_summary, envir = .GlobalEnv)
  
  # Return a list of results
  return(list(
    summary = df_summary,
    mean_PDSI_dry = mean_PDSI_dry,
    mean_PDSI_wet = mean_PDSI_wet,
    species = all_species,
    minor_bays = minor_bays
  ))
}

result <- process_PDSI_data(loodf_AB_Pred, "loodf_AB_Pred_PDSI_forplotting")

# Get mean PDSI values for Severely Dry and Severely Wet years
mean_PDSI_dry_AB <- mean(loodf_AB_Pred_PDSI_forplotting$mean_est[loodf_AB_Pred_PDSI_forplotting$Var2 == "PDSI" & loodf_AB_Pred_PDSI_forplotting$weather == "Severely Dry"], na.rm = TRUE)
mean_PDSI_wet_AB <- mean(loodf_AB_Pred_PDSI_forplotting$mean_est[loodf_AB_Pred_PDSI_forplotting$Var2 == "PDSI" & loodf_AB_Pred_PDSI_forplotting$weather == "Severely Wet"], na.rm = TRUE)


# Remove Salinity-related rows
exported_df_filtered  <- loodf_AB_Pred_PDSI_forplotting%>%
  filter(!Var2 %in% c("PDSI", "Salinity_CopanoBay", "Salinity_MesquiteBay", "Salinity_AransasBay")) %>%
  mutate(
    Var2 = recode(Var2,
                  "BullShark_AransasBay" = "Bull Shark (Aransas Bay)",
                  "AllMenhaden_AransasBay" = "Menhaden (Aransas Bay)",
                  "AlligatorGar_AransasBay" = "Alligator Gar (Aransas Bay)",
                  "AllMullet_AransasBay" = "Mullet (Aransas Bay)",
                  "BullShark_CopanoBay" = "Bull Shark (Copano Bay)",
                  "AllMenhaden_CopanoBay" = "Menhaden (Copano Bay)",
                  "AlligatorGar_CopanoBay" = "Alligator Gar (Copano Bay)",
                  "AllMullet_CopanoBay" = "Mullet (Copano Bay)",
                  "BullShark_MesquiteBay" = "Bull Shark (Mesquite Bay)",
                  "AllMenhaden_MesquiteBay" = "Menhaden (Mesquite Bay)",
                  "AlligatorGar_MesquiteBay" = "Alligator Gar (Mesquite Bay)",
                  "AllMullet_MesquiteBay" = "Mullet (Mesquite Bay)"
    )
  )


# Create the bar plot
PDSI_plot <- ggplot(exported_df_filtered, aes(x = Var2, y = mean_est, fill = weather)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +
  
  # Create error bars extending only in the correct direction
  geom_errorbar(aes(
    ymin = ifelse(mean_est > 0, mean_est, mean_est - se_est),  
    ymax = ifelse(mean_est > 0, mean_est + se_est, mean_est)   
  ), width = 0.25, position = position_dodge(0.7)) +  
  
  # Color for weather categories
  scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) +
  
  theme_bw() +  
  
  # Axis and label customization
  labs(
    x = NULL,  
    y = "Mean Standardized CPUE",  
    fill = NULL  
  ) +
  
  # Additional theme adjustments
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    plot.title = element_blank(),  
    legend.position = c(0.85, 0.90),  
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 16),
    legend.background = element_blank(),
    legend.key.height = unit(1, "cm")  
  ) +
  # Annotate the mean PDSI values for both Severely Dry and Severely Wet
  annotate("text", 
           x = 9.6,  # Adjust x position as needed to be right under the legend
           y = max(exported_df_filtered$mean_est) - 0.23,  # Adjust y position to place below the legend
           label = paste("Severely Dry Mean PDSI = ", round(mean_PDSI_dry_AB, 1)),
           size = 4.5, hjust = 0, vjust = 0.5, color = "black") +
  annotate("text", 
           x = 9.6,  # Adjust x position for the second annotation
           y = max(exported_df_filtered$mean_est) - 0.37,  # Place below the first annotation
           label = paste("Severely Wet Mean PDSI = ", round(mean_PDSI_wet_AB, 1)),
           size = 4.5, hjust = 0, vjust = 0.5, color = "black") 

print(PDSI_plot)


