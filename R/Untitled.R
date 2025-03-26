create_species_plot <- function(df, species_col, minor_bays_col, species, minor_bays) {
  
  # Step 1: Filter the dataframe to keep only rows where Var2 is 'PDSI'
  df_filtered <- df %>%
    filter(Var2 == "PDSI") %>%
    arrange(desc(obs))
  
  # Step 2: Identify the 3 years with the highest and lowest observed PDSI values
  highest_years <- head(df_filtered, 3) 
  lowest_years <- tail(df_filtered, 3)   
  
  # Get the years for the highest and lowest PDSI
  highest_years_list <- unique(highest_years$Var1)
  lowest_years_list <- unique(lowest_years$Var1)
  
  # Combine the highest and lowest years
  years_to_keep <- c(highest_years_list, lowest_years_list)
  
  # Step 3: Filter the original dataframe to keep rows for these years across all variables (3 highest and 3 lowest PDSI years)
  df_filtered <- df %>%
    filter(Var1 %in% years_to_keep)
  
  top_bottom_years <- df_filtered %>%
    filter(Var2 == "PDSI") %>%
    select(Var1, obs) %>%
    arrange(desc(obs)) %>%
    slice(c(1:3, (n() - 2):n()))  
  
  # Create the 'weather' column based on the PDSI values of those years
  df_filtered <- df_filtered %>%
    mutate(weather = ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs < 0], "Severely Dry", 
                            ifelse(Var1 %in% top_bottom_years$Var1[top_bottom_years$obs > 0], "Severely Wet", NA)))
  
  # Calculate mean and standard error of 'est' for each variable (Var2) and weather condition
  df_mean <- df_filtered %>%
    group_by(Var2, weather) %>%
    summarise(
      mean_est = mean(est, na.rm = TRUE),  
      se_est = sd(est, na.rm = TRUE) / sqrt(n())  
    )
  
  # Modify the variable names directly in the dataframe to include desired formatting
  df_mean_filtered <- df_mean %>%
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
  
  # Filter the dataset to include only the selected species and minor bays
  df_mean_filtered <- df_mean_filtered %>%
    filter(grepl(paste(species, collapse = "|"), Var2)) %>%
    filter(grepl(paste(minor_bays, collapse = "|"), Var2))
  
  # Get mean PDSI values for Severely Dry and Severely Wet years
  mean_PDSI_dry <- mean(df_filtered$obs[df_filtered$Var2 == "PDSI" & df_filtered$weather == "Severely Dry"], na.rm = TRUE)
  mean_PDSI_wet <- mean(df_filtered$obs[df_filtered$Var2 == "PDSI" & df_filtered$weather == "Severely Wet"], na.rm = TRUE)
  
  # Create the bar plot with error bars extending only in the correct direction
  plot <- ggplot(df_mean_filtered, aes(x = Var2, y = mean_est, fill = weather)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7, color = "black", alpha = 0.8) +  # Adjust alpha and black border
    geom_errorbar(aes(
      ymin = ifelse(mean_est > 0, mean_est, mean_est - se_est),  
      ymax = ifelse(mean_est > 0, mean_est + se_est, mean_est)   
    ), width = 0.25, position = position_dodge(0.7)) +
    scale_fill_manual(values = c("Severely Dry" = "#e9c46a", "Severely Wet" = "#264653")) + 
    theme_bw() +
    labs(
      x = NULL,  
      y = "Mean Standardized CPUE",  
      fill = NULL
    ) +
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
    annotate("text", 
             x = 9.6,  
             y = max(df_mean_filtered$mean_est) - 0.23,  
             label = paste("Severely Dry Mean PDSI = ", round(mean_PDSI_dry, 1)),
             size = 4.5, hjust = 0, vjust = 0.5, color = "black") +
    annotate("text", 
             x = 9.6,  
             y = max(df_mean_filtered$mean_est) - 0.37,  
             label = paste("Severely Wet Mean PDSI = ", round(mean_PDSI_wet, 1)),
             size = 4.5, hjust = 0, vjust = 0.5, color = "black") 
  
  # Save the plot
  ggsave("species_bay_PDSI_plot_test.png", plot, dpi = 200, bg = "white", width = 2500, height = 1500, units = "px")
  
  # Return the plot
  return(plot)
}

# Example call to the function










































































































































































































































































































