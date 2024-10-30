# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggpubr)


#### ANNUAL MEANS 

# Step 1: Create a "Year" column in both dataframes by extracting it from the "Date" column
NERR_WQ_AB <- MANERR_WQ_AransasBay_AB %>%
  mutate(Year = as.numeric(format(as.Date(Date), "%Y")))

TPWD_Abiotics_AB  <- AransasBay_Abiotics %>%
  mutate(Year = as.numeric(format(as.Date(Date), "%Y")))

# Step 2: Calculate the mean Temp and Salinity per Year for each dataframe
mean_NERR_AB <- NERR_WQ_AB %>%
  group_by(Year) %>%
  summarise(mean_temp_NERR = mean(Temp, na.rm = TRUE),
            mean_salinity_NERR = mean(Salinity, na.rm = TRUE))

mean_TPWD_AB  <- TPWD_Abiotics_AB  %>%
  group_by(Year) %>%
  summarise(mean_temp_TPWD = mean(Temp, na.rm = TRUE),
            mean_salinity_TPWD = mean(Salinity, na.rm = TRUE))

# Step 3: Merge the two dataframes on the "Year" column
AnnualTemps <- merge(mean_TPWD_AB, mean_NERR_AB, by = "Year")

# Step 4: Filter the new dataframe to only include data from the year 2007 onwards
AnnualTemps <- AnnualTemps %>%
  filter(Year >= 2007)

# Step 5: Calculate the Pearson and Spearman correlation coefficients for Temp and Salinity
# Temperature correlations
temp_pearson <- cor(AnnualTemps$mean_temp_TPWD, AnnualTemps$mean_temp_NERR, use = "complete.obs", method = "pearson")
temp_spearman <- cor(AnnualTemps$mean_temp_TPWD, AnnualTemps$mean_temp_NERR, use = "complete.obs", method = "spearman")

# Salinity correlations
salinity_pearson <- cor(AnnualTemps$mean_salinity_TPWD, AnnualTemps$mean_salinity_NERR, use = "complete.obs", method = "pearson")
salinity_spearman <- cor(AnnualTemps$mean_salinity_TPWD, AnnualTemps$mean_salinity_NERR, use = "complete.obs", method = "spearman")


# Step 6: Plot the correlation between mean Temp values
ann_temp_plot <- ggplot(AnnualTemps, aes(x = mean_temp_AB, y = mean_temp_BS)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "grey55") +
  ggtitle(paste("Annual Temp Correlation (Pearson: ", round(temp_pearson, 2), ", Spearman: ", round(temp_spearman, 2), ")")) +
  xlab("MANERR Temp") +
  ylab("Bag Seine Temp") +
  theme_bw()
ann_temp_plot

# Step 7: Plot the correlation between mean Salinity values
ann_salinity_plot <- ggplot(AnnualTemps, aes(x = mean_salinity_AB, y = mean_salinity_BS)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "grey55") +
  ggtitle(paste("Annual Salinity Correlation (Pearson: ", round(salinity_pearson, 2), ", Spearman: ", round(salinity_spearman, 2), ")")) +
  xlab("MANERR Salinity") +
  ylab("Bag Seine Salinity") +
  theme_bw()

# Display the plots
print(temp_plot)
print(salinity_plot)

### MONTHLY MEANS

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Create a "Month.Year" column in both dataframes by extracting month and year from the "Date" column
NERR_WQ_AB <- MANERR_WQ_AransasBay_AB %>%
  mutate(Month.Year = format(as.Date(Date), "%m.%y"))

TPWD_Abiotics_AB  <- AransasBay_Abiotics %>%
  mutate(Month.Year = format(as.Date(Date), "%m.%y"))


# Step 2: Calculate the mean Temp and Salinity per Year for each dataframe
mean_NERR_AB <- NERR_WQ_AB %>%
  group_by(Month.Year) %>%
  summarise(mean_temp_NERR = mean(Temp, na.rm = TRUE),
            mean_salinity_NERR = mean(Salinity, na.rm = TRUE))

mean_TPWD_AB  <- TPWD_Abiotics_AB  %>%
  group_by(Month.Year) %>%
  summarise(mean_temp_TPWD = mean(Temp, na.rm = TRUE),
            mean_salinity_TPWD = mean(Salinity, na.rm = TRUE))

# Step 3: Merge the two dataframes on the "Month.Year" column
MonthlyTemps <- merge(mean_TPWD_AB, mean_NERR_AB, by = "Month.Year")

# Step 4: Calculate the Pearson and Spearman correlation coefficients for Temp and Salinity
# Temperature correlations
temp_pearson <- cor(MonthlyTemps$mean_temp_TPWD, MonthlyTemps$mean_temp_NERR, use = "complete.obs", method = "pearson")
temp_spearman <- cor(MonthlyTemps$mean_temp_TPWD, MonthlyTemps$mean_temp_NERR, use = "complete.obs", method = "spearman")

# Salinity correlations
salinity_pearson <- cor(MonthlyTemps$mean_salinity_AB, MonthlyTemps$mean_salinity_BS, use = "complete.obs", method = "pearson")
salinity_spearman <- cor(MonthlyTemps$mean_salinity_AB, MonthlyTemps$mean_salinity_BS, use = "complete.obs", method = "spearman")

# Output correlation coefficients
print(paste("Temperature Pearson Correlation Coefficient:", temp_pearson))
print(paste("Temperature Spearman Correlation Coefficient:", temp_spearman))
print(paste("Salinity Pearson Correlation Coefficient:", salinity_pearson))
print(paste("Salinity Spearman Correlation Coefficient:", salinity_spearman))

# Step 5: Plot the correlation between mean Temp values
monthly_temp_plot <- ggplot(MonthlyTemps, aes(x = mean_temp_AB, y = mean_temp_BS)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "grey55") +
  ggtitle(paste("Monthly Temp Correlation (Pearson: ", round(temp_pearson, 2), ", Spearman: ", round(temp_spearman, 2), ")")) +
  xlab("MANERR Temp") +
  ylab("Bag Seine Temp") +
  theme_bw()

# Step 6: Plot the correlation between mean Salinity values
monthly_salinity_plot <- ggplot(MonthlyTemps, aes(x = mean_salinity_AB, y = mean_salinity_BS)) +
  geom_point(color = "blue") +
  geom_smooth(method = "lm", se = FALSE, color = "grey55") +
  ggtitle(paste("Monthly Salinity Correlation (Pearson: ", round(salinity_pearson, 2), ", Spearman: ", round(salinity_spearman, 2), ")")) +
  xlab("MANERR Salinity") +
  ylab("Bag Seine Salinity") +
  theme_bw()

combined_plot <- ggarrange(ann_temp_plot, ann_salinity_plot, monthly_temp_plot, monthly_salinity_plot, ncol = 2, nrow = 2)
ggsave("ABcorrelations.png", plot = combined_plot, width = 12, height = 8)


