

library(dplyr)
library(tidyr)

# merging temp and salinity data from three surveys
AransasBay_Abiotics <- bind_rows(
  AransasBay_BS_shrimp %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  AransasBay_GN_cpue %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE),
  AransasBay_BT %>%
    select(MINOR_AREA_CODE, STATION_CODE, START_TEMPERATURE_NUM, START_SALINITY_NUM, DATE)
)

# function to get data in mean annual values per minor bay in wide format for dsem
processforDSEM <- function(data, species, minor_area_mapping, date_column, result_df) {
  # Step 1: Assign values to the new column MINOR_BAY based on MINOR_AREA_CODE
  data <- data %>%
    mutate(MINOR_BAY = case_when(
      MINOR_AREA_CODE %in% minor_area_mapping$Copano_Bay ~ "Copano Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$Aransas_Bay ~ "Aransas Bay",
      MINOR_AREA_CODE %in% minor_area_mapping$Mesquite_Bay ~ "Mesquite Bay",
      TRUE ~ NA_character_
    ))
  
  # Step 2: Extract YEAR from the specified date column
  data$YEAR <- substr(data[[date_column]], 1, 4)
  
  # Step 3: Calculate the mean value of each species per YEAR and MINOR_BAY
  wide_data <- data %>%
    group_by(YEAR, MINOR_BAY) %>%
    summarize(across(all_of(species), ~ mean(.x, na.rm = TRUE), .names = "{.col}"), .groups = "drop") %>%
    pivot_wider(
      names_from = MINOR_BAY,
      values_from = all_of(species),
      names_glue = "{.value}_{MINOR_BAY}"  # Ensure consistent naming: species_minor_bay
    )
  
  # Step 4: Merge the output with the existing result dataframe
  result_df <- result_df %>%
    full_join(wide_data, by = "YEAR")
  
  return(result_df)
}

# define mapping of MINOR_AREA_CODEs to MINOR_BAY
minor_area_mapping <- list(
  Copano_Bay = c(270, 317, 293, 120, 240, 310),
  Aransas_Bay = c(227, 20, 13, 280, 285, 95, 94),
  Mesquite_Bay = c(143, 70, 315, 250, 90)
)

# make initial df
AransasBay_Sciaenid_Wide <- data.frame(YEAR = character())

# use function to get values for 5 fish species from gill net surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_GN_cpue,
  species = c("HardheadCatfish", "GafftopsailCatfish", "RedDrum", "BlackDrum", "SpottedSeatrout"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# use function to get values for brown shrimp from bag seine surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_BS_shrimp,
  species = c("BrownShrimp"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# use function to get values for juvenile blue crabs from bag seine surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_BS_bluecrabs,
  species = c("BlueCrabSmall"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# use function to get values for temp and salinity from all gill net, bag seine and bottom trawl surveys
AransasBay_Sciaenid_Wide <- processforDSEM(
  data = AransasBay_Abiotics,
  species = c("START_TEMPERATURE_NUM", "START_SALINITY_NUM"),
  minor_area_mapping = minor_area_mapping,
  date_column = "DATE",
  result_df = AransasBay_Sciaenid_Wide
)

# the df now has all the needed data (I think?). we still need to transform/standardize the values, but first lets clean up the df (remove and rename some columns)
AransasBay_Sciaenid_Wide <- AransasBay_Sciaenid_Wide %>%
  select(-BrownShrimp_NA, -START_TEMPERATURE_NUM_NA, -START_SALINITY_NUM_NA, -BlueCrabSmall_NA)

AransasBay_Sciaenid_Wide <- AransasBay_Sciaenid_Wide %>%
  rename_with(
    .fn = ~ gsub("^START_TEMPERATURE_NUM_", "Temperature_", .),
    .cols = starts_with("START_TEMPERATURE_NUM_")
  ) %>%
  rename_with(
    .fn = ~ gsub("^START_SALINITY_NUM_", "Salinity_", .),
    .cols = starts_with("START_SALINITY_NUM_")
  )

# quick inspection of some temporal trends, using catfish because I have an idea of what these trends should look like based on zach's recent work 
library(ggplot2)

time_series_data <- AransasBay_Sciaenid_Wide %>%
  select(YEAR, starts_with("GafftopsailCatfish"), starts_with("HardheadCatfish")) %>%
  pivot_longer(
    cols = -YEAR,
    names_to = c("Species", "Bay"),
    names_sep = "_",  # Assumes the column names are in "Species_Bay" format
    values_to = "Mean_CPUE"
  )

ggplot(time_series_data, aes(x = YEAR, y = Mean_CPUE, color = Bay, group = interaction(Species, Bay))) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ Species, scales = "free_y") +
  labs(
    title = "Time Series of Gafftopsail Catfish and Hardhead Catfish CPUE",
    x = "Year",
    y = "Mean CPUE",
    color = "Bay"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# remove spaces from column names
colnames(AransasBay_Sciaenid_Wide) <- gsub(" ", "", colnames(AransasBay_Sciaenid_Wide))


# okay now lets transform the values and create a new df
AransasBay_Sciaenid_Wide_Trans <- AransasBay_Sciaenid_Wide %>%
  # Step 1: Apply log transformation
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ log(. + 1)  # Log transformation with a small shift to avoid log(0)
  )) %>%
  # Step 2: Apply standardization (mean = 0, SD = 1)
  mutate(across(
    .cols = where(is.numeric) & !contains("YEAR"),
    .fns = ~ (. - mean(.)) / sd(.)
  ))

#now data is ready for dsem. lets export into an excel file 
