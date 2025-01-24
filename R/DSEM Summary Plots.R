
library(ggplot2)
library(dplyr)

###################################################################

#### here, i want to take all of the dsem coefficients from the annual time step dsems and ...
### compare them to the seasonal time step dsems
#### i do that by presenting box plots of mean coefficient values (with SE) for each unique predictor-response combo
#### the means come from four values: two bay systems, with and without lags within each major bay system...
### although i acknowledge grouping lagged and non lagged values can be problematic 
#### this is repeated for each trophic system (sciaenid and keystone predator)
#### the hypothesis here is that annual time steps will produce different (in terms of direction and magnitude)...
### relationship than seasonal time steps

### all code/analyses from lists produced from the two other scripts: "DSEM Prelim Fitting (Annual)" and "DSEM Prelim Fitting (Seasonal)"

#### predator model function to get df of dsem coefficients
convert_to_df_pred <- function(coef_matrix) {
  if (!is.list(coef_matrix)) {
    stop("Input must be a list.")
  }
  matrix <- do.call(cbind, coef_matrix)
  rownames(matrix) <- c("Salinity", "StripedMullet", "BullShark", "AlligatorGar", "Temperature", "GulfMenhaden")
  as.data.frame(as.table(matrix))
}



#### apply to annual time step (keystone predator system)
df1_pred <- convert_to_df_pred(coef_matrix_AB_Pred_NoLag)
df2_pred <- convert_to_df_pred(coef_matrix_AB_Pred_YesLag)
df3_pred <- convert_to_df_pred(coef_matrix_GB_Pred_NoLag)
df4_pred <- convert_to_df_pred(coef_matrix_GB_Pred_YesLag)

# combine all data frames
all_df_pred <- rbind(df1_pred, df2_pred, df3_pred, df4_pred)
colnames(all_df_pred) <- c("Predictor", "Response", "Value")

# aggregate to calculate the mean and standard error
predplot_df <- all_df_pred %>%
  group_by(Predictor, Response) %>%
  summarize(
    Mean = mean(Value),
    StdError = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# filter out rows where Response is "Temperature" or "Salinity"
predplot_df_filt <-predplot_df %>%
  filter(!Response %in% c("Temperature", "Salinity")) %>%
  filter(!Predictor %in% c("BullShark", "AlligatorGar"))

# reorder the Predictor column
predplot_df_filt$Predictor <- factor(
  predplot_df_filt$Predictor,
  levels = c("Salinity", "Temperature", "StripedMullet", "GulfMenhaden")
)

# remove rows for relationships that do not exist
predplot_df_filt <- predplot_df_filt %>%
  filter(Mean != 0)

# create the bar graph with error bars and updated styling
AnnualPredSummaryPlot<-ggplot(predplot_df_filt, aes(x = Predictor, y = Mean, fill = Response)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),   color = "grey55") + 
  geom_errorbar(
    aes(ymin = Mean - StdError, ymax = Mean + StdError),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) + 
  labs(
    x = "Predictor",
    y = "Coefficient ",
    fill = "Response",
    title = "Keystone Predator System - Annual Timestep"
  ) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0),
        panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 2), "lines") 
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_discrete(
    expand = expansion(mult = c(0, 0)) 
  ) +
  geom_vline(
    xintercept = seq(0.5, length(unique(predplot_df_filt$Predictor)), by = 1),
    color = "black",
    linetype = "solid"
  ) 



###################################################################



#### sciaenid model function to get df of dsem coefficients
convert_to_df_sciaenid<- function(coef_matrix) {
  if (!is.list(coef_matrix)) {
    stop("Input must be a list.")
  }
  matrix <- do.call(cbind, coef_matrix)
  rownames(matrix) <- c("Salinity", "BrownShrimp", "SpottedSeatrout","Temperature", "SouthernFlounder", "RedDrum",  "BlueCrabSmall")
  as.data.frame(as.table(matrix))
}



#### apply to annual time step (sciaenid system)
df1_sciaenid <- convert_to_df_sciaenid(coef_matrix_AB_Sciaenid_NoLag)
df2_sciaenid <- convert_to_df_sciaenid(coef_matrix_AB_Sciaenid_YesLag)
df3_sciaenid <- convert_to_df_sciaenid(coef_matrix_GB_Sciaenid_NoLag)
df4_sciaenid <- convert_to_df_sciaenid(coef_matrix_GB_Sciaenid_YesLag)

# combine all data frames
all_df_sciaenid <- rbind(df1_sciaenid, df2_sciaenid, df3_sciaenid, df4_sciaenid)
colnames(all_df_sciaenid) <- c("Predictor", "Response", "Value")

# aggregate to calculate the mean and standard error
sciaenidplot_df <- all_df_sciaenid %>%
  group_by(Predictor, Response) %>%
  summarize(
    Mean = mean(Value),
    StdError = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# filter out rows where Response is "Temperature" or "Salinity"
sciaenidplot_df_filt <-sciaenidplot_df %>%
  filter(!Response %in% c("Temperature", "Salinity")) %>%
  filter(!Predictor %in% c("SouthernFlounder", "RedDrum", "SpottedSeatrout"))

# reorder the Predictor column
sciaenidplot_df_filt$Predictor <- factor(
  sciaenidplot_df_filt$Predictor,
  levels = c("Salinity", "Temperature", "BlueCrabSmall", "BrownShrimp")
)

sciaenidplot_df_filt <- sciaenidplot_df_filt %>%
  filter(Mean != 0)

# create the bar graph with error bars and updated styling
AnnualSciaenidSummaryPlot<-ggplot(sciaenidplot_df_filt, aes(x = Predictor, y = Mean, fill = Response)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),   color = "grey55") + 
  geom_errorbar(
    aes(ymin = Mean - StdError, ymax = Mean + StdError),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) + 
  labs(
    x = "Predictor",
    y = "Coefficient ",
    fill = "Response",
    title = "Sciaenid System - Annual Timestep"
  ) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0),
        panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 2), "lines") 
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_discrete(
    expand = expansion(mult = c(0, 0)) 
  ) +
  geom_vline(
    xintercept = seq(0.5, length(unique(sciaenidplot_df_filt$Predictor)), by = 1),
    color = "black",
    linetype = "solid"
  ) 



###########################################################



#### apply to seasonal time step (sciaenid system)
df1_sciaenid_Season <- convert_to_df_sciaenid(coef_matrix_AB_Sciaenid_NoLag_Season)
df2_sciaenid_Season <- convert_to_df_sciaenid(coef_matrix_AB_Sciaenid_YesLag_Season)
df3_sciaenid_Season <- convert_to_df_sciaenid(coef_matrix_GB_Sciaenid_NoLag_Season)
df4_sciaenid_Season <- convert_to_df_sciaenid(coef_matrix_GB_Sciaenid_YesLag_Season)

# combine all data frames
all_df_sciaenid_Season <- rbind(df1_sciaenid_Season, df2_sciaenid_Season, df3_sciaenid_Season, df4_sciaenid_Season)
colnames(all_df_sciaenid_Season) <- c("Predictor", "Response", "Value")

# aggregate to calculate the mean and standard error
sciaenidplot_df_Season <- all_df_sciaenid_Season %>%
  group_by(Predictor, Response) %>%
  summarize(
    Mean = mean(Value),
    StdError = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# filter out rows where Response is "Temperature" or "Salinity"
sciaenidplot_df_filt_Season <-sciaenidplot_df_Season %>%
  filter(!Response %in% c("Temperature", "Salinity")) %>%
  filter(!Predictor %in% c("SouthernFlounder", "RedDrum", "SpottedSeatrout"))

# reorder the Predictor column
sciaenidplot_df_filt_Season$Predictor <- factor(
  sciaenidplot_df_filt_Season$Predictor,
  levels = c("Salinity", "Temperature", "BlueCrabSmall", "BrownShrimp")
)

sciaenidplot_df_filt_Season<- sciaenidplot_df_filt_Season %>%
  filter(Mean != 0)

# create the bar graph with error bars and updated styling
SeasonalSciaenidSummaryPlot<-ggplot(sciaenidplot_df_filt_Season, aes(x = Predictor, y = Mean, fill = Response)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),   color = "grey55") +
  geom_errorbar(
    aes(ymin = Mean - StdError, ymax = Mean + StdError),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) + 
  labs(
    x = "Predictor",
    y = "Coefficient ",
    fill = "Response",
    title = "Sciaenid System - Seasonal Timestep"
  ) +
  theme_bw() +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(angle = 0),
        panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 2), "lines") 
  ) +
  scale_fill_brewer(palette = "Set3") + 
  scale_x_discrete(
    expand = expansion(mult = c(0, 0))
  ) +
  geom_vline(
    xintercept = seq(0.5, length(unique(sciaenidplot_df_filt_Season$Predictor)), by = 1),
    color = "black",
    linetype = "solid"
  ) 



###########################################################



#### apply to seasonal time step (keystone predator system)
df1_pred_Season <- convert_to_df_pred(coef_matrix_AB_Pred_NoLag_Season)
df2_pred_Season <- convert_to_df_pred(coef_matrix_AB_Pred_YesLag_Season)
df3_pred_Season <- convert_to_df_pred(coef_matrix_GB_Pred_NoLag_Season)
df4_pred_Season <- convert_to_df_pred(coef_matrix_GB_Pred_YesLag_Season)

# combine all data frames
all_df_pred_Season <- rbind(df1_pred_Season, df2_pred_Season, df3_pred_Season, df4_pred_Season)
colnames(all_df_pred_Season) <- c("Predictor", "Response", "Value")

# aggregate to calculate the mean and standard error
predplot_df_Season <- all_df_pred_Season %>%
  group_by(Predictor, Response) %>%
  summarize(
    Mean = mean(Value),
    StdError = sd(Value) / sqrt(n()),
    .groups = "drop"
  )

# filter out rows where Response is "Temperature" or "Salinity"
predplot_df_filt_Season <-predplot_df_Season %>%
  filter(!Response %in% c("Temperature", "Salinity")) %>%
  filter(!Predictor %in% c("BullShark", "AlligatorGar"))

# reorder the Predictor column
predplot_df_filt_Season$Predictor <- factor(
  predplot_df_filt_Season$Predictor,
  levels = c("Salinity", "Temperature", "StripedMullet", "GulfMenhaden")
)

predplot_df_filt_Season <- predplot_df_filt_Season %>%
  filter(Mean != 0)

# create the bar graph with error bars and updated styling
SeasonalPredSummaryPlot<-ggplot(predplot_df_filt_Season, aes(x = Predictor, y = Mean, fill = Response)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9),   color = "grey55") + 
  geom_errorbar(
    aes(ymin = Mean - StdError, ymax = Mean + StdError),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) + 
  labs(
    x = "Predictor",
    y = "Coefficient ",
    fill = "Response",
    title = "Keystone Predator System - Seasonal Timestep"
  ) +
  theme_bw() +
  theme(text = element_text(size = 14),
    axis.text.x = element_text(angle = 0),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 2), "lines") 
  ) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_discrete(
    expand = expansion(mult = c(0, 0)) 
  ) +
  geom_vline(
    xintercept = seq(0.5, length(unique(predplot_df_filt_Season$Predictor)), by = 1),
    color = "black",
    linetype = "solid"
  ) 


#### combine and export plots
library(cowplot)

CombinedPredPlot <- plot_grid(
  SeasonalPredSummaryPlot, 
  AnnualPredSummaryPlot, 
  ncol = 1,  
  align = "v"  
)

ggsave(
  filename = "CombinedPredPlot.png",
  plot = CombinedPredPlot,
  width = 10,
  height = 12,
  dpi = 300
)

CombinedSciaenidPlot <- plot_grid(
  SeasonalSciaenidSummaryPlot, 
  AnnualSciaenidSummaryPlot, 
  ncol = 1, 
  align = "v"  
)

ggsave(
  filename = "CombinedSciaenidPlot.png",
  plot = CombinedSciaenidPlot,
  width = 10,
  height = 12,
  dpi = 300
)