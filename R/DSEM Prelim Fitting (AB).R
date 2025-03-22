
library(dsem)
library(ggplot2)
library(ggpubr)
library(ggraph)
library(phylopath)
library(dplyr)
library(ggdag)
library(readxl)
library(ggraph)
library(qgraph)
library(DHARMa)


###################################################################


# First System and Major Bay: Sciaenid Trophic System and Aransas Bay 

# get data into a time series object from all numeric columns except YEAR
AB_Sciaenid_TS <- ts(
  AransasBay_Sciaenid_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Sciaenid_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Sciaenid_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the sciaenid system for aransas bay for zero and 1 year lags

# full bottom up
semAB_Sciaenid_fullbottomup <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Prey predicting Sciaenid
Atlanticcroaker_MesquiteBay -> RedDrum_MesquiteBay, 0, b1
Atlanticcroaker_AransasBay -> RedDrum_AransasBay, 0, b2
Atlanticcroaker_CopanoBay -> RedDrum_CopanoBay, 0, b3

BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, b4
BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, b5
BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, b6

Atlanticcroaker_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b7
Atlanticcroaker_AransasBay -> SpottedSeatrout_AransasBay, 0, b8
Atlanticcroaker_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b9

BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b10
BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, b11
BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

fit_semAB_Sciaenid_fullbottomup = dsem(sem = semAB_Sciaenid_fullbottomup,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_fullbottomup)

# extract the coefficients
get_part_Sciaenid = function(x){
  vars = c("Salinity","SpottedSeatrout","PDSI","Atlanticcroaker","RedDrum", "BlueCrabSmall")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# function to make path diagram
plot_qgraph_panel <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram") {
  
  # Extract the column names (variable names) from the time series object
  abbrev_names <- colnames(data_ts)
  
  # Check if the time series object has column names
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  
  # Set plot margins (bottom, left, top, right) to add cushion space
  par(mar = c(25, 25, 25, 25))  
  
  # Plot the graph
  qgraph(coef_matrix,         
         layout = "groups",              # Layout style
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Set positive relationships to blue
         negCol = "red3",                 # Optionally set negative relationships to red
         labels = abbrev_names,          # Use the automatically extracted variable names
         title = plot_title,             # Use the specified plot title
         title.cex = 1,                  # Increase font size of the title
         label.cex = .8,                  # Increase font size of the variable names
         vsize = 10,
         label.scale = FALSE,
         shape = "ellipse")
}


# get the coefficient into a matrix for the path diagrams
coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_fullbottomup, lag=0))
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_fullbottomup, lag=1))

# create an empty data frame with specified column names
df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          Atlanticcroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrabSmall = numeric(),
                          stringsAsFactors = FALSE)

png("~/Desktop/AB_Sciaenid.png", width = 1200, height = 2000, res = 200)
par(mfrow = c(2, 1))
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag$coef,  plot_title = "Sciaenid-AB-No Lag")
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag$coef, plot_title = "Sciaenid-AB-1 Year Lag")
par(mfrow = c(1, 1))
dev.off()

# no DD bottom up model
semAB_Sciaenid_bottomup_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Prey predicting Sciaenid
Atlanticcroaker_MesquiteBay -> RedDrum_MesquiteBay, 0, b1
Atlanticcroaker_AransasBay -> RedDrum_AransasBay, 0, b2
Atlanticcroaker_CopanoBay -> RedDrum_CopanoBay, 0, b3

BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, b4
BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, b5
BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, b6

Atlanticcroaker_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b7
Atlanticcroaker_AransasBay -> SpottedSeatrout_AransasBay, 0, b8
Atlanticcroaker_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b9

BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, b10
BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, b11
BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

"

fit_semAB_Sciaenid_bottomup_noDD = dsem(sem = semAB_Sciaenid_bottomup_noDD,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_bottomup_noDD)

# abiotic only model
semAB_Sciaenid_abiotic <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12


# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12
"

fit_semAB_Sciaenid_abiotic = dsem(sem = semAB_Sciaenid_abiotic,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_abiotic)

# Full top down model
semAB_Sciaenid_fulltopdown <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Sciaenid predicting prey
RedDrum_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b1
RedDrum_AransasBay -> Atlanticcroaker_AransasBay, 0, b2
RedDrum_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b3

RedDrum_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b4
RedDrum_AransasBay -> BlueCrabSmall_AransasBay, 0, b5
RedDrum_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b6

SpottedSeatrout_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b7
SpottedSeatrout_AransasBay -> Atlanticcroaker_AransasBay, 0, b8
SpottedSeatrout_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b9

SpottedSeatrout_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b10
SpottedSeatrout_AransasBay -> BlueCrabSmall_AransasBay, 0, b11
SpottedSeatrout_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Sciaenid_fulltopdown = dsem(sem = semAB_Sciaenid_fulltopdown,
                          tsdata = AB_Sciaenid_TS,
                          control = dsem_control(
                            quiet = TRUE))

summary(fit_semAB_Sciaenid_fulltopdown)



# top down no DD model
semAB_Sciaenid_topdown_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Sciaenid predicting prey
RedDrum_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b1
RedDrum_AransasBay -> Atlanticcroaker_AransasBay, 0, b2
RedDrum_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b3

RedDrum_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b4
RedDrum_AransasBay -> BlueCrabSmall_AransasBay, 0, b5
RedDrum_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b6

SpottedSeatrout_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, b7
SpottedSeatrout_AransasBay -> Atlanticcroaker_AransasBay, 0, b8
SpottedSeatrout_CopanoBay -> Atlanticcroaker_CopanoBay, 0, b9

SpottedSeatrout_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, b10
SpottedSeatrout_AransasBay -> BlueCrabSmall_AransasBay, 0, b11
SpottedSeatrout_CopanoBay -> BlueCrabSmall_CopanoBay, 0, b12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Sciaenid_topdown_noDD  = dsem(sem = semAB_Sciaenid_topdown_noDD,
                               tsdata = AB_Sciaenid_TS,
                               control = dsem_control(
                                 quiet = TRUE))
summary(fit_semAB_Sciaenid_topdown_noDD)


# No Prey Model
semAB_Sciaenid_notrophics <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> Atlanticcroaker_MesquiteBay, 0, p4
PDSI -> Atlanticcroaker_AransasBay, 0, p5
PDSI -> Atlanticcroaker_CopanoBay, 0, p6

PDSI -> BlueCrabSmall_MesquiteBay, 0, p7
PDSI -> BlueCrabSmall_AransasBay, 0, p8
PDSI -> BlueCrabSmall_CopanoBay, 0, p9

PDSI -> RedDrum_MesquiteBay, 0, p10
PDSI -> RedDrum_AransasBay, 0, p11
PDSI -> RedDrum_CopanoBay, 0, p12

PDSI -> SpottedSeatrout_MesquiteBay, 0, p13
PDSI -> SpottedSeatrout_AransasBay, 0, p14
PDSI -> SpottedSeatrout_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> Atlanticcroaker_MesquiteBay, 0, s1
Salinity_AransasBay -> Atlanticcroaker_AransasBay, 0, s2
Salinity_CopanoBay -> Atlanticcroaker_CopanoBay, 0, s3

Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, s4
Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s5
Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, s6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, s7
Salinity_AransasBay -> RedDrum_AransasBay, 0, s8
Salinity_CopanoBay -> RedDrum_CopanoBay, 0, s9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, s10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, s11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, s12

# Lagged relationships for Sciaenid (Lag = 1)
PDSI -> RedDrum_MesquiteBay, 1, l1
PDSI -> RedDrum_AransasBay, 1, l2
PDSI -> RedDrum_CopanoBay, 1, l3

PDSI -> SpottedSeatrout_MesquiteBay, 1, l4
PDSI -> SpottedSeatrout_AransasBay, 1, l5
PDSI -> SpottedSeatrout_CopanoBay, 1, l6

Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, l7
Salinity_AransasBay -> RedDrum_AransasBay, 1, l8
Salinity_CopanoBay -> RedDrum_CopanoBay, 1, l9

Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l10
Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, l11
Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l12

RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, l25
RedDrum_AransasBay -> RedDrum_AransasBay, 1, l26
RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, l27

SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, l28
SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, l29
SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, l30
"

# Fit the No Prey  model
fit_semAB_Sciaenid_notrophics = dsem(sem = semAB_Sciaenid_notrophics,
                              tsdata = AB_Sciaenid_TS,
                              control = dsem_control(
                                quiet = TRUE))

summary(fit_semAB_Sciaenid_notrophics)

# Compare AICs of all six models
AIC(fit_semAB_Sciaenid_abiotic)
AIC(fit_semAB_Sciaenid_notrophics)
AIC(fit_semAB_Sciaenid_bottomup_noDD)
AIC(fit_semAB_Sciaenid_topdown_noDD)
AIC(fit_semAB_Sciaenid_fullbottomup)
AIC(fit_semAB_Sciaenid_fulltopdown)


# function extract the coefficients
get_part_Sciaenid = function(x){
  vars = c("Salinity","SpottedSeatrout","PDSI","Atlanticcroaker","RedDrum", "BlueCrabSmall")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# function to make path diagram
plot_qgraph_panel <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram") {
  
  # Extract the column names (variable names) from the time series object
  abbrev_names <- colnames(data_ts)
  
  # Check if the time series object has column names
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  
  # Set plot margins (bottom, left, top, right) to add cushion space
  par(mar = c(25, 25, 25, 25))  
  
  # Plot the graph
  qgraph(coef_matrix,         
         layout = "groups",              # Layout style
         edge.labels = TRUE,             # Display edge labels
         posCol = "navy",                # Set positive relationships to blue
         negCol = "red3",                 # Optionally set negative relationships to red
         labels = abbrev_names,          # Use the automatically extracted variable names
         title = plot_title,             # Use the specified plot title
         title.cex = 1,                  # Increase font size of the title
         label.cex = .8,                  # Increase font size of the variable names
         vsize = 10,
         label.scale = FALSE,
         shape = "ellipse")
}


# get the coefficient into a matrix for the path diagrams, for the lowest AIC Model
coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=0))
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fit_semAB_Sciaenid_notrophics, lag=1))

# create an empty data frame with specified column names
df_Sciaenid <- data.frame(Salinity = numeric(),
                          SpottedSeatrout = numeric(),
                          PDSI = numeric(),
                          Atlanticcroaker = numeric(),
                          RedDrum = numeric(),
                          BlueCrabSmall = numeric(),
                          stringsAsFactors = FALSE)

png("~/Desktop/AB_Sciaenid.png", width = 1200, height = 2000, res = 200)
par(mfrow = c(2, 1))
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag$coef,  plot_title = "Sciaenid-AB-No Lag")
plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag$coef, plot_title = "Sciaenid-AB-1 Year Lag")
par(mfrow = c(1, 1))
dev.off()




###################################################################




# Repeating the above but for the...

# Second System and Major Bay: Keystone Predator Trophic System and Aransas Bay 

#get data into a time series object from all numeric columns except YEAR
AB_Pred_TS <- ts(
  AransasBay_Pred_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Pred_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Pred_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# Full bottom up model
semAB_Pred_fullbottomup <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Prey predicting predators
AllMullet_MesquiteBay -> BullShark_MesquiteBay, 0, b1
AllMullet_AransasBay -> BullShark_AransasBay, 0, b2
AllMullet_CopanoBay -> BullShark_CopanoBay, 0, b3

AllMenhaden_MesquiteBay -> BullShark_MesquiteBay, 0, b4
AllMenhaden_AransasBay -> BullShark_AransasBay, 0, b5
AllMenhaden_CopanoBay -> BullShark_CopanoBay, 0, b6

AllMullet_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b7
AllMullet_AransasBay -> AlligatorGar_AransasBay, 0, b8
AllMullet_CopanoBay -> AlligatorGar_CopanoBay, 0, b9

AllMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b10
AllMenhaden_AransasBay -> AlligatorGar_AransasBay, 0, b11
AllMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 0, b12

# Lagged relationships (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

AllMullet_MesquiteBay -> AllMullet_MesquiteBay, 1, l13
AllMullet_AransasBay -> AllMullet_AransasBay, 1, l14
AllMullet_CopanoBay -> AllMullet_CopanoBay, 1, l15

AllMenhaden_MesquiteBay -> AllMenhaden_MesquiteBay, 1, l16
AllMenhaden_AransasBay -> AllMenhaden_AransasBay, 1, l17
AllMenhaden_CopanoBay -> AllMenhaden_CopanoBay, 1, l18

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30
"

fit_semAB_Pred_fullbottomup = dsem(sem = semAB_Pred_fullbottomup,
                      tsdata = AB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_semAB_Pred_fullbottomup)


# Abiotic Only Model
semAB_Pred_abiotic <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12
"

fit_semAB_Pred_abiotic = dsem(sem = semAB_Pred_abiotic,
                                  tsdata = AB_Pred_TS,
                                  control = dsem_control(
                                    quiet = TRUE))
summary(fit_semAB_Pred_abiotic)


# No Density Dependence bottom up Model
semAB_Pred_bottomupnoDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Prey predicting predators
AllMullet_MesquiteBay -> BullShark_MesquiteBay, 0, b1
AllMullet_AransasBay -> BullShark_AransasBay, 0, b2
AllMullet_CopanoBay -> BullShark_CopanoBay, 0, b3

AllMenhaden_MesquiteBay -> BullShark_MesquiteBay, 0, b4
AllMenhaden_AransasBay -> BullShark_AransasBay, 0, b5
AllMenhaden_CopanoBay -> BullShark_CopanoBay, 0, b6

AllMullet_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b7
AllMullet_AransasBay -> AlligatorGar_AransasBay, 0, b8
AllMullet_CopanoBay -> AlligatorGar_CopanoBay, 0, b9

AllMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 0, b10
AllMenhaden_AransasBay -> AlligatorGar_AransasBay, 0, b11
AllMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888
"

# Fit the No Density Dependence bottom up  model
fit_semAB_Pred_bottomupnoDD = dsem(sem = semAB_Pred_bottomupnoDD,
                        tsdata = AB_Pred_TS,
                        control = dsem_control(
                          quiet = TRUE))
summary(fit_semAB_Pred_bottomupnoDD)


# No trophics Model
semAB_Pred_notrophics <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30
"

# Fit the No trophics SEM model
fit_semAB_Pred_notrophics= dsem(sem = semAB_Pred_notrophics,
                          tsdata = AB_Pred_TS,
                          control = dsem_control(
                            quiet = TRUE))
summary(fit_semAB_Pred_notrophics)


# Full model for top down
semAB_Pred_fulltopdown <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Predators predicting prey
BullShark_MesquiteBay -> AllMullet_MesquiteBay, 0, b1
BullShark_AransasBay -> AllMullet_AransasBay, 0, b2
BullShark_CopanoBay -> AllMullet_CopanoBay, 0, b3

BullShark_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b4
BullShark_AransasBay -> AllMenhaden_AransasBay, 0, b5
BullShark_CopanoBay -> AllMenhaden_CopanoBay, 0, b6

AlligatorGar_MesquiteBay -> AllMullet_MesquiteBay, 0, b7
AlligatorGar_AransasBay -> AllMullet_AransasBay, 0, b8
AlligatorGar_CopanoBay -> AllMullet_CopanoBay, 0, b9

AlligatorGar_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b10
AlligatorGar_AransasBay -> AllMenhaden_AransasBay, 0, b11
AlligatorGar_CopanoBay -> AllMenhaden_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 1)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888

BullShark_MesquiteBay -> BullShark_MesquiteBay, 1, l25
BullShark_AransasBay -> BullShark_AransasBay, 1, l26
BullShark_CopanoBay -> BullShark_CopanoBay, 1, l27

AlligatorGar_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l28
AlligatorGar_AransasBay -> AlligatorGar_AransasBay, 1, l29
AlligatorGar_CopanoBay -> AlligatorGar_CopanoBay, 1, l30

AllMullet_MesquiteBay -> AllMullet_MesquiteBay, 1, l13
AllMullet_AransasBay -> AllMullet_AransasBay, 1, l14
AllMullet_CopanoBay -> AllMullet_CopanoBay, 1, l15

AllMenhaden_MesquiteBay -> AllMenhaden_MesquiteBay, 1, l16
AllMenhaden_AransasBay -> AllMenhaden_AransasBay, 1, l17
AllMenhaden_CopanoBay -> AllMenhaden_CopanoBay, 1, l18
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Pred_fulltopdown = dsem(sem = semAB_Pred_fulltopdown,
                      tsdata = AB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_semAB_Pred_fulltopdown)

# top down no density dependence model
semAB_Pred_topdown_noDD <- "
PDSI -> Salinity_MesquiteBay, 0, p1
PDSI -> Salinity_AransasBay, 0, p2
PDSI -> Salinity_CopanoBay, 0, p3

PDSI -> AllMullet_MesquiteBay, 0, p4
PDSI -> AllMullet_AransasBay, 0, p5
PDSI -> AllMullet_CopanoBay, 0, p6

PDSI -> AllMenhaden_MesquiteBay, 0, p7
PDSI -> AllMenhaden_AransasBay, 0, p8
PDSI -> AllMenhaden_CopanoBay, 0, p9

PDSI -> BullShark_MesquiteBay, 0, p10
PDSI -> BullShark_AransasBay, 0, p11
PDSI -> BullShark_CopanoBay, 0, p12

PDSI -> AlligatorGar_MesquiteBay, 0, p13
PDSI -> AlligatorGar_AransasBay, 0, p14
PDSI -> AlligatorGar_CopanoBay, 0, p15

# Salinity as a predictor
Salinity_MesquiteBay -> AllMullet_MesquiteBay, 0, s1
Salinity_AransasBay -> AllMullet_AransasBay, 0, s2
Salinity_CopanoBay -> AllMullet_CopanoBay, 0, s3

Salinity_MesquiteBay -> AllMenhaden_MesquiteBay, 0, s4
Salinity_AransasBay -> AllMenhaden_AransasBay, 0, s5
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 0, s6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, s7
Salinity_AransasBay -> BullShark_AransasBay, 0, s8
Salinity_CopanoBay -> BullShark_CopanoBay, 0, s9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, s10
Salinity_AransasBay -> AlligatorGar_AransasBay, 0, s11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, s12

# Predators predicting prey
BullShark_MesquiteBay -> AllMullet_MesquiteBay, 0, b1
BullShark_AransasBay -> AllMullet_AransasBay, 0, b2
BullShark_CopanoBay -> AllMullet_CopanoBay, 0, b3

BullShark_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b4
BullShark_AransasBay -> AllMenhaden_AransasBay, 0, b5
BullShark_CopanoBay -> AllMenhaden_CopanoBay, 0, b6

AlligatorGar_MesquiteBay -> AllMullet_MesquiteBay, 0, b7
AlligatorGar_AransasBay -> AllMullet_AransasBay, 0, b8
AlligatorGar_CopanoBay -> AllMullet_CopanoBay, 0, b9

AlligatorGar_MesquiteBay -> AllMenhaden_MesquiteBay, 0, b10
AlligatorGar_AransasBay -> AllMenhaden_AransasBay, 0, b11
AlligatorGar_CopanoBay -> AllMenhaden_CopanoBay, 0, b12

# Lagged relationships for predators (Lag = 2)
PDSI -> BullShark_MesquiteBay, 1, l1
PDSI -> BullShark_AransasBay, 1, l2
PDSI -> BullShark_CopanoBay, 1, l3

PDSI -> AlligatorGar_MesquiteBay, 1, l4
PDSI -> AlligatorGar_AransasBay, 1, l5
PDSI -> AlligatorGar_CopanoBay, 1, l6

Salinity_MesquiteBay -> BullShark_MesquiteBay, 1, l7
Salinity_AransasBay -> BullShark_AransasBay, 1, l8
Salinity_CopanoBay -> BullShark_CopanoBay, 1, l9

Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 1, l10
Salinity_AransasBay -> AlligatorGar_AransasBay, 1, l11
Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, l12

PDSI -> AllMullet_MesquiteBay, 1, l133
PDSI -> AllMullet_AransasBay, 1, l144
PDSI -> AllMullet_CopanoBay, 1, l155

PDSI -> AllMenhaden_MesquiteBay, 1, l166
PDSI -> AllMenhaden_AransasBay, 1, l177
PDSI -> AllMenhaden_CopanoBay, 1, l188

Salinity_MesquiteBay -> AllMullet_MesquiteBay, 1, l1333
Salinity_AransasBay -> AllMullet_AransasBay, 1, l1444
Salinity_CopanoBay -> AllMullet_CopanoBay, 1, l1555

Salinity_MesquiteBay ->  AllMenhaden_MesquiteBay, 1, l1666
Salinity_AransasBay -> AllMenhaden_AransasBay, 1, l7777
Salinity_CopanoBay -> AllMenhaden_CopanoBay, 1, l1888
"

# Fit the model for prey prediction in Aransas Bay
fit_semAB_Pred_topdown_noDD = dsem(sem = semAB_Pred_topdown_noDD ,
                           tsdata = AB_Pred_TS,
                           control = dsem_control(
                             quiet = TRUE))
summary(fit_semAB_Pred_topdown_noDD)

AIC(fit_semAB_Pred_abiotic)
AIC(fit_semAB_Pred_notrophics)
AIC(fit_semAB_Pred_topdown_noDD)
AIC(fit_semAB_Pred_bottomupnoDD)
AIC(fit_semAB_Pred_fulltopdown)
AIC(fit_semAB_Pred_fullbottomup)



get_part_pred = function(x){
  vars = c("Salinity","AllMullet","BullShark","AlligatorGar","PDSI", "AllMenhaden")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

coef_matrix_GB_Pred_NoLag <- get_part_pred(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=0))
coef_matrix_GB_Pred_YesLag <- get_part_pred(as_fitted_DAG(fit_semAB_Pred_fulltopdown, lag=1))

df_Pred_GB <- data.frame(Salinity = numeric(),
                         AllMullet = numeric(),
                         BullShark = numeric(),
                         AlligatorGar = numeric(),
                         PDSI = numeric(),
                         AllMenhaden = numeric(),
                         stringsAsFactors = FALSE)

par(mfrow = c(1, 2))
plot1 <- plot_qgraph_panel(df_Pred_GB, coef_matrix_GB_Pred_NoLag$coef,  plot_title = "Keystone Predator-GB-No Lag")
plot2 <- plot_qgraph_panel(df_Pred_GB, coef_matrix_GB_Pred_YesLag$coef, plot_title = "Keystone Predator-GB-1 Year Lag")
par(mfrow = c(1, 1))




###################################################################




#Testing my voice-to-text app. 