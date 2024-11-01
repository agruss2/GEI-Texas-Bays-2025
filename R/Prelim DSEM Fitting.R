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

# set up sem structure for a time lag of zero and one year
sem1= "
    BlueCrabSmall ->  RedDrum, 0, a
   Temp ->  RedDrum, 0, b
   Temp ->  BlueCrabSmall, 0, c
  BlueCrabSmall ->  BlackDrum, 0, d
   Temp ->  BlackDrum, 0, e
    BlueCrabSmall -> SpottedSeatrout, 0, f
   Temp ->  SpottedSeatrout, 0, g
    BlueCrabSmall -> HardheadCatfish, 0, h
   Temp ->  HardheadCatfish, 0, i
      BlueCrabSmall -> GafftopsailCatfish, 0, j
   Temp ->  GafftopsailCatfish, 0, k
  BrownShrimp ->  RedDrum, 0, l
   Temp ->  BrownShrimp, 0, m
  BrownShrimp ->  BlackDrum, 0, n
  BrownShrimp -> SpottedSeatrout, 0, o
  BrownShrimp -> HardheadCatfish, 0, p
  BrownShrimp -> GafftopsailCatfish, 0, q
   Salinity ->  RedDrum, 0, r
    Salinity ->  BlueCrabSmall, 0, s
    Salinity ->  BlackDrum, 0, t
    Salinity ->  SpottedSeatrout, 0, u
    Salinity ->  HardheadCatfish, 0, v
    Salinity ->  GafftopsailCatfish, 0, w
    Salinity ->  BrownShrimp, 0, x
    BlueCrabSmall ->  RedDrum, 1, aa
   Temp ->  RedDrum, 1, bb
   Temp ->  BlueCrabSmall, 1, cc
   BlueCrabSmall -> BlueCrabSmall, 1, dd
   RedDrum -> RedDrum, 1, ee
  BlueCrabSmall ->  BlackDrum, 1, ff
   Temp ->  BlackDrum, 1, gg
   BlackDrum -> BlackDrum, 1, hh
    BlueCrabSmall -> SpottedSeatrout, 1, ii
   Temp ->  SpottedSeatrout, 1, jj
   SpottedSeatrout -> SpottedSeatrout, 1, kk
    BlueCrabSmall -> HardheadCatfish, 1, ll
   Temp ->  HardheadCatfish, 1, mm
  HardheadCatfish -> HardheadCatfish, 1, nn
      BlueCrabSmall -> GafftopsailCatfish, 1, oo
   Temp ->  GafftopsailCatfish, 1, pp
 GafftopsailCatfish-> GafftopsailCatfish, 1, qq
  BrownShrimp ->  RedDrum, 1, rr
   Temp ->  BrownShrimp, 1, ss
   BrownShrimp -> BrownShrimp, 1, tt
  BrownShrimp ->  BlackDrum, 1, uu
  BrownShrimp -> SpottedSeatrout, 1, vv
  BrownShrimp -> HardheadCatfish, 1, ww
  BrownShrimp -> GafftopsailCatfish, 1, xx
   Salinity ->  RedDrum, 1, yy
    Salinity ->  BlueCrabSmall, 1, zz
    Salinity ->  BlackDrum, 1, aaa
    Salinity ->  SpottedSeatrout, 1, bbb
    Salinity ->  HardheadCatfish, 1, ccc
    Salinity ->  GafftopsailCatfish, 1, ddd
    Salinity ->  BrownShrimp, 1, eee
"

### note: these data are already log and z score transformed for annual means (as opposed to raw values)

# get variables of interest in a time series 
data_GB <- ts(GalvestonBay_GN_YearMean_log[,(c("BlueCrabSmall","Temp", "Salinity", "RedDrum", "SpottedSeatrout", "BlackDrum", "HardheadCatfish", "GafftopsailCatfish", "BrownShrimp"))], 
              start = c(min(GalvestonBay_GN_YearMean_log$YEAR))) 

data_AB <- ts(AransasBay_GN_YearMean_log[,(c("BlueCrabSmall","Temp", "Salinity", "RedDrum", "SpottedSeatrout", "BlackDrum", "HardheadCatfish", "GafftopsailCatfish", "BrownShrimp"))], 
              start = c(min(GalvestonBay_GN_YearMean_log$YEAR))) 


# fit the dsem (model) for GB
fitGB = dsem( sem = sem1,
             tsdata = data_GB,
             family = rep("fixed", ncol(data_GB)),             
             estimate_delta0 = TRUE,
             control = dsem_control(
               quiet = TRUE))
summary(fitGB)

# extract the coefficient matrix for GB 
coef_matrix_GBlag1 <- as_fitted_DAG(fitGB, lag = 1)$coef
coef_matrix_GBlag0 <- as_fitted_DAG(fitGB, lag = 0)$coef

# fit the dsem (model) for AB
fitAB = dsem( sem = sem1,
              tsdata = data_AB,
              family = rep("fixed", ncol(data_AB)),             
              estimate_delta0 = TRUE,
              control = dsem_control(
                quiet = TRUE))
summary(fitAB)

# extract the coefficient matrix for AN
coef_matrix_ABlag1 <- as_fitted_DAG(fitAB, lag = 1)$coef
coef_matrix_ABlag0 <- as_fitted_DAG(fitAB, lag = 0)$coef


# convert the matrix to a data frame for easier manipulation
coef_df <- as.data.frame(as.table(coef_matrix))


### FUNCTION to make path diagram and show in Plots panel

plot_qgraph_panel <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram") {
  
  # extract the column names (variable names) from the time series object
  abbrev_names <- colnames(data_ts)
  
  # check if the time series object has column names
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  
  # set plot margins (bottom, left, top, right) to add cushion space
  par(mar = c(25, 25, 25, 25))  
  
  # plot the graph
  qgraph(coef_matrix,         
         layout = "spring",              
         edge.labels = TRUE,             
         posCol = "navy",                
         negCol = "red3",                 
         labels = abbrev_names,          
         title = plot_title,             
         title.cex = 1.5,                 
         label.cex = .8,                  
         vsize = 10,
         label.scale = FALSE,
         shape = "ellipse")
}


### FUNCTION to make path diagram and export to desktop
export_qgraph <- function(data_ts, coef_matrix, plot_title = "Fitted DSEM Path Diagram", file_name) {
  
  # extract the column names (variable names) from the time series object
  abbrev_names <- colnames(data_ts)
  
  # check if the time series object has column names
  if (is.null(abbrev_names)) {
    stop("The time series object does not have column names.")
  }
  
  png(filename = file_name, width = 1200, height = 800)
  on.exit(dev.off(), add = TRUE)  
  
  # set plot margins (bottom, left, top, right) to add cushion space
  par(mar = c(25, 25, 25, 25))
  
  # plot the graph to the file
  qgraph(coef_matrix,         
         layout = "spring",              
         edge.labels = TRUE,            
         posCol = "navy",                
         negCol = "red3",                 
         labels = abbrev_names,          
         title = plot_title,             
         title.cex = 2.5,                 
         label.cex = 1.5,                
         vsize = 10,
         label.scale = FALSE,
         shape = "ellipse")
}

# make and export path diagrams for GB (one for lag, one for no lag)
plot_qgraph_panel(data_GB, coef_matrix_GBlag0, plot_title = "Galveston Bay - Sciaenid Food Web - No Lag")
export_qgraph(data_GB, coef_matrix_GBlag0, plot_title = "Galveston Bay - Sciaenid Food Web - No Lag", 
              file_name = "Galveston Bay - Sciaenid Food Web - No Lag.png")

plot_qgraph_panel(data_GB, coef_matrix_GBlag1, plot_title = "Galveston Bay - Sciaenid Food Web - 1 Year Lag")
export_qgraph(data_GB, coef_matrix_GBlag1, plot_title = "Galveston Bay - Sciaenid Food Web - 1 Year Lag", 
              file_name = "Galveston Bay - Sciaenid Food Web - 1 Year Lag.png")

# make and export path diagrams for AB (one for lag, one for no lag)
plot_qgraph_panel(data_AB, coef_matrix_ABlag0, plot_title = "Aransas Bay - Sciaenid Food Web - No Lag")
export_qgraph(data_AB, coef_matrix_ABlag0, plot_title = "Aransas Bay - Sciaenid Food Web - No Lag", 
              file_name = "Aransas Bay - Sciaenid Food Web - No Lag.png")

plot_qgraph_panel(data_AB, coef_matrix_ABlag1, plot_title = "Aransas Bay - Sciaenid Food Web - 1 Year Lag")
export_qgraph(data_AB, coef_matrix_ABlag1, plot_title = "Aransas Bay - Sciaenid Food Web - 1 Year Lag", 
              file_name = "Aransas Bay - Sciaenid Food Web - 1 Year Lag.png")


