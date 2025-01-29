
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

devtools::install_github("James-Thorson-NOAA/dsem", ref="main")


###################################################################

# fitting 4 DSEMs (two trophic systems and two major bays) at annual time steps 

# all analyses/code come from 4 excel files:
# AransasBay_Sciaenid_Wide_Trans 
# AransasBay_Pred_Wide_Trans
# GalvestonBay_Sciaenid_Wide_Trans 
# GalvestonBay_Pred_Wide_Trans


# first system and Major Bay: Sciaenid Trophic System and Aransas Bay 

# get data into a time series object from all numeric columns except YEAR
AB_Sciaenid_TS <- ts(
  AransasBay_Sciaenid_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Sciaenid_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Sciaenid_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the sciaenid system for aransas bay for zero and 1 year lags
# the 3 seperate chunks are just to keep track of the 3 minor bays... ignore my inconsistent letter identifiers
semAB_Sciaenid <- "
    BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, a
    Temperature_AransasBay -> RedDrum_AransasBay, 0, b
    Temperature_AransasBay -> BlueCrabSmall_AransasBay, 0, c
    BlueCrabSmall_AransasBay -> SouthernFlounder_AransasBay, 0, d
    Temperature_AransasBay -> SouthernFlounder_AransasBay, 0, e
    BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, f
    Temperature_AransasBay -> SpottedSeatrout_AransasBay, 0, g
    BrownShrimp_AransasBay -> RedDrum_AransasBay, 0, l
    Temperature_AransasBay -> BrownShrimp_AransasBay, 0, m
    BrownShrimp_AransasBay -> SouthernFlounder_AransasBay, 0, n
    BrownShrimp_AransasBay -> SpottedSeatrout_AransasBay, 0, o
    Salinity_AransasBay -> RedDrum_AransasBay, 0, r
    Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s
    Salinity_AransasBay -> SouthernFlounder_AransasBay, 0, t
    Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, u
    Salinity_AransasBay -> BrownShrimp_AransasBay, 0, x
    BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 1, aa
    Temperature_AransasBay -> RedDrum_AransasBay, 1, bb
    Temperature_AransasBay -> BlueCrabSmall_AransasBay, 1, cc
    BlueCrabSmall_AransasBay -> SouthernFlounder_AransasBay, 1, ff
    Temperature_AransasBay -> SouthernFlounder_AransasBay, 1, gg
    BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 1, ii
    Temperature_AransasBay -> SpottedSeatrout_AransasBay, 1, jj
    BrownShrimp_AransasBay -> RedDrum_AransasBay, 1, rr
    Temperature_AransasBay -> BrownShrimp_AransasBay, 1, ss
    BrownShrimp_AransasBay -> SouthernFlounder_AransasBay, 1, uu
    BrownShrimp_AransasBay -> SpottedSeatrout_AransasBay, 1, vv
    Salinity_AransasBay -> RedDrum_AransasBay, 1, yy
    Salinity_AransasBay -> BlueCrabSmall_AransasBay, 1, zz
    Salinity_AransasBay -> SouthernFlounder_AransasBay, 1, aaa
    Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, bbb
    Salinity_AransasBay -> BrownShrimp_AransasBay, 1, eee

    BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, abcd
    Temperature_CopanoBay -> RedDrum_CopanoBay, 0, efgh
    Temperature_CopanoBay -> BlueCrabSmall_CopanoBay, 0, ijkl
    BlueCrabSmall_CopanoBay -> SouthernFlounder_CopanoBay, 0, mnop
    Temperature_CopanoBay -> SouthernFlounder_CopanoBay, 0, qrst
    BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, uvwx
    Temperature_CopanoBay -> SpottedSeatrout_CopanoBay, 0, yzab
    BrownShrimp_CopanoBay -> RedDrum_CopanoBay, 0, stuv
    Temperature_CopanoBay -> BrownShrimp_CopanoBay, 0, wxyz
    BrownShrimp_CopanoBay -> SouthernFlounder_CopanoBay, 0, abcd
    BrownShrimp_CopanoBay -> SpottedSeatrout_CopanoBay, 0, efgh
    Salinity_CopanoBay -> RedDrum_CopanoBay, 0, qrst
    Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, uvwx
    Salinity_CopanoBay -> SouthernFlounder_CopanoBay, 0, yzab
    Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, cdef
    Salinity_CopanoBay -> BrownShrimp_CopanoBay, 0, opqr
    BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 1, stuv
    Temperature_CopanoBay -> RedDrum_CopanoBay, 1, wxyz
    Temperature_CopanoBay -> BlueCrabSmall_CopanoBay, 1, abcd
    BlueCrabSmall_CopanoBay -> SouthernFlounder_CopanoBay, 1, mnop
    Temperature_CopanoBay -> SouthernFlounder_CopanoBay, 1, qrst
    BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 1, yzab
    Temperature_CopanoBay -> SpottedSeatrout_CopanoBay, 1, cdef
    BrownShrimp_CopanoBay -> RedDrum_CopanoBay, 1, ijkl
    Temperature_CopanoBay -> BrownShrimp_CopanoBay, 1, mnop
    BrownShrimp_CopanoBay -> SouthernFlounder_CopanoBay, 1, uvwx
    BrownShrimp_CopanoBay -> SpottedSeatrout_CopanoBay, 1, yzab
    Salinity_CopanoBay -> RedDrum_CopanoBay, 1, klmn
    Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 1, opqr
    Salinity_CopanoBay -> SouthernFlounder_CopanoBay, 1, stuv
    Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, wxyz
    Salinity_CopanoBay -> BrownShrimp_CopanoBay, 1, ijkl

    BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, abcde
    Temperature_MesquiteBay -> RedDrum_MesquiteBay, 0, fghij
    Temperature_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, klmno
    BlueCrabSmall_MesquiteBay -> SouthernFlounder_MesquiteBay, 0, pqrst
    Temperature_MesquiteBay -> SouthernFlounder_MesquiteBay, 0, uvwxy
    BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, zabcd
    Temperature_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, efghi
    BrownShrimp_MesquiteBay -> RedDrum_MesquiteBay, 0, ghijk
    Temperature_MesquiteBay -> BrownShrimp_MesquiteBay, 0, lmnop
    BrownShrimp_MesquiteBay -> SouthernFlounder_MesquiteBay, 0, qrstu
    BrownShrimp_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, vwxyz
    Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, klmno
    Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, pqrst
    Salinity_MesquiteBay -> SouthernFlounder_MesquiteBay, 0, uvwxy
    Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, zabcd
    Salinity_MesquiteBay -> BrownShrimp_MesquiteBay, 0, opqrs
    BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 1, tuvwz
    Temperature_MesquiteBay -> RedDrum_MesquiteBay, 1, abcdf
    Temperature_MesquiteBay -> BlueCrabSmall_MesquiteBay, 1, ghijk
    BlueCrabSmall_MesquiteBay -> SouthernFlounder_MesquiteBay, 1, vwxyz
    Temperature_MesquiteBay -> SouthernFlounder_MesquiteBay, 1, abcde
    BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, klmno
    Temperature_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, pqrst
    BrownShrimp_MesquiteBay -> RedDrum_MesquiteBay, 1, ghijk
    Temperature_MesquiteBay -> BrownShrimp_MesquiteBay, 1, lmnop
    BrownShrimp_MesquiteBay -> BrownShrimp_MesquiteBay, 1, qrstu
    BrownShrimp_MesquiteBay -> SouthernFlounder_MesquiteBay, 1, vwxyz
    BrownShrimp_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, abcde
    Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, pqrst
    Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 1, uvwxy
    Salinity_MesquiteBay -> SouthernFlounder_MesquiteBay, 1, zabcd
    Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, efghi
    Salinity_MesquiteBay -> BrownShrimp_MesquiteBay, 1, tuvwz

"

# fit the dsem model for the sciaenid system for aransas bay 
fitAB_Sciaenid = dsem( sem = semAB_Sciaenid,
              tsdata = AB_Sciaenid_TS,
              control = dsem_control(
                quiet = TRUE))
summary(fitAB_Sciaenid)

# i am not sure what this code is doing
# i think it takes the time series for each site (in our case, three minor bays) and combines their output/fit
# see fig 3 of thorson et al. 2024 "Note that the path coefficient estimates shown here are specified as identical across
# all 12 sites that are treated as separate time series for each species (48 time series total)" 
get_part_Sciaenid = function(x){
  vars = c("Salinity","BrownShrimp","SpottedSeatrout","Temperature","SouthernFlounder","RedDrum", "BlueCrabSmall")
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
  par(mar = c(25, 25, 25, 25))  # Adjust the values to add more space if needed
  
  # Plot the graph
  qgraph(coef_matrix,         
         layout = "spring",              # Layout style
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


# extract the coefficient matrix 
coef_matrix_AB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fitAB_Sciaenid, lag=0))
coef_matrix_AB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fitAB_Sciaenid, lag=1))

# create an empty data frame with specified column names
df_Sciaenid <- data.frame(Salinity = numeric(),
                 BrownShrimp = numeric(),
                 SpottedSeatrout = numeric(),
                 Temperature = numeric(),
                 SouthernFlounder = numeric(),
                 RedDrum = numeric(),
                 BlueCrabSmall = numeric(),
                 stringsAsFactors = FALSE)

# make path diagrams
png("sciaenidAB-path-annual.png", width = 2000, height = 1300, res = 300)
par(mfrow = c(1, 2))
plot1 <- plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_NoLag$coef,  plot_title = "Sciaenid-AB-No Lag (Annual)")
plot2 <- plot_qgraph_panel(df_Sciaenid, coef_matrix_AB_Sciaenid_YesLag$coef, plot_title = "Sciaenid-AB-Yes Lag (Annual)")
par(mfrow = c(1, 1))
dev.off()




###################################################################




# Repeating the above but for the...
# Second System and Major Bay: Keystone Predator Trophic System and Aransas Bay 

# get data into a time series object from all numeric columns except YEAR
AB_Pred_TS <- ts(
  AransasBay_Pred_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(AransasBay_Pred_Wide_Trans$YEAR)),  
  end = c(max(AransasBay_Pred_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the keystone predator system for aransas bay for zero and 1 year lags
# the 3 separate chunks are just to keep track of the 3 minor bays... ignore my inconsistent letter identifiers
semAB_Pred <- "
    Temperature_AransasBay -> StripedMullet_AransasBay, 0, a
    Temperature_AransasBay -> BullShark_AransasBay, 0, b
    Temperature_AransasBay -> AlligatorGar_AransasBay, 0, c
    Salinity_AransasBay -> StripedMullet_AransasBay, 0, d
    Salinity_AransasBay -> BullShark_AransasBay, 0, e
    Salinity_AransasBay -> AlligatorGar_AransasBay, 0, f
    StripedMullet_AransasBay -> BullShark_AransasBay, 0, g
    StripedMullet_AransasBay -> AlligatorGar_AransasBay, 0, h
    Temperature_MesquiteBay -> StripedMullet_MesquiteBay , 1, i
    Temperature_MesquiteBay  -> BullShark_MesquiteBay , 1, j
    Temperature_MesquiteBay  -> AlligatorGar_MesquiteBay , 1, k
    Salinity_MesquiteBay  -> StripedMullet_MesquiteBay , 1, l
    Salinity_MesquiteBay  -> BullShark_MesquiteBay , 1, m
    Salinity_MesquiteBay  -> AlligatorGar_MesquiteBay , 1, n
    StripedMullet_MesquiteBay  -> BullShark_MesquiteBay , 1, o
    StripedMullet_MesquiteBay  -> AlligatorGar_MesquiteBay , 1, p
    Temperature_CopanoBay -> StripedMullet_CopanoBay, 0, q
    Temperature_CopanoBay -> BullShark_CopanoBay, 0, r
    Temperature_CopanoBay -> AlligatorGar_CopanoBay, 0, s
    Salinity_CopanoBay -> StripedMullet_CopanoBay, 0, t
    Salinity_CopanoBay -> BullShark_CopanoBay, 0, u
    Salinity_CopanoBay -> AlligatorGar_CopanoBay, 0, v
    StripedMullet_CopanoBay -> BullShark_CopanoBay, 0, w
    StripedMullet_CopanoBay -> AlligatorGar_CopanoBay, 0, x
    GulfMenhaden_AransasBay -> AlligatorGar_AransasBay, 0, aw
    GulfMenhaden_AransasBay -> BullShark_AransasBay, 0, ay
    GulfMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 0, ax
    GulfMenhaden_CopanoBay -> BullShark_CopanoBay, 0, az
    GulfMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 0, ba
    GulfMenhaden_MesquiteBay -> BullShark_MesquiteBay, 0, bb
    Temperature_CopanoBay -> GulfMenhaden_CopanoBay, 0, bc
    Salinity_CopanoBay -> GulfMenhaden_CopanoBay, 0, bd
    Temperature_AransasBay -> GulfMenhaden_AransasBay, 0, be
    Salinity_AransasBay -> GulfMenhaden_AransasBay, 0, bf
    Temperature_MesquiteBay -> GulfMenhaden_MesquiteBay, 0, bg
    Salinity_MesquiteBay -> GulfMenhaden_MesquiteBay, 0, bh

    Temperature_CopanoBay -> StripedMullet_CopanoBay, 1, y
    Temperature_CopanoBay -> BullShark_CopanoBay, 1, z
    Temperature_CopanoBay -> AlligatorGar_CopanoBay, 1, aa
    Salinity_CopanoBay -> StripedMullet_CopanoBay, 1, ab
    Salinity_CopanoBay -> BullShark_CopanoBay, 1, ac
    Salinity_CopanoBay -> AlligatorGar_CopanoBay, 1, ad
    StripedMullet_CopanoBay -> BullShark_CopanoBay, 1, ae
    StripedMullet_CopanoBay -> AlligatorGar_CopanoBay, 1, af
    Temperature_MesquiteBay -> StripedMullet_MesquiteBay, 0, ag
    Temperature_MesquiteBay -> BullShark_MesquiteBay, 0, ah
    Temperature_MesquiteBay -> AlligatorGar_MesquiteBay, 0, ai
    Salinity_MesquiteBay -> StripedMullet_MesquiteBay, 0, aj
    Salinity_MesquiteBay -> BullShark_MesquiteBay, 0, ak
    Salinity_MesquiteBay -> AlligatorGar_MesquiteBay, 0, al
    StripedMullet_MesquiteBay -> BullShark_MesquiteBay, 0, am
    StripedMullet_MesquiteBay -> AlligatorGar_MesquiteBay, 0, an
    Temperature_AransasBay -> StripedMullet_AransasBay, 1, ao
    Temperature_AransasBay -> BullShark_AransasBay, 1, ap
    Temperature_AransasBay-> AlligatorGar_AransasBay, 1, aq
    Salinity_AransasBay -> StripedMullet_AransasBay, 1, ar
    Salinity_AransasBay -> BullShark_AransasBay, 1, as
    Salinity_AransasBay -> AlligatorGar_AransasBay, 1, at
    StripedMullet_AransasBay -> BullShark_AransasBay, 1, au
    StripedMullet_AransasBay -> AlligatorGar_AransasBay, 1, av
    GulfMenhaden_AransasBay -> AlligatorGar_AransasBay, 1, bi
    GulfMenhaden_AransasBay -> BullShark_AransasBay, 1, bj
    GulfMenhaden_CopanoBay -> AlligatorGar_CopanoBay, 1, bk
    GulfMenhaden_CopanoBay -> BullShark_CopanoBay, 1, bl
    GulfMenhaden_MesquiteBay -> AlligatorGar_MesquiteBay, 1, bm
    GulfMenhaden_MesquiteBay -> BullShark_MesquiteBay, 1, bn
    Temperature_CopanoBay -> GulfMenhaden_CopanoBay, 1, bo
    Salinity_CopanoBay -> GulfMenhaden_CopanoBay, 1, bp
    Temperature_AransasBay -> GulfMenhaden_AransasBay, 1, bq
    Salinity_AransasBay -> GulfMenhaden_AransasBay, 1, br
    Temperature_MesquiteBay -> GulfMenhaden_MesquiteBay, 1, bs
    Salinity_MesquiteBay -> GulfMenhaden_MesquiteBay, 1, bt

"

# fit the dsem model for the keystone predator system for aransas bay 
fit_AB_Pred_TS = dsem(sem = semAB_Pred,
                      tsdata = AB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_AB_Pred_TS)

# i am not totally sure what this code is doing
# i think it takes the time series for each site (in our case, three minor bays) and combines their output/fit
# see fig 3 of thorson et al. 2024 "Note that the path coefficient estimates shown here are specified as identical across
# all 12 sites that are treated as separate time series for each species (48 time series total)" 
get_part_pred = function(x){
  vars = c("Salinity","StripedMullet","BullShark","AlligatorGar","Temperature", "GulfMenhaden")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# extract the coefficient matrix 
coef_matrix_AB_Pred_NoLag <- get_part_pred(as_fitted_DAG(fit_AB_Pred_TS, lag=0))
coef_matrix_AB_Pred_YesLag <- get_part_pred(as_fitted_DAG(fit_AB_Pred_TS, lag=1))

# create an empty data frame with specified column names
df_Pred <- data.frame(Salinity = numeric(),
                      StripedMullet = numeric(),
                      BullShark = numeric(),
                      AlligatorGar = numeric(),
                      Temperature = numeric(),
                      GulfMenhaden = numeric(),
                      stringsAsFactors = FALSE)

png("predAB-path-annual.png", width = 2000, height = 1300, res = 300)
par(mfrow = c(1, 2))
plot1 <- plot_qgraph_panel(df_Pred, coef_matrix_AB_Pred_NoLag$coef,  plot_title = "Keystone Predator-AB-No Lag (Annual)")
plot2 <- plot_qgraph_panel(df_Pred, coef_matrix_AB_Pred_YesLag$coef, plot_title = "Keystone Predator-AB-Yes Lag (Annual)")
par(mfrow = c(1, 1))
dev.off()




###################################################################




# Repeating the above but for the...
# Third System and Major Bay: Sciaenid Predator Trophic System and Galveston Bay 

# get data into a time series object from all numeric columns except YEAR
GB_Sciaenid_TS <- ts(
  GalvestonBay_Sciaenid_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(GalvestonBay_Sciaenid_Wide_Trans$YEAR)),  
  end = c(max(GalvestonBay_Sciaenid_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

semGB_Sciaenid <- "     BlueCrabSmall_GalvestonBay -> RedDrum_GalvestonBay, 0, a
    Temperature_GalvestonBay -> RedDrum_GalvestonBay, 0, b
    Temperature_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, c
    BlueCrabSmall_GalvestonBay -> SouthernFlounder_GalvestonBay, 0, d
    Temperature_GalvestonBay -> SouthernFlounder_GalvestonBay, 0, e
    BlueCrabSmall_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, f
    Temperature_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, g
    BrownShrimp_GalvestonBay -> RedDrum_GalvestonBay, 0, l
    Temperature_GalvestonBay -> BrownShrimp_GalvestonBay, 0, m
    BrownShrimp_GalvestonBay -> SouthernFlounder_GalvestonBay, 0, n
    BrownShrimp_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, o
    Salinity_GalvestonBay -> RedDrum_GalvestonBay, 0, r
    Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 0, s
    Salinity_GalvestonBay -> SouthernFlounder_GalvestonBay, 0, t
    Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 0, u
    Salinity_GalvestonBay -> BrownShrimp_GalvestonBay, 0, x
    BlueCrabSmall_GalvestonBay -> RedDrum_GalvestonBay, 1, aa
    Temperature_GalvestonBay -> RedDrum_GalvestonBay, 1, bb
    Temperature_GalvestonBay -> BlueCrabSmall_GalvestonBay, 1, cc
    BlueCrabSmall_GalvestonBay -> SouthernFlounder_GalvestonBay, 1, ff
    Temperature_GalvestonBay -> SouthernFlounder_GalvestonBay, 1, gg
    BlueCrabSmall_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, ii
    Temperature_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, jj
    BrownShrimp_GalvestonBay -> RedDrum_GalvestonBay, 1, rr
    Temperature_GalvestonBay -> BrownShrimp_GalvestonBay, 1, ss
    BrownShrimp_GalvestonBay -> SouthernFlounder_GalvestonBay, 1, uu
    BrownShrimp_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, vv
    Salinity_GalvestonBay -> RedDrum_GalvestonBay, 1, yy
    Salinity_GalvestonBay -> BlueCrabSmall_GalvestonBay, 1, zz
    Salinity_GalvestonBay -> SouthernFlounder_GalvestonBay, 1, aaa
    Salinity_GalvestonBay -> SpottedSeatrout_GalvestonBay, 1, bbb
    Salinity_GalvestonBay -> BrownShrimp_GalvestonBay, 1, eee
    BlueCrabSmall_TrinityBay -> RedDrum_TrinityBay, 0, abcd
    Temperature_TrinityBay -> RedDrum_TrinityBay, 0, efgh
    Temperature_TrinityBay -> BlueCrabSmall_TrinityBay, 0, ijkl
    BlueCrabSmall_TrinityBay -> SouthernFlounder_TrinityBay, 0, mnop
    Temperature_TrinityBay -> SouthernFlounder_TrinityBay, 0, qrst
    BlueCrabSmall_TrinityBay -> SpottedSeatrout_TrinityBay, 0, uvwx
    Temperature_TrinityBay -> SpottedSeatrout_TrinityBay, 0, yzab
    BrownShrimp_TrinityBay -> RedDrum_TrinityBay, 0, stuv
    Temperature_TrinityBay -> BrownShrimp_TrinityBay, 0, wxyz
    BrownShrimp_TrinityBay -> SouthernFlounder_TrinityBay, 0, abcd
    BrownShrimp_TrinityBay -> SpottedSeatrout_TrinityBay, 0, efgh
    Salinity_TrinityBay -> RedDrum_TrinityBay, 0, qrst
    Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 0, uvwx
    Salinity_TrinityBay -> SouthernFlounder_TrinityBay, 0, yzab
    Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 0, cdef
    Salinity_TrinityBay -> BrownShrimp_TrinityBay, 0, opqr
    BlueCrabSmall_TrinityBay -> RedDrum_TrinityBay, 1, stuv
    Temperature_TrinityBay -> RedDrum_TrinityBay, 1, wxyz
    Temperature_TrinityBay -> BlueCrabSmall_TrinityBay, 1, abcd
    BlueCrabSmall_TrinityBay -> SouthernFlounder_TrinityBay, 1, mnop
    Temperature_TrinityBay -> SouthernFlounder_TrinityBay, 1, qrst
    BlueCrabSmall_TrinityBay -> SpottedSeatrout_TrinityBay, 1, yzab
    Temperature_TrinityBay -> SpottedSeatrout_TrinityBay, 1, cdef
    BrownShrimp_TrinityBay -> RedDrum_TrinityBay, 1, ijkl
    Temperature_TrinityBay -> BrownShrimp_TrinityBay, 1, mnop
    BrownShrimp_TrinityBay -> SouthernFlounder_TrinityBay, 1, uvwx
    BrownShrimp_TrinityBay -> SpottedSeatrout_TrinityBay, 1, yzab
    Salinity_TrinityBay -> RedDrum_TrinityBay, 1, klmn
    Salinity_TrinityBay -> BlueCrabSmall_TrinityBay, 1, opqr
    Salinity_TrinityBay -> SouthernFlounder_TrinityBay, 1, stuv
    Salinity_TrinityBay -> SpottedSeatrout_TrinityBay, 1, wxyz
    Salinity_TrinityBay -> BrownShrimp_TrinityBay, 1, ijkl
    BlueCrabSmall_WestBay -> RedDrum_WestBay, 0, abcde
    Temperature_WestBay -> RedDrum_WestBay, 0, fghij
    Temperature_WestBay -> BlueCrabSmall_WestBay, 0, klmno
    BlueCrabSmall_WestBay -> SouthernFlounder_WestBay, 0, pqrst
    Temperature_WestBay -> SouthernFlounder_WestBay, 0, uvwxy
    BlueCrabSmall_WestBay -> SpottedSeatrout_WestBay, 0, zabcd
    Temperature_WestBay -> SpottedSeatrout_WestBay, 0, efghi
    BrownShrimp_WestBay -> RedDrum_WestBay, 0, ghijk
    Temperature_WestBay -> BrownShrimp_WestBay, 0, lmnop
    BrownShrimp_WestBay -> SouthernFlounder_WestBay, 0, qrstu
    BrownShrimp_WestBay -> SpottedSeatrout_WestBay, 0, vwxyz
    Salinity_WestBay -> RedDrum_WestBay, 0, klmno
    Salinity_WestBay -> BlueCrabSmall_WestBay, 0, pqrst
    Salinity_WestBay -> SouthernFlounder_WestBay, 0, uvwxy
    Salinity_WestBay -> SpottedSeatrout_WestBay, 0, zabcd
    Salinity_WestBay -> BrownShrimp_WestBay, 0, opqrs
    BlueCrabSmall_WestBay -> RedDrum_WestBay, 1, tuvwz
    Temperature_WestBay -> RedDrum_WestBay, 1, abcdf
    Temperature_WestBay -> BlueCrabSmall_WestBay, 1, ghijk
    BlueCrabSmall_WestBay -> SouthernFlounder_WestBay, 1, vwxyz
    Temperature_WestBay -> SouthernFlounder_WestBay, 1, abcde
    BlueCrabSmall_WestBay -> SpottedSeatrout_WestBay, 1, klmno
    Temperature_WestBay -> SpottedSeatrout_WestBay, 1, pqrst
    BrownShrimp_WestBay -> RedDrum_WestBay, 1, ghijk
    Temperature_WestBay -> BrownShrimp_WestBay, 1, lmnop
    BrownShrimp_WestBay -> BrownShrimp_WestBay, 1, qrstu
    BrownShrimp_WestBay -> SouthernFlounder_WestBay, 1, vwxyz
    BrownShrimp_WestBay -> SpottedSeatrout_WestBay, 1, abcde
    Salinity_WestBay -> RedDrum_WestBay, 1, pqrst
    Salinity_WestBay -> BlueCrabSmall_WestBay, 1, uvwxy
    Salinity_WestBay -> SouthernFlounder_WestBay, 1, zabcd
    Salinity_WestBay -> SpottedSeatrout_WestBay, 1, efghi
    Salinity_WestBay -> BrownShrimp_WestBay, 1, lmnop
"

# fit the dsem model for the sciaenid system for aransas bay 
fitGB_Sciaenid = dsem( sem = semGB_Sciaenid,
                       tsdata = GB_Sciaenid_TS,
                       control = dsem_control(
                         quiet = TRUE))
summary(fitGB_Sciaenid)

# extract the coefficient matrix 
coef_matrix_GB_Sciaenid_NoLag <- get_part_Sciaenid(as_fitted_DAG(fitGB_Sciaenid, lag=0))
coef_matrix_GB_Sciaenid_YesLag <- get_part_Sciaenid(as_fitted_DAG(fitGB_Sciaenid, lag=1))

# create an empty data frame with specified column names
df_Sciaenid_GB <- data.frame(Salinity = numeric(),
                             BrownShrimp = numeric(),
                             SpottedSeatrout = numeric(),
                             Temperature = numeric(),
                             SouthernFlounder = numeric(),
                             RedDrum = numeric(),
                             BlueCrabSmall = numeric(),
                             stringsAsFactors = FALSE)

# make path diagrams
png("sciaenidGB-path-annual.png", width = 2000, height = 1300, res = 300)
par(mfrow = c(1, 2))
plot1 <- plot_qgraph_panel(df_Sciaenid_GB, coef_matrix_GB_Sciaenid_NoLag$coef,  plot_title = "Sciaenid-GB-No Lag (Annual)")
plot2 <- plot_qgraph_panel(df_Sciaenid_GB, coef_matrix_GB_Sciaenid_YesLag$coef, plot_title = "Sciaenid-GB-Yes Lag (Annual)")
par(mfrow = c(1, 1))
dev.off()




###################################################################



# Repeating the above but for the...
# Fourth System and Major Bay: Keystone Predator Trophic System and Galveston Bay 

# get data into a time series object from all numeric columns except YEAR
GB_Pred_TS <- ts(
  GalvestonBay_Pred_Wide_Trans %>%
    select(-YEAR),  
  start = c(min(GalvestonBay_Pred_Wide_Trans$YEAR)),  
  end = c(max(GalvestonBay_Pred_Wide_Trans$YEAR)),
  frequency = 1  # assuming yearly data
)

# specify the dsem model for the keystone predator system for aransas bay for zero and 1 year lags
# the 3 separate chunks are just to keep track of the 3 minor bays... ignore my inconsistent letter identifiers
semGB_Pred <- "
Temperature_GalvestonBay -> StripedMullet_GalvestonBay, 0, a
Temperature_GalvestonBay -> BullShark_GalvestonBay, 0, b
Temperature_GalvestonBay -> AlligatorGar_GalvestonBay, 0, c
Salinity_GalvestonBay -> StripedMullet_GalvestonBay, 0, d
Salinity_GalvestonBay -> BullShark_GalvestonBay, 0, e
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 0, f
StripedMullet_GalvestonBay -> BullShark_GalvestonBay, 0, g
StripedMullet_GalvestonBay -> AlligatorGar_GalvestonBay, 0, h
Temperature_WestBay -> StripedMullet_WestBay, 1, i
Temperature_WestBay -> BullShark_WestBay, 1, j
Temperature_WestBay -> AlligatorGar_WestBay, 1, k
Salinity_WestBay -> StripedMullet_WestBay, 1, l
Salinity_WestBay -> BullShark_WestBay, 1, m
Salinity_WestBay -> AlligatorGar_WestBay, 1, n
StripedMullet_WestBay -> BullShark_WestBay, 1, o
StripedMullet_WestBay -> AlligatorGar_WestBay, 1, p
Temperature_TrinityBay -> StripedMullet_TrinityBay, 0, q
Temperature_TrinityBay -> BullShark_TrinityBay, 0, r
Temperature_TrinityBay -> AlligatorGar_TrinityBay, 0, s
Salinity_TrinityBay -> StripedMullet_TrinityBay, 0, t
Salinity_TrinityBay -> BullShark_TrinityBay, 0, u
Salinity_TrinityBay -> AlligatorGar_TrinityBay, 0, v
StripedMullet_TrinityBay -> BullShark_TrinityBay, 0, w
StripedMullet_TrinityBay -> AlligatorGar_TrinityBay, 0, x
GulfMenhaden_GalvestonBay -> AlligatorGar_GalvestonBay, 0, aw
GulfMenhaden_GalvestonBay -> BullShark_GalvestonBay, 0, ay
GulfMenhaden_TrinityBay -> AlligatorGar_TrinityBay, 0, ax
GulfMenhaden_TrinityBay -> BullShark_TrinityBay, 0, az
GulfMenhaden_WestBay -> AlligatorGar_WestBay, 0, ba
GulfMenhaden_WestBay -> BullShark_WestBay, 0, bb
Temperature_TrinityBay -> GulfMenhaden_TrinityBay, 0, bc
Salinity_TrinityBay -> GulfMenhaden_TrinityBay, 0, bd
Temperature_GalvestonBay -> GulfMenhaden_GalvestonBay, 0, be
Salinity_GalvestonBay -> GulfMenhaden_GalvestonBay, 0, bf
Temperature_WestBay -> GulfMenhaden_WestBay, 0, bg
Salinity_WestBay -> GulfMenhaden_WestBay, 0, bh
Temperature_TrinityBay -> StripedMullet_TrinityBay, 1, y
Temperature_TrinityBay -> BullShark_TrinityBay, 1, z
Temperature_TrinityBay -> AlligatorGar_TrinityBay, 1, aa
Salinity_TrinityBay -> StripedMullet_TrinityBay, 1, ab
Salinity_TrinityBay -> BullShark_TrinityBay, 1, ac
Salinity_TrinityBay -> AlligatorGar_TrinityBay, 1, ad
StripedMullet_TrinityBay -> BullShark_TrinityBay, 1, ae
StripedMullet_TrinityBay -> AlligatorGar_TrinityBay, 1, af
Temperature_WestBay -> StripedMullet_WestBay, 0, ag
Temperature_WestBay -> BullShark_WestBay, 0, ah
Temperature_WestBay -> AlligatorGar_WestBay, 0, ai
Salinity_WestBay -> StripedMullet_WestBay, 0, aj
Salinity_WestBay -> BullShark_WestBay, 0, ak
Salinity_WestBay -> AlligatorGar_WestBay, 0, al
StripedMullet_WestBay -> BullShark_WestBay, 0, am
StripedMullet_WestBay -> AlligatorGar_WestBay, 0, an
Temperature_GalvestonBay -> StripedMullet_GalvestonBay, 1, ao
Temperature_GalvestonBay -> BullShark_GalvestonBay, 1, ap
Temperature_GalvestonBay -> AlligatorGar_GalvestonBay, 1, aq
Salinity_GalvestonBay -> StripedMullet_GalvestonBay, 1, ar
Salinity_GalvestonBay -> BullShark_GalvestonBay, 1, as
Salinity_GalvestonBay -> AlligatorGar_GalvestonBay, 1, at
StripedMullet_GalvestonBay -> BullShark_GalvestonBay, 1, au
StripedMullet_GalvestonBay -> AlligatorGar_GalvestonBay, 1, av
GulfMenhaden_GalvestonBay -> AlligatorGar_GalvestonBay, 1, bi
GulfMenhaden_GalvestonBay -> BullShark_GalvestonBay, 1, bj
GulfMenhaden_TrinityBay -> AlligatorGar_TrinityBay, 1, bk
GulfMenhaden_TrinityBay -> BullShark_TrinityBay, 1, bl
GulfMenhaden_WestBay -> AlligatorGar_WestBay, 1, bm
GulfMenhaden_WestBay -> BullShark_WestBay, 1, bn
Temperature_TrinityBay -> GulfMenhaden_TrinityBay, 1, bo
Salinity_TrinityBay -> GulfMenhaden_TrinityBay, 1, bp
Temperature_GalvestonBay -> GulfMenhaden_GalvestonBay, 1, bq
Salinity_GalvestonBay -> GulfMenhaden_GalvestonBay, 1, br
Temperature_WestBay -> GulfMenhaden_WestBay, 1, bs
Salinity_WestBay -> GulfMenhaden_WestBay, 1, bt
"

# fit the dsem model for the keystone predator system for aransas bay 
fit_GB_Pred_TS = dsem(sem = semGB_Pred,
                      tsdata = GB_Pred_TS,
                      control = dsem_control(
                        quiet = TRUE))
summary(fit_GB_Pred_TS)

# i am not totally sure what this code is doing
# i think it takes the time series for each site (in our case, three minor bays) and combines their output/fit
# see fig 3 of thorson et al. 2024 "Note that the path coefficient estimates shown here are specified as identical across
# all 12 sites that are treated as separate time series for each species (48 time series total)" 
get_part_pred = function(x){
  vars = c("Salinity","StripedMullet","BullShark","AlligatorGar","Temperature", "GulfMenhaden")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# extract the coefficient matrix 
coef_matrix_GB_Pred_NoLag <- get_part_pred(as_fitted_DAG(fit_GB_Pred_TS, lag=0))
coef_matrix_GB_Pred_YesLag <- get_part_pred(as_fitted_DAG(fit_GB_Pred_TS, lag=1))

# create an empty data frame with specified column names
df_Pred_GB <- data.frame(Salinity = numeric(),
                      StripedMullet = numeric(),
                      BullShark = numeric(),
                      AlligatorGar = numeric(),
                      Temperature = numeric(),
                      GulfMenhaden = numeric(),
                      stringsAsFactors = FALSE)

png("predGB-path-annual.png", width = 2000, height = 1300, res = 300)
par(mfrow = c(1, 2))
plot1 <- plot_qgraph_panel(df_Pred_GB, coef_matrix_GB_Pred_NoLag$coef,  plot_title = "Keystone Predator-GB-No Lag (Annual)")
plot2 <- plot_qgraph_panel(df_Pred_GB, coef_matrix_GB_Pred_YesLag$coef, plot_title = "Keystone Predator-GB-Yes Lag (Annual)")
par(mfrow = c(1, 1))
dev.off()


