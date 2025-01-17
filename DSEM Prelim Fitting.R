
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
sem1 <- "
    BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 0, a
    Temperature_AransasBay -> RedDrum_AransasBay, 0, b
    Temperature_AransasBay -> BlueCrabSmall_AransasBay, 0, c
    BlueCrabSmall_AransasBay -> BlackDrum_AransasBay, 0, d
    Temperature_AransasBay -> BlackDrum_AransasBay, 0, e
    BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 0, f
    Temperature_AransasBay -> SpottedSeatrout_AransasBay, 0, g
    BlueCrabSmall_AransasBay -> HardheadCatfish_AransasBay, 0, h
    Temperature_AransasBay -> HardheadCatfish_AransasBay, 0, i
    BlueCrabSmall_AransasBay -> GafftopsailCatfish_AransasBay, 0, j
    Temperature_AransasBay -> GafftopsailCatfish_AransasBay, 0, k
    BrownShrimp_AransasBay -> RedDrum_AransasBay, 0, l
    Temperature_AransasBay -> BrownShrimp_AransasBay, 0, m
    BrownShrimp_AransasBay -> BlackDrum_AransasBay, 0, n
    BrownShrimp_AransasBay -> SpottedSeatrout_AransasBay, 0, o
    BrownShrimp_AransasBay -> HardheadCatfish_AransasBay, 0, p
    BrownShrimp_AransasBay -> GafftopsailCatfish_AransasBay, 0, q
    Salinity_AransasBay -> RedDrum_AransasBay, 0, r
    Salinity_AransasBay -> BlueCrabSmall_AransasBay, 0, s
    Salinity_AransasBay -> BlackDrum_AransasBay, 0, t
    Salinity_AransasBay -> SpottedSeatrout_AransasBay, 0, u
    Salinity_AransasBay -> HardheadCatfish_AransasBay, 0, v
    Salinity_AransasBay -> GafftopsailCatfish_AransasBay, 0, w
    Salinity_AransasBay -> BrownShrimp_AransasBay, 0, x
    BlueCrabSmall_AransasBay -> RedDrum_AransasBay, 1, aa
    Temperature_AransasBay -> RedDrum_AransasBay, 1, bb
    Temperature_AransasBay -> BlueCrabSmall_AransasBay, 1, cc
    BlueCrabSmall_AransasBay -> BlueCrabSmall_AransasBay, 1, dd
    RedDrum_AransasBay -> RedDrum_AransasBay, 1, ee
    BlueCrabSmall_AransasBay -> BlackDrum_AransasBay, 1, ff
    Temperature_AransasBay -> BlackDrum_AransasBay, 1, gg
    BlackDrum_AransasBay -> BlackDrum_AransasBay, 1, hh
    BlueCrabSmall_AransasBay -> SpottedSeatrout_AransasBay, 1, ii
    Temperature_AransasBay -> SpottedSeatrout_AransasBay, 1, jj
    SpottedSeatrout_AransasBay -> SpottedSeatrout_AransasBay, 1, kk
    BlueCrabSmall_AransasBay -> HardheadCatfish_AransasBay, 1, ll
    Temperature_AransasBay -> HardheadCatfish_AransasBay, 1, mm
    HardheadCatfish_AransasBay -> HardheadCatfish_AransasBay, 1, nn
    BlueCrabSmall_AransasBay -> GafftopsailCatfish_AransasBay, 1, oo
    Temperature_AransasBay -> GafftopsailCatfish_AransasBay, 1, pp
    GafftopsailCatfish_AransasBay -> GafftopsailCatfish_AransasBay, 1, qq
    BrownShrimp_AransasBay -> RedDrum_AransasBay, 1, rr
    Temperature_AransasBay -> BrownShrimp_AransasBay, 1, ss
    BrownShrimp_AransasBay -> BrownShrimp_AransasBay, 1, tt
    BrownShrimp_AransasBay -> BlackDrum_AransasBay, 1, uu
    BrownShrimp_AransasBay -> SpottedSeatrout_AransasBay, 1, vv
    BrownShrimp_AransasBay -> HardheadCatfish_AransasBay, 1, ww
    BrownShrimp_AransasBay -> GafftopsailCatfish_AransasBay, 1, xx
    Salinity_AransasBay -> RedDrum_AransasBay, 1, yy
    Salinity_AransasBay -> BlueCrabSmall_AransasBay, 1, zz
    Salinity_AransasBay -> BlackDrum_AransasBay, 1, aaa
    Salinity_AransasBay -> SpottedSeatrout_AransasBay, 1, bbb
    Salinity_AransasBay -> HardheadCatfish_AransasBay, 1, ccc
    Salinity_AransasBay -> GafftopsailCatfish_AransasBay, 1, ddd
    Salinity_AransasBay -> BrownShrimp_AransasBay, 1, eee

    BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 0, abcd
    Temperature_CopanoBay -> RedDrum_CopanoBay, 0, efgh
    Temperature_CopanoBay -> BlueCrabSmall_CopanoBay, 0, ijkl
    BlueCrabSmall_CopanoBay -> BlackDrum_CopanoBay, 0, mnop
    Temperature_CopanoBay -> BlackDrum_CopanoBay, 0, qrst
    BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 0, uvwx
    Temperature_CopanoBay -> SpottedSeatrout_CopanoBay, 0, yzab
    BlueCrabSmall_CopanoBay -> HardheadCatfish_CopanoBay, 0, cdef
    Temperature_CopanoBay -> HardheadCatfish_CopanoBay, 0, ghij
    BlueCrabSmall_CopanoBay -> GafftopsailCatfish_CopanoBay, 0, klmn
    Temperature_CopanoBay -> GafftopsailCatfish_CopanoBay, 0, opqr
    BrownShrimp_CopanoBay -> RedDrum_CopanoBay, 0, stuv
    Temperature_CopanoBay -> BrownShrimp_CopanoBay, 0, wxyz
    BrownShrimp_CopanoBay -> BlackDrum_CopanoBay, 0, abcd
    BrownShrimp_CopanoBay -> SpottedSeatrout_CopanoBay, 0, efgh
    BrownShrimp_CopanoBay -> HardheadCatfish_CopanoBay, 0, ijkl
    BrownShrimp_CopanoBay -> GafftopsailCatfish_CopanoBay, 0, mnop
    Salinity_CopanoBay -> RedDrum_CopanoBay, 0, qrst
    Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 0, uvwx
    Salinity_CopanoBay -> BlackDrum_CopanoBay, 0, yzab
    Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 0, cdef
    Salinity_CopanoBay -> HardheadCatfish_CopanoBay, 0, ghij
    Salinity_CopanoBay -> GafftopsailCatfish_CopanoBay, 0, klmn
    Salinity_CopanoBay -> BrownShrimp_CopanoBay, 0, opqr
    BlueCrabSmall_CopanoBay -> RedDrum_CopanoBay, 1, stuv
    Temperature_CopanoBay -> RedDrum_CopanoBay, 1, wxyz
    Temperature_CopanoBay -> BlueCrabSmall_CopanoBay, 1, abcd
    BlueCrabSmall_CopanoBay -> BlueCrabSmall_CopanoBay, 1, efgh
    RedDrum_CopanoBay -> RedDrum_CopanoBay, 1, ijkl
    BlueCrabSmall_CopanoBay -> BlackDrum_CopanoBay, 1, mnop
    Temperature_CopanoBay -> BlackDrum_CopanoBay, 1, qrst
    BlackDrum_CopanoBay -> BlackDrum_CopanoBay, 1, uvwx
    BlueCrabSmall_CopanoBay -> SpottedSeatrout_CopanoBay, 1, yzab
    Temperature_CopanoBay -> SpottedSeatrout_CopanoBay, 1, cdef
    SpottedSeatrout_CopanoBay -> SpottedSeatrout_CopanoBay, 1, ghij
    BlueCrabSmall_CopanoBay -> HardheadCatfish_CopanoBay, 1, klmn
    Temperature_CopanoBay -> HardheadCatfish_CopanoBay, 1, opqr
    HardheadCatfish_CopanoBay -> HardheadCatfish_CopanoBay, 1, stuv
    BlueCrabSmall_CopanoBay -> GafftopsailCatfish_CopanoBay, 1, wxyz
    Temperature_CopanoBay -> GafftopsailCatfish_CopanoBay, 1, abcd
    GafftopsailCatfish_CopanoBay -> GafftopsailCatfish_CopanoBay, 1, efgh
    BrownShrimp_CopanoBay -> RedDrum_CopanoBay, 1, ijkl
    Temperature_CopanoBay -> BrownShrimp_CopanoBay, 1, mnop
    BrownShrimp_CopanoBay -> BrownShrimp_CopanoBay, 1, qrst
    BrownShrimp_CopanoBay -> BlackDrum_CopanoBay, 1, uvwx
    BrownShrimp_CopanoBay -> SpottedSeatrout_CopanoBay, 1, yzab
    BrownShrimp_CopanoBay -> HardheadCatfish_CopanoBay, 1, cdef
    BrownShrimp_CopanoBay -> GafftopsailCatfish_CopanoBay, 1, ghij
    Salinity_CopanoBay -> RedDrum_CopanoBay, 1, klmn
    Salinity_CopanoBay -> BlueCrabSmall_CopanoBay, 1, opqr
    Salinity_CopanoBay -> BlackDrum_CopanoBay, 1, stuv
    Salinity_CopanoBay -> SpottedSeatrout_CopanoBay, 1, wxyz
    Salinity_CopanoBay -> HardheadCatfish_CopanoBay, 1, abcd
    Salinity_CopanoBay -> GafftopsailCatfish_CopanoBay, 1, efgh
    Salinity_CopanoBay -> BrownShrimp_CopanoBay, 1, ijkl

    BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 0, abcde
    Temperature_MesquiteBay -> RedDrum_MesquiteBay, 0, fghij
    Temperature_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, klmno
    BlueCrabSmall_MesquiteBay -> BlackDrum_MesquiteBay, 0, pqrst
    Temperature_MesquiteBay -> BlackDrum_MesquiteBay, 0, uvwxy
    BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, zabcd
    Temperature_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, efghi
    BlueCrabSmall_MesquiteBay -> HardheadCatfish_MesquiteBay, 0, jklmn
    Temperature_MesquiteBay -> HardheadCatfish_MesquiteBay, 0, opqrs
    BlueCrabSmall_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 0, tuvwz
    Temperature_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 0, abcdf
    BrownShrimp_MesquiteBay -> RedDrum_MesquiteBay, 0, ghijk
    Temperature_MesquiteBay -> BrownShrimp_MesquiteBay, 0, lmnop
    BrownShrimp_MesquiteBay -> BlackDrum_MesquiteBay, 0, qrstu
    BrownShrimp_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, vwxyz
    BrownShrimp_MesquiteBay -> HardheadCatfish_MesquiteBay, 0, abcde
    BrownShrimp_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 0, fghij
    Salinity_MesquiteBay -> RedDrum_MesquiteBay, 0, klmno
    Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 0, pqrst
    Salinity_MesquiteBay -> BlackDrum_MesquiteBay, 0, uvwxy
    Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 0, zabcd
    Salinity_MesquiteBay -> HardheadCatfish_MesquiteBay, 0, efghi
    Salinity_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 0, jklmn
    Salinity_MesquiteBay -> BrownShrimp_MesquiteBay, 0, opqrs
    BlueCrabSmall_MesquiteBay -> RedDrum_MesquiteBay, 1, tuvwz
    Temperature_MesquiteBay -> RedDrum_MesquiteBay, 1, abcdf
    Temperature_MesquiteBay -> BlueCrabSmall_MesquiteBay, 1, ghijk
    BlueCrabSmall_MesquiteBay -> BlueCrabSmall_MesquiteBay, 1, lmnop
    RedDrum_MesquiteBay -> RedDrum_MesquiteBay, 1, qrstu
    BlueCrabSmall_MesquiteBay -> BlackDrum_MesquiteBay, 1, vwxyz
    Temperature_MesquiteBay -> BlackDrum_MesquiteBay, 1, abcde
    BlackDrum_MesquiteBay -> BlackDrum_MesquiteBay, 1, fghij
    BlueCrabSmall_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, klmno
    Temperature_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, pqrst
    SpottedSeatrout_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, uvwxy
    BlueCrabSmall_MesquiteBay -> HardheadCatfish_MesquiteBay, 1, zabcd
    Temperature_MesquiteBay -> HardheadCatfish_MesquiteBay, 1, efghi
    HardheadCatfish_MesquiteBay -> HardheadCatfish_MesquiteBay, 1, jklmn
    BlueCrabSmall_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 1, opqrs
    Temperature_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 1, tuvwz
    GafftopsailCatfish_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 1, abcdf
    BrownShrimp_MesquiteBay -> RedDrum_MesquiteBay, 1, ghijk
    Temperature_MesquiteBay -> BrownShrimp_MesquiteBay, 1, lmnop
    BrownShrimp_MesquiteBay -> BrownShrimp_MesquiteBay, 1, qrstu
    BrownShrimp_MesquiteBay -> BlackDrum_MesquiteBay, 1, vwxyz
    BrownShrimp_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, abcde
    BrownShrimp_MesquiteBay -> HardheadCatfish_MesquiteBay, 1, fghij
    BrownShrimp_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 1, klmno
    Salinity_MesquiteBay -> RedDrum_MesquiteBay, 1, pqrst
    Salinity_MesquiteBay -> BlueCrabSmall_MesquiteBay, 1, uvwxy
    Salinity_MesquiteBay -> BlackDrum_MesquiteBay, 1, zabcd
    Salinity_MesquiteBay -> SpottedSeatrout_MesquiteBay, 1, efghi
    Salinity_MesquiteBay -> HardheadCatfish_MesquiteBay, 1, jklmn
    Salinity_MesquiteBay -> GafftopsailCatfish_MesquiteBay, 1, opqrs
    Salinity_MesquiteBay -> BrownShrimp_MesquiteBay, 1, tuvwz
"

# fit the dsem model for the sciaenid system for aransas bay 
fitAB = dsem( sem = sem1,
              tsdata = AB_Sciaenid_TS,
              family = rep("fixed", ncol(AB_Sciaenid_TS)),             
              estimate_delta0 = TRUE,
              control = dsem_control(
                quiet = TRUE))
summary(fitAB)

# i am not sure what this code is doing
# i think it takes the time series for each site (in our case, three minor bays) and combines their output/fit
# see fig 3 of thorson et al. 2024 "Note that the path coefficient estimates shown here are specified as identical across
# all 12 sites that are treated as separate time series for each species (48 time series total)" 
get_part = function(x){
  vars = c("Salinity","BrownShrimp","GafftopsailCatfish","HardheadCatfish","SpottedSeatrout","Temperature","BlackDrum","RedDrum", "BlueCrabSmall")
  index = sapply( vars, FUN=\(y) grep(y,rownames(x$coef))[1] )
  x$coef = x$coef[index,index]
  dimnames(x$coef) = list( vars, vars )
  return(x)
}

# extract model fit values for and set up path diagrams
p1 = plot( get_part(as_fitted_DAG(fitAB)), type="width", show.legend=FALSE)
p1$layers[[1]]$mapping$edge_width = 0.9
p2 = plot( get_part(as_fitted_DAG(fitAB, what="p_value" )), type="width",
           show.legend=FALSE, colors=c('black', 'black'))
p2$layers[[1]]$mapping$edge_width = 0.1

# arranging  path diagrams of coefficients and p values
ggarrange(p1 + scale_x_continuous(expand = c(0.3, 0)),
          p2 + scale_x_continuous(expand = c(0.3, 0)),
          labels = c("Simultaneous effects", "Two-sided p-value"),
          ncol = 1, nrow = 2)

# results are underwhelming... p value counting issues aside, there is only one significant relationship and no pattern in terms of positive or negative effects
# not shown above, i inspected output/path diagrams for separately for nonlagged and lagged effects and not much changed
# consider having two values per year, a spring mean and fall mean, not just one annual mean (gill net surveys are only done in the spring and fall)