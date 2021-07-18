rm(list=ls(all=TRUE))

library(Hmisc)
library(Rcpp)
library(locfit)
library(numDeriv)
library(nprobust)
library(locpol)
library(KernSmooth)
library(locfit)
library(ggplot2)
library(reshape)
library(grid)


n      = 500
sim    = 5000
kernel = "epa"
p      = 1
q      = 2
deriv  = 0
vce    = "hc3"
level  = 5
qz     = qnorm((100-level/2)/100)
h.seq  = seq(0.1,0.7,0.02)
eval   = c(-2/3,-1/3,0,1/3,2/3)

model = 1       # select from 1 to 6
c     = eval[1] # select from 1 to 5

# Generate Output (generated in output folder)
source("functions/simuls_lp_grid.R") # to run simulations over bandwidth grid defined by h.seq
source("functions/simuls_lp_bws.R")  # to run simulations using population and data-driven bandwidths
source("functions/simuls_lp_comp.R") # to run simulations comparing performance at different data-driven bandwidths 
source("functions/simuls_lp_vce.R")  # to run simulations comparing performance at different vce choices

# Generate Plots (generated in plots folder)
source("functions/plots_lp.R")

# Generate Tables (generated in plots folder)
source("functions/tables_lp.R")



