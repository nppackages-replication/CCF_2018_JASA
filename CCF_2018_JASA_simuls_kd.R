rm(list=ls(all=TRUE))

library("nor1mix")
library(Hmisc)
library(nprobust)

setwd ("C:/Users/nsc19/Dropbox/2017/nprobust/simuls/jasa/replications")

n     = 500
kernel= "epa"
alpha = 0.05
p     = 2
deriv = 0
q     = qnorm(1-alpha/2)
rho   = 1 
sim   = 5000
h.seq = seq(0.02,0.8,0.02)
rho.seq = seq(0.2,2,0.2)
eval = c(-2,-1,0,1,2)

model = 1       # select from 1 to 6
c     = eval[1] # select from 1 to 5

# Generate Output (generated in output folder)
source("functions/simuls_kd_grid.R") # to run simulations over bandwidth grid defined by h.seq
source("functions/simuls_kd_bws.R")  # to run simulations using population and data-driven bandwidths
source("functions/simuls_kd_3d.R")   # to run simulations comparing performance at different data-driven bandwidths 

# Generate Plots (generated in plots folder)
source("functions/plots_kd.R")

# Generate Tables (generated in plots folder)
source("functions/tables_kd.R")
