library(coda)
library(deSolve)
library(Rcpp)
library(plyr)
library(reshape2)
library(ggplot2)
library(devtools)
#devtools::install_github("jameshay218/zikaProj")
#library(zikaProj)
setwd("~/Documents/zikaProj")
load_all()
#############################################
# Script to generate some simulation results
#############################################
pars <- setupListPars(10)
pars[[2]][11] <- 10
pars[[2]][1] <- 20000000
pars[[3]][17] <- 50 # Set R0 lower
# R0 = 3.857297
print(r0.calc(pars[[3]],1000000,20000000))

# Mean head circumference at SD for girls by WHO growth standards:
# <http://www.who.int/childgrowth/standards/second_set/hcfa_girls_0_13_zscores.pdf?ua=1>
pars[[3]][5] <- 35
pars[[3]][6] <- 1.1844
pars[[3]][4] <- 1.2 # also set SD for microcephaly infants to 1.2
pars[[3]][1] <- 30
# -2 SD from mean is 31.5 cm, which we will use as diagnosis threshold
threshold <- 32

pars[[1]] <- seq(0,(2+1/12),by=1/365)
pars[[3]][8] <- 1/12
pars[[3]][9] <- 1.25

paramTable <- setupParTable(pars[[3]])
setwd("~/zikaTmp2")
iterations <- 5000
adaptive <- 5000

pernambuco_dat <- read.csv("~/Dropbox/Zika/Data/pernambuco_dat.csv", header=FALSE)

pars1 <- pars[[3]]
pars1[3] <- 25
pars1[4] <- 4
pars1[17] <- 77

y1 <- run_metropolis_MCMC(pars1,iterations,unname(as.matrix(pernambuco_dat)),pars[[1]],pars[[2]],paramTable,0.44,500,1,0,adaptive,"realdat",500,TRUE,threshold)
  