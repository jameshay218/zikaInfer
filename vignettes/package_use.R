## ----setup, echo=FALSE---------------------------------------------------
library(zikaProj)
library(lazymcmc)
library(kableExtra)
library(knitr)
options(knitr.table.format = "html")
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

## ------------------------------------------------------------------------
data(exampleParTab)

## Extract parameter values to solve model. Note that 
## these parameters MUST be named
pars <- exampleParTab$values
names(pars) <- exampleParTab$names

## ------------------------------------------------------------------------
## Generate starting population sizes based on human population size
## and mosquito density
y0s <- generate_y0s(pars["N_H"], pars["density"],iniI=10)
## Times to solve model over. Note that the unit of time is in days
ts <- seq(0,3003,by=1)
y <- solveSEIRModel_rlsoda(ts, y0s, pars, TRUE)
plot(y[,"I_H"],type='l',col="red",
     xlab="Time (days)", ylab="ZIKV incidence in humans")

