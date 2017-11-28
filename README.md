# zikaProj
> Infer teratogenic congenital syndrome risk from transmission and incidence data

## Introduction
This repository contains code and methodology to allow the user to estimate a gestational risk profile for congenital diseases caused by teratogenic pathogens (namely microcephaly caused by Zika virus (ZIKV) infection). The full model description can be found in the accompanying vignette here:

[Model methodology](https://jameshay218.github.io/driftSim/inst/doc/methods.html)

The package itself can be split into a number of sections:

1. Functions required to solve the transmission (SEIR) model
2. Functions required to generate the gestational-week-varying risk curve
3. Functions to calculate the posterior probability of model parameters given a set of data and an MCMC algorithm to estimate posterior distributions
4. Functions for post-MCMC statistics and analysis
5. Plotting functions

A full working R script is provided below, but users should refer to the [extended usage](https://jameshay218.github.io/driftSim/inst/doc/package_use.html) vignette for a full breakdown of the work flow.

## Installation
Installation should be straightforward, though there are two dependencies on non-CRAN packages which can be found here: 

[`rlsoda`](https://github.com/richfitz/rlsoda)
[`lazymcmc`](https://github.com/jameshay218/lazymcmc)

Note that the [`rlsoda`](https://github.com/richfitz/rlsoda) package is just used to speed up solving the ODE model. This could otherwise be solved used the `lsoda` package in base R. Users wishing to change this can set the solver using the `solver` argument when calling the posterior function (see `?create_posterior`).

To install the package itself, simply type the following:
```r
devtools::install_github("jameshay218/zikaProj")
```

## Usage - simulated data
The script below can be used to generate simulated data and re-estimate the process-generating parameters. This script is included to assure the reader that the inference framework works correctly.
```r
library(zikaProj)

## Read in some example parameters
## A table like this is needed to control prior bounds and which parameters are 
## fixed during MCMC.
data(exampleParTab)

## Times to solve model over
ts <- seq(0,3003,by=1)

## Generate simulated data with known parameters
simDat <- generate_multiple_data(ts, parTab=exampleParTab, weeks=FALSE, 
                                  dataRangeMicro=c(600,2000),dataRangeInc=c(600,1500), 
                                  noise=FALSE,peakTimeRange=60)
          
## Extract incidence and microcephaly data
microDat <- simDat[[1]]
incDat <- simDat[[2]]
peakTimes <- simDat[[3]]

## Generate random starting points for MCMC chain
## Note that we have to constrain the starting points for R0 and the epidemic seed time
## such that trajectory is near expected peak time
startTab <- generateStartingParTab(exampleParTab, peakTimes, restrictedR0=TRUE,"")

## MCMC control parameters
mcmcPars <- c("adaptive_period"=20000,"iterations"=50000,"opt_freq"=1000,"thin"=10,"save_block"=100,"popt"=0.44)

## Run MCMC chain using univariate sampler
result <- lazymcmc::run_MCMC(parTab=startTab, data=microDat, mcmcPars=mcmcPars,filename="test",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=NULL,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=NULL)

## Use these results to get covariance matrix to run multivariate sampler
chain <- read.csv(result$file)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],]
covMat <- cov(chain[,2:(ncol(chain)-1)])
startTab$values <- get_best_pars(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
mcmcPars["popt"] <- 0.234

## Run MCMC chain using multivariate sampler
final <- lazymcmc::run_MCMC(parTab=startTab, data=microDat, mcmcPars=mcmcPars,filename="test",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=mvrPars,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=NULL)
                   
chain <- read.csv(final$file)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],(which(exampleParTab$fixed == 0) + 1)]

## Compare posterior densities to simulated point values
pars <- exampleParTab$values
names(pars) <- exampleParTab$names
print(pars)

plot(coda::as.mcmc(chain))
```

## Usage - real data
The script below can be used to estimate posterior densities for data from Bahia, Brazil.
```r
library(zikaProj)

## Read in some example parameters
## A table like this is needed to control prior bounds and which parameters are 
## fixed during MCMC.
data(exampleParTab)

## Times to solve model over
ts <- seq(0,3003,by=1)

## Read in data
## Note that "2013-01-01" is considered to be day 0
data(bahiaMicroDat)
data(bahiaIncDat)
data(bahiaPeakTimes)

## Generate random starting points for MCMC chain
## Note that we have to constrain the starting points for R0 and the epidemic seed time
## such that trajectory is near expected peak time
startTab <- generateStartingParTab(exampleParTab, bahiaPeakTimes, restrictedR0=TRUE,"")

## MCMC control parameters
mcmcPars <- c("adaptive_period"=20000,"iterations"=50000,"opt_freq"=1000,"thin"=10,"save_block"=100,"popt"=0.44)

## Run MCMC chain using univariate sampler
result <- lazymcmc::run_MCMC(parTab=startTab, data=bahiaMicroDat, mcmcPars=mcmcPars,filename="test",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=NULL,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=bahiaIncDat,peakTimes=NULL)

## Use these results to get covariance matrix to run multivariate sampler
chain <- read.csv(result$file)
chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],]
covMat <- cov(chain[,2:(ncol(chain)-1)])
startTab$values <- get_best_pars(chain)
mvrPars <- list(covMat,2.38/sqrt(nrow(startTab[startTab$fixed==0,])),w=0.8)
mcmcPars["popt"] <- 0.234

## Run MCMC chain using multivariate sampler
final <- lazymcmc::run_MCMC(parTab=startTab, data=bahiaMicroDat, mcmcPars=mcmcPars,filename="test",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=mvrPars,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=bahiaIncDat,peakTimes=NULL)
                   
chain <- read.csv(final$file)

## Look at inferred risk curve
plot_random_microceph_curves(chain,100)

chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],(which(exampleParTab$fixed == 0) + 1)]
plot(coda::as.mcmc(chain))

```
## License

GPL-2 © [James Hay &lt;james.hay13@imperial.ac.uk&gt;](https://github.com/).
