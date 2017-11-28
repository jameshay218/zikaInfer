## ----setup, echo=FALSE, message=FALSE, warning=FALSE---------------------
library(zikaProj)
library(lazymcmc)
library(kableExtra)
library(knitr)
library(zoo)
options(knitr.table.format = "html")
knitr::opts_chunk$set(fig.width=6, fig.height=4) 

## ----message=FALSE, warning=FALSE----------------------------------------
data(exampleParTab)

## Extract parameter values to solve model. Note that 
## these parameters MUST be named
pars <- exampleParTab$values
names(pars) <- exampleParTab$names

## ----message=FALSE, warning=FALSE----------------------------------------
## Generate starting population sizes based on human population size
## and mosquito density
y0s <- generate_y0s(pars["N_H"], pars["density"],iniI=10)
## Times to solve model over. Note that the unit of time is in days
ts <- seq(0,3003,by=1)
y <- solveSEIRModel_rlsoda(ts, y0s, pars, TRUE)
plot(y[,"I_H"],type='l',col="red",
     xlab="Time (days)", ylab="ZIKV incidence in humans")

## ----message=FALSE, warning=FALSE----------------------------------------
## Calculate R0
print(r0.calc(pars))

## Find the mosquito density required for a desired R0
print(density.calc(pars, R0=3))

print(calculate_AR(r0=3))

## ------------------------------------------------------------------------
risk <- generate_micro_curve(pars)
plot(risk,type='l',col="blue",
     xlab="Gestational time at infection (days)",ylab="Prob of congenital abnormality")

## ----message=FALSE, warning=FALSE----------------------------------------
probM <- generate_probM(y[,"I_M"],pars["N_H"],risk,pars["b"],pars["p_MH"],pars["baselineProb"],1)
plot(probM,type='l',col="green",
     xlab="Time (days)",ylab="Proportion of microcephaly affected births")

## ----message=FALSE, warning=FALSE----------------------------------------
library(ggplot2)

## Generate simulated data with known parameters
simDat <- generate_multiple_data(ts, parTab=exampleParTab, weeks=FALSE, 
                                  dataRangeMicro=c(600,2000),dataRangeInc=c(600,1500), 
                                  noise=FALSE,peakTimeRange=60)
microDat <- simDat[[1]]
incDat <- simDat[[2]]
peakTimes <- simDat[[3]]

## Save simulated data
write.table(microDat,"sim_microDat.csv",sep=",",row.names=FALSE)
write.table(incDat,"sim_incDat.csv",sep=",",row.names=FALSE)

## Check the generated data
ggplot(incDat) + 
  geom_line(aes(x=startDay,y=inc/N_H), col="red") +
  geom_line(data=microDat,aes(x=startDay,y=microCeph*0.001/births),col="blue") + 
  facet_wrap(~local) +
  ylab("Per capita ZIKV incidence (red)") +
  xlab("Time (days)") +
  scale_y_continuous(sec.axis=sec_axis(~.*1000,name="Per birth microcephaly\n incidence (blue)")) +
  theme_bw()

## ----  message=FALSE, warning=FALSE--------------------------------------
library(lazymcmc)
## Generate random starting points for MCMC chain
## Note that we have to constrain the starting points for R0 and the epidemic seed time
## such that trajectory is near expected peak time
startTab <- generateStartingParTab(exampleParTab, peakTimes, restrictedR0=TRUE,"")

## MCMC control parameters
mcmcPars <- c("adaptive_period"=20000,"iterations"=50000,"opt_freq"=1000,"thin"=10,"save_block"=100,"popt"=0.44)

## Run MCMC chain using univariate sampler
result <- lazymcmc::run_MCMC(parTab=startTab, data=microDat, mcmcPars=mcmcPars,filename="test_univariate",
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
## Run two chains to calculate gelman diagnostics
final <- lazymcmc::run_MCMC(parTab=startTab, data=microDat, mcmcPars=mcmcPars,filename="test_2_multivariate",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=mvrPars,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=NULL)
finalB <- lazymcmc::run_MCMC(parTab=startTab, data=microDat, mcmcPars=mcmcPars,filename="test_1_multivariate",
                    CREATE_POSTERIOR_FUNC = create_posterior, mvrPars=mvrPars,PRIOR_FUNC=NULL,
                    OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=NULL)

## ----message=FALSE, warning=FALSE----------------------------------------
library(data.table)
library(coda)
# Function to read in MCMC chains if possible and then calculate convergence diagnostics
read_and_check <- function(wd, parTab, burnin){
  chain <- lazymcmc::load_mcmc_chains(wd, parTab, TRUE,1,burnin, TRUE,FALSE,FALSE)
  if(length(chain[[1]]) > 1){
    ess <- ess_diagnostics(chain[[1]],200)
    gelman <- gelman_diagnostics(chain[[1]], 1.1)
    minESS <- min(ess$ESS)
    whichMinESS <- names(which.min(ess$ESS))
    maxGelman <- gelman$WorstGelman[1]
    maxGelmanName <- gelman$WorstGelman[2]
    mpsrf <- gelman$WorstGelman[3]
    rerun <- gelman$Rerun | minESS < 200
  } else {
    minESS <- NA
    whichMinESS <- NA
    maxGelman <- NA
    maxGelmanName <- NA
    mpsrf <- NA
    rerun <- TRUE
  }
  return(list(minESS,whichMinESS,maxGelman,maxGelmanName,mpsrf,rerun))
}

dir <- getwd() ## The working directory that the MCMC chain was run
## The parameter table used in this run will 
## be saved in the working directory
res <- read_and_check(dir, startTab, mcmcPars["adaptive_period"])
print(res)

## ------------------------------------------------------------------------
plot_random_microceph_curves(chain,100)
indiv_model_fit("sim_microDat.csv","sim_incDat.csv", "bahia","Simulated data",0.001,200,0.03, 500, FALSE, TRUE, startTab, chain, FALSE, FALSE)

