#######################################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## Check all MCMC chains to see if they have converged based on gelman diagnostics and ESS from the coda package.
## Just specify the working directory with all of the chains at the top, and edit the save location at the bottom
## of this script
library(data.table)
library(lazymcmc)
library(zikaProj)
library(coda)

## Location of all the MCMC chains
chainwd <- "~/net/home/zika/outputs/"
setwd(chainwd)

minESS <- NULL
minESSName <- NULL
maxGelman <- NULL
maxGelmanName <- NULL
mpsrf <- NULL
reruns <- NULL
allModels <- NULL
runNames <- NULL
allStates <- NULL

## Recursively get all directories in this location - make full directory locations
allDirs <- list.dirs()
allDirs <- allDirs[allDirs != "."]
allDirs <- paste0(chainwd, allDirs)

## Function to read in MCMC chains if possible and then calculate convergence diagnostics
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

## Go through each directory
for(dir in allDirs){
  print(dir)
  ## If this is a forecasting analysis, run slightly differently
    if(length(grep("forecast",dir)) == 0){
      parTab <- read_inipars(dir)
      res <- read_and_check(dir, parTab, 750000)
    } else {
      parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv")
      res <- read_and_check(dir, parTab, 50000)
    }
    minESS <- c(minESS, res[[1]])
    minESSName <- c(minESSName, res[[2]])
    
    maxGelman <- c(maxGelman, res[[3]])
    maxGelmanName <- c(maxGelmanName, res[[4]])
    mpsrf <- c(mpsrf, res[[5]])
    reruns <- c(reruns, res[[6]])
}

convergence <- data.frame(allDirs, minESS,minESSName,maxGelman,maxGelmanName,mpsrf,reruns)
colnames(convergence) <- c("runName","minESS","whichMinESS","maxGelman","whichMaxGelman","mpsrf","rerun")
write.table(convergence, "~/Documents/Zika/convergenceCheck.csv",sep=",",row.names=FALSE)
