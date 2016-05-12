library(coda)
library(deSolve)
library(Rcpp)
library(plyr)
library(reshape2)
library(ggplot2)
#devtools::install_github("jameshay218/zikaProj")
library(zikaProj)
library(doParallel)
library(doMC)
library(MASS)
registerDoMC(6)
## Number of mosquitoes per person
densityM <- 3
iterations <- 50000
adaptive <- 10000

setwd("~/tmpZika4")
allDat <- read.csv("~/Dropbox/Zika/Data/allDat.csv")
places <- unique(allDat$local)
#places <- places[1:2]

endI <- length(places)
places <- places[1:2] 

for(i in 1:endI){
  writeLines(c(""),paste("log",i,".txt",sep=""))
}
#####
covMat <- matrix(rep(1,nrow(paramTable)^2),nrow(paramTable),nrow(paramTable))
mvrPars <- list(covMat,0.8,1)
#####


for(i in 1:length(places)){
#results <- foreach(i=1:endI, .packages=c("zikaProj","coda","deSolve","Rcpp","plyr","reshape2") ) %dopar% {
 #sink(paste("log",i,".txt",sep=""),append=TRUE)
  filename <- as.character(places[i])
#  cat(paste("Run: ",i,sep=""))
  #'  Setup data
  tmpDat <- allDat[allDat$local == places[i],]
  start <- seq(tmpDat[1,"start"],(nrow(tmpDat)/12) - 1/12,by=1/12)
  end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
  buckets <- cbind(start,end)
  N_H <- tmpDat[1,"NH"]
  L_H <- tmpDat[1,"LH"]
  tmpDat <- tmpDat[,c("microCeph","births")]
  tmpDat$births <- tmpDat$births - tmpDat$microCeph
  
  pars <- setupListPars(duration=end[length(end)], N_H = N_H, N_M = N_H*densityM)  
  pars[[1]][1] <- 4
  pars[[3]]["sampFreq"] <- 30
  pars[[3]]["probMicro"] <- 0.2
  pars[[3]]["baselineProb"] <- 0.01
  pars[[3]]["epiStart"] <- 1.2
  pars[[3]]["L_H"] <- L_H
  pars[[3]]["constSeed"] <- 0
  
  paramTable <- setupParTable(pars[[3]])
  paramTable[c(3,4),"fixed"] <- 1
  paramTable[c(7,8,10,18),"fixed"] <- 0
  paramTable[10,"upper_bounds"]
  paramTable[10,"upper_bounds"] <- 3.9
  paramTable[10,"lower_bounds"] <- pars[[1]][2]
  
  startPars <- pars[[3]]
  startPars["probMicro"] <- runif(1,0,0.2)
  startPars["baselineProb"] <- runif(1,0,0.1)
  startPars["b"] <- runif(1,25,125)
  startPars["epiStart"] <- runif(1,0.8,1.2)
  
  tryCatch(y <- run_metropolis_MCMC(startPars,iterations,tmpDat,pars[[1]],pars[[2]],paramTable,0.44,500,1,0,adaptive,filename,500,FALSE,NULL, buckets,mvrPars),finally=print("Wow"))
  
}
allStats <- NULL
for(i in 1:length(places)){
  setwd("~/tmpZika3")
  y <- paste(places[i],"_chain.csv",sep="")
  filename <- places[i]
  chain <- read.csv(y)
  chain <- chain[adaptive:nrow(chain),c("probMicro","baselineProb","epiStart","b","r0","lnlike")]
  setwd("~/Dropbox/Zika/results2")
  pdf(paste(filename,".pdf",sep=""))
  plot(as.mcmc(chain))
  dev.off()
  
  stats <- summary(as.mcmc(chain))
  toSave <- cbind(stats$statistics, stats$quantiles)
  toSave <- as.data.frame(toSave)
  toSave <- cbind(toSave,"local"=rep(as.character(places[i]),nrow(toSave)))
  
  allStats[[i]] <- toSave
}
setwd("~/Dropbox/Zika/results2")
finalStats <- do.call("rbind",allStats)
write.csv(finalStats,"finalstats2.csv")