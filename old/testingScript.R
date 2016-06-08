library(coda)
library(deSolve)
library(Rcpp)
library(plyr)
library(reshape2)
library(ggplot2)
#devtools::install_github("jameshay218/zikaProj")
#library(zikaProj)

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
pars[[3]][5] <- 33.8787
pars[[3]][6] <- 1.1844
pars[[3]][4] <- 1.2 # also set SD for microcephaly infants to 1.2

# -2 SD from mean is 31.5 cm, which we will use as diagnosis threshold
threshold <- 31.5

pars[[1]] <- seq(0,(2+1/12),by=1/365)
pars[[3]][8] <- 1/12
pars[[3]][9] <- 1.25

sampling <- c(0.9, 0.5, 0.1, 0.05, 0.01)
filenames <- c("high","half","low","lowish","verylow")
pars1 <- pars[[3]]
paramTable <- setupParTable(pars[[3]])
setwd("~/zikaTmp2")
iterations <- 100000
adaptive <- 10000
for(i in 1:length(sampling)){
  pars[[3]][2] <- sampling[i]
  dat <- zika.sim(pars)
  
  # Set starting params
  pars1 <- pars[[3]]
  pars1[3] <- 25
  pars1[4] <- 4
  pars1[16] <- 77

  countedDat <- cbind(apply(dat,1,function(x) length(which(x < threshold))),apply(dat,1,length))
 
  filename1 <- paste(filenames[i],"_dat",sep="")
  filename2 <- paste(filenames[i],"_count",sep="")
  y <- run_metropolis_MCMC(pars1,iterations,dat,pars[[1]],pars[[2]],paramTable,0.44,500,1,0,adaptive,filename1,500,TRUE,NULL)
  y1 <- run_metropolis_MCMC(pars1,iterations,countedDat,pars[[1]],pars[[2]],paramTable,0.44,500,1,0,adaptive,filename2,500,TRUE,threshold)
  
  chain1 <- read.csv(y)
  chain1 <- chain1[adaptive:nrow(chain1),c(4,5,17)]
  png(paste(filename1,".png",sep=""))
  plot(as.mcmc(chain1))
  dev.off()
  
  chain2 <- read.csv(y1)
  chain2 <- chain2[adaptive:nrow(chain2),c(4,5,17)]
  png(paste(filename2,".png",sep=""))
  plot(as.mcmc(chain2))
  dev.off()
}
###############################################################################
  