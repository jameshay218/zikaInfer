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
iterations <- 20000
adaptive <- 5000

pernambuco_dat <- read.csv("~/Dropbox/Zika/Data/pernambuco_dat.csv", header=FALSE)

pars1 <- pars[[3]]
pars1[3] <- 25
pars1[4] <- 4
pars1[17] <- 77

y1 <- run_metropolis_MCMC(pars1,iterations,greb,pars[[1]],pars[[2]],paramTable,0.44,500,1,0,adaptive,"realdat",500,TRUE,NULL,NULL)
  
p <- NULL
alphas <- cbind(alphas,1-alphas)
for(i in 1:nrow(alphas)){
  p[i] <- alphas[i,1]*pnorm(threshold,mus[1],sds[1],1,0) + alphas[i,2]*pnorm(threshold,mus[2],sds[2],1,0)
}

testing <- seq(0,1.5,by=1/120)
tmp <- NULL
tmpAlphas <- NULL
for(i in 1:length(testing)){
  pars1 <- pars[[3]]
  pars1[9] <- testing[i]
  y <- solveModel(list(pars[[1]],pars[[2]],pars1))
  tmpAlphas[[i]] <- alphas <- calculate_alphas_buckets(
    as.matrix(unname(y[,c("times","If","Sf","Ef","Rf")])),
    pars1[2],
    as.matrix(unname(dat[,c("start","end")])))
  tmp[i] <- likelihood_threshold(
    unname(as.matrix(dat[,c("microCeph","births")])),
    unname(cbind(alphas,1-alphas)),c(pars1[3],pars1[5]),
    c(pars1[4],pars1[6]),
  32)
}

pars <- setupListPars(duration=2+1/12,N_H=10000000,N_M = 200000000)
pars[[3]][1] <- 30
pars[[2]][11] <- 100
pars[[3]][8] <- 3
pars[[3]][9] <- 1.4
pars[[3]][5] <- 35
pars[[3]][6] <- 1.18
pars[[3]][7] <- 0.2
dat <- zika.sim(pars)

testing <- seq(0,200,by=1)
tmp <- NULL
countedDat <- cbind(apply(dat,1,function(x) length(which(x < threshold))),apply(dat,1,length))
countedDat <- as.data.frame(countedDat)
colnames(countedDat) <- c("microCeph","births")
countedDat$start <- seq(0,pars[[1]][1]-pars[[3]][1]/365,by=pars[[3]][1]/365)
countedDat$end <- seq(pars[[3]][1]/365,pars[[1]][1],by=pars[[3]][1]/365)

testing <- seq(0,0.2,by=0.001)
tmp <- NULL
for(i in 1:length(testing)){
print(i)
 pars1 <- pars[[3]]
 pars1["baselineProb"] <- testing[i]
 tmp[i] <- posterior(pars[[1]],pars[[2]],pars1,greb,NULL,NULL)
}