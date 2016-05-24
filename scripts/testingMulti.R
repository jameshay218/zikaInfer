library(zikaProj)
library(devtools)
library(deSolve)
library(MASS)
library(ggplot2)

setwd("~/Documents/zikaProj")
load_all()


realDat <- read.csv("~/Dropbox/Zika/Data/allDat20.04.16.csv")
parTab <- setupParTable(version=3,realDat=realDat)

#places <- c("all","pernambuco","bahia","saopaulo")

states_with_data <- c(
                      "all",
                      "pernambuco",
                      "bahia",
                      "saopaulo",
                      "paraiba",
                      "maranhao",
                      "ceara",
                      "sergipe",
                      "riodejaneiro",
                      "piaui",
                      "riograndedonorte",
                      "minasgerais",
                      "matogrosso",
                      "alagoas",
                      "para",
                      "acre",
                      "goias",
                      "espiritosanto",
                      "tocantins")

states_with_data <- c("all","pernambuco","bahia")

#parTab <- parTab[parTab$local %in% places,]
parTab <- parTab[parTab$local %in% states_with_data,]
parTab[parTab$names=="b","fixed"] <- 1

#parTab[parTab$fixed == 0,"steps"] <- c(0.0003, 0.005, 0.06, 0.0005, 0.0005, 0.01, 0.0005, 0.01, 0.00025, 0.008)

#pars <- setupListPars(version=3)
#tmpVal <- parTab$values
#parTab$values <- testDatPars

testDat <- generate_multiple_data(pars[[1]],parTab,NULL)
#parTab$values <- tmpVal
adaptive_period <- 50000
burnin <- 20000

places <- states_with_data[states_with_data != "all"]
realDat <- realDat[realDat$local %in% places,4:9]
#places <- places[places != "all"]

#peakTimes <- data.frame("start"=c(280,350,225),"end"=c(320,390,250),"local"=c("pernambuco","bahia","saopaulo"))
#peakTimes$local <- as.character(peakTimes$local)

peakTimes <- data.frame("start"=450,"end"=560,"local"=as.character(unique(realDat$local)))
peakTimes$local <- as.character(peakTimes$local)
realDat <- realDat[,colnames(testDat)]
parTab[parTab$names=="L_H","values"] <- parTab[parTab$names=="L_H","values"]*365
parTab[parTab$names=="epiStart","values"] <- 360
parTab[parTab$names=="epiStart","lower_bounds"] <- 100
parTab[parTab$names=="epiStart","upper_bounds"] <- 400
parTab[parTab$names=="density","lower_bounds"] <- 1.61
parTab[parTab$names=="density","values"] <- 8

parTab[parTab$names=="mean","fixed"] <- 0
parTab[parTab$names=="var","fixed"] <- 0
parTab[parTab$names=="var","values"] <- 20

setwd("~/Documents/zikaProj/results/intial2")

omg <- run_metropolis_MCMC(
        iterations=50000,
        data=realDat,
        t_pars=pars[[1]],
        y0s=NULL,
        N_H=10000,
        N_M=30000,
        version=3,
        param_table=parTab,
        opt_freq=1000,
        thin=1,
        burnin=burnin,
        adaptive_period=adaptive_period,
        filename="real",
        save_block=500,
        threshold=NULL,
        buckets=NULL,
        mvrPars=NULL,
        incDat=NULL,
        allPriors=NULL,
        peakTimes=peakTimes,
        VERBOSE=FALSE
)

chain <- read.csv("real_chain.csv")
chain <- chain[(adaptive_period+burnin):nrow(chain),]
covMat <- cov(chain[,2:(ncol(chain)-1)])      
startPars <- as.numeric(chain[which.max(chain[,ncol(chain)]),2:(ncol(chain)-1)])
names(startPars) <- parTab$names
mvrPars <- list(covMat,0.8,1)

parTab$values <- startPars

omg <- run_metropolis_MCMC(
  iterations=250000,
  data=realDat,
  t_pars=pars[[1]],
  y0s=NULL,
  N_H=10000,
  N_M=30000,
  version=3,
  param_table=parTab,
  popt=0.234,
  opt_freq=1000,
  thin=1,
  burnin=burnin,
  adaptive_period=adaptive_period,
  filename="real1",
  save_block=500,
  threshold=NULL,
  buckets=NULL,
  mvrPars=mvrPars,
  incDat=NULL,
  allPriors=NULL,
  peakTimes=peakTimes,
  VERBOSE=FALSE
)

chain <- read.csv("real1_chain.csv")
chain <- chain[(adaptive_period+burnin):nrow(chain),]
#chain <- chain[!is.na(chain),]
indices <- which(parTab$fixed == 0)
pdf("traces.pdf")
plot(as.mcmc(chain[,c(indices+1,ncol(chain))]))
dev.off()


