realDat <- read.csv("~/Dropbox/Zika/Data/allDat20.04.16.csv")
parTab <- setupParTable(version=3,realDat=realDat)

places <- c("all","pernambuco","bahia","saopaulo")

testDatPars <- c(2.000000e-05, 0.000000e+00, 1.400000e+01, 4.000000e+00, 4.000000e+00, 
                 5.000000e+00, 2.500000e-01, 5.000000e-01, 5.000000e-01, 0.000000e+00, 1.200000e+01,
                 1.000000e+01, 1.000000e+00, 7.000000e+00, 7.430000e+01, 9.277727e+06, 4.100000e+00,
                 2.500000e+01, 7.310000e+01, 1.512637e+07, 3.200000e+00, 5.000000e+00, 7.890000e+01,
                 4.403530e+07, 7.000000e+00, 6.000000e+01)

parTab <- parTab[parTab$local %in% places,]
parTab[parTab$names=="b","fixed"] <- 1

parTab[parTab$fixed == 0,"steps"] <- c(0.0003, 0.005, 0.06, 0.0005, 0.0005, 0.01, 0.0005, 0.01, 0.00025, 0.008)

pars <- setupListPars(version=3)
tmpVal <- parTab$values
parTab$values <- testDatPars

testDat <- generate_multiple_data(pars[[1]],parTab,NULL)
parTab$values <- tmpVal
adaptive_period <- 50000
burnin <- 10000

places <- places[places != "all"]

peakTimes <- data.frame("start"=c(280,350,225),"end"=c(320,390,250),"local"=c("pernambuco","bahia","saopaulo"))
peakTimes$local <- as.character(peakTimes$local)

omg <- run_metropolis_MCMC(
        iterations=100000,
        data=testDat,
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
        filename="test",
        save_block=500,
        threshold=NULL,
        buckets=NULL,
        mvrPars=NULL,
        incDat=NULL,
        allPriors=NULL,
        peakTimes=peakTimes,
        VERBOSE=FALSE
)

chain <- read.csv(omg)
chain <- chain[adaptive_period+burnin:nrow(chain),]
covMat <- cov(greb[,2:(ncol(greb)-1)])      
startPars <- as.numeric(chain[which.max(chain[,ncol(chain)]),2:(ncol(chain)-1)])
names(startPars) <- parTab$names
mvrPars <- list(covMat,0.8,1)

parTab$values <- startPars

omg <- run_metropolis_MCMC(
  iterations=100000,
  data=testDat,
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
  filename="test1",
  save_block=500,
  threshold=NULL,
  buckets=NULL,
  mvrPars=mvrPars,
  incDat=NULL,
  allPriors=NULL,
  peakTimes=peakTimes,
  VERBOSE=FALSE
)

chain <- read.csv(omg)
chain <- chain[(adaptive_period+burnin):nrow(chain),]
#chain <- chain[!is.na(chain),]
indices <- which(parTab$fixed == 0)
plot(as.mcmc(chain[,c(indices+1,ncol(chain))]))