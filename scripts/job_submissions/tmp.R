
peakTimes <- data.frame(local=c("bahia","pernambuco","riograndedonorte","pernambuco_confirmed","riograndedonorte_confirmed"),start=c(855,804,862,804,862)-30, end = c(855,804,862,804,862)+30)


stateNames="riograndedonorte"
microDat=microDat
incDat=NULL
preParTab=NULL
allowablePars=NULL
incStates=NULL
mcmcPars=mcmcPars
mcmcPars2=mcmcPars2
allPriors=priors
useInc=FALSE
usePeakTimes=TRUE
covMat=NULL
peakTimeRange=NULL
peakTime=NULL
sim=FALSE
predict=FALSE
microChain=NULL
prePeakTimes=peakTimes
normLik=FALSE
stateWeights=NULL
extra_unfixed=extra_unfixed
sharedProb = TRUE
mosquito_lifespan=5
version=1
chainNo =1
runName="test"

