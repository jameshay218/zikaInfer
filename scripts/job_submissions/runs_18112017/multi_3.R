three_states <- c("bahia","pernambuco","riograndedonorte")

microDatFile <- "~/Documents/Zika/Data/brazil/brazil_microCeph_28022017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
microDat <- microDat[microDat$startDay < 1214,]
incDat <- NULL

## Generate peak times
peakTimes <- data.frame(local=three_states,start=c(855,804,862)-30, end = c(855,804,862)+30)

####################
## Report data runs
####################
maxChain <- 5
versions <- c(1,2,3)
priors <- NULL
extra_unfixed <- NULL

combos <- expand.grid(runName = "multi_3", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsA <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                         stateNames=three_states,microDat=microDat,incDat=incDat,preParTab=NULL,
                                         allowablePars=NULL,incStates=NULL,
                                         mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                         useInc=FALSE,usePeakTimes=TRUE,covMat=NULL,
                                         peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                         microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                         stateWeights=NULL,extra_unfixed=extra_unfixed,mosquito_lifespan=5,
                                         do_call=TRUE,timeout=0)
priors <- function(values){
    names(values) <- names
    d_eh <- dnorm(values["D_EH"],5,0.94,1)
    d_em <- dnorm(values["D_EM"], 8.2, 1.12,1)
    d_ih <- dnorm(values["D_IH"], 4.7, 0.51,1)
    
    return(sum(d_eh, d_em, d_ih))  
}
extra_unfixed <- c("D_EH","D_EM","D_IH")
combos <- expand.grid(runName = "multi_3_varied", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsB <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                             stateNames=three_states,microDat=microDat,incDat=incDat,preParTab=NULL,
                             allowablePars=NULL,incStates=NULL,
                             mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                             useInc=FALSE,usePeakTimes=TRUE,covMat=NULL,
                             peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                             microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                             stateWeights=NULL,extra_unfixed=extra_unfixed,do_call=TRUE,timeout=0)

jobs_multi3 <- list(jobsA,jobsB)
