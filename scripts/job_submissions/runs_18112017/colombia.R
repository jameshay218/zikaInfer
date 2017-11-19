microDatFile <- "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2016.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

microDat$local <- "colombia"
incDat$local <- "colombia"

## Generate peak times
peakTimes <- data.frame(local="colombia",start=1128-30, end=1128+30)

maxChain <- 5
versions <- c(1,2,3)
priors <- NULL
extra_unfixed <- NULL

##########


combos <- expand.grid(runName = "colombia_inc", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsA <-  queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                              stateNames="colombia",
                              microDat=microDat,incDat=incDat,preParTab=NULL,
                              allowablePars=NULL,incStates=NULL,
                              mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                              useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                              peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                              microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                              stateWeights=NULL,extra_unfixed=extra_unfixed,
                              do_call=TRUE,timeout=0)


combos <- expand.grid(runName = "colombia_peak", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsB <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                              stateNames="colombia",
                              microDat=microDat,incDat=incDat,preParTab=NULL,
                              allowablePars=NULL,incStates=NULL,
                              mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                              useInc=FALSE,usePeakTimes=TRUE,covMat=NULL,
                              peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                              microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                              stateWeights=NULL,extra_unfixed=extra_unfixed,
                              do_call=TRUE,timeout=0)


combos <- expand.grid(runName = "colombia_inc_month", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_microceph_monthly.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
jobsC <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)

jobs_colombia <- list(jobsA,jobsB,jobsC)
