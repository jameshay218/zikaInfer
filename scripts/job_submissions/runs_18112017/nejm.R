local <- "northeast"
microDatFile <- "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
microDat$local <- "northeast"
incDatFile <- "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
incDat$local <- "northeast"
incDat <- incDat[incDat$startDay < 1040,]
versions <- c(1,2,3)
maxChain <- 5
priors <- NULL
extra_unfixed <- NULL
combos <- expand.grid(runName = "northeast", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobs_nejm <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                  stateNames="northeast",
                                  microDat=microDat,incDat=incDat,preParTab=NULL,
                                  allowablePars=NULL,incStates=NULL,
                                  mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                  useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                  peakTimeRange=120,peakTime=855,sim=FALSE, predict=FALSE,
                                  microChain=NULL,prePeakTimes=NULL,normLik=FALSE,
                                  stateWeights=NULL,extra_unfixed=extra_unfixed,
                                  do_call=TRUE,timeout=0)




local <- "northeast"
microDatFile <- "~/Documents/Zika/Data/brazil/deOliveira2017/microCeph_deOliveira2017_clean.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
microDat$local <- "northeast"
incDatFile <- "~/Documents/Zika/Data/brazil/deOliveira2017/zikv_inc_deOliveira2017_clean.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
incDat$local <- "northeast"
incDat <- incDat[incDat$startDay < 1034,]
versions <- c(2,3)
maxChain <- 5
priors <- NULL
extra_unfixed <- NULL
combos <- expand.grid(runName = "northeast_deOliveira2017", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobs_deOliveira <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                        stateNames="northeast",
                                        microDat=microDat,incDat=incDat,preParTab=NULL,
                                        allowablePars=NULL,incStates=NULL,
                                        mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                        useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                        peakTimeRange=120,peakTime=855,sim=FALSE, predict=FALSE,
                                        microChain=NULL,prePeakTimes=NULL,normLik=FALSE,
                                        stateWeights=NULL,extra_unfixed=extra_unfixed,
                                        do_call=TRUE,timeout=0)


