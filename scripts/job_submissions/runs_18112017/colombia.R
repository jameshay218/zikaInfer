
## Generate peak times
peakTimes <- data.frame(local="colombia",start=1128-30, end=1128+30)

maxChain <- 5
versions <- c(1,2)
priors <- NULL
extra_unfixed <- NULL

##########
## PAHO suspected data
microDatFile <- "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly_2017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
microDat$local <- "colombia"
incDat$local <- "colombia"

combos <- expand.grid(runName = "colombia_inc_suspected", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
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
## PAHO without incidence
microDatFile <- "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly_2017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
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

##########
## PAHO confirmed data
microDatFile <- "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly_2017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
microDat$local <- "colombia"
incDat$local <- "colombia"

combos <- expand.grid(runName = "colombia_inc_confirmed", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsC <-  queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                              stateNames="colombia",
                              microDat=microDat,incDat=incDat,preParTab=NULL,
                              allowablePars=NULL,incStates=NULL,
                              mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                              useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                              peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                              microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                              stateWeights=NULL,extra_unfixed=extra_unfixed,
                              do_call=TRUE,timeout=0)


## Cuevas microcephaly, PAHO suspected ZIKV
combos <- expand.grid(runName = "colombia_inc_month", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_microceph_monthly.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
microDat$local <- "colombia"
incDat$local <- "colombia"
jobsD <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)


## Cuevas microcephaly, PAHO confirmed ZIKV
combos <- expand.grid(runName = "colombia_inc_month_confirmed", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_microceph_monthly.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
microDat$local <- "colombia"
incDat$local <- "colombia"
jobsE <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)




## Epidemiological report suspected microcephaly, PAHO suspected ZIKV
combos <- expand.grid(runName = "colombia_epi_reports_suspected_zikv_suspected", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_reports_microceph_2017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)
microDat <- microDat[2:nrow(microDat),]
microDat$local <- "colombia"
incDat$local <- "colombia"
jobsF <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)



## Epidemiological report confirmed microcephaly, PAHO suspected ZIKV
combos <- expand.grid(runName = "colombia_epi_reports_micro_conf_inc_suspected", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_reports_microceph_2017_confirmed.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

microDat$local <- "colombia"
incDat$local <- "colombia"
jobsG <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)

## Epidemiological report confirmed microcephaly, PAHO confirmed ZIKV
combos <- expand.grid(runName = "colombia_epi_reports_micro_conf_inc_confirmed", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_reports_microceph_2017_confirmed.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

microDat$local <- "colombia"
incDat$local <- "colombia"
jobsH <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)


## Epidemiological report suspected microcephaly, PAHO confirmed ZIKV
combos <- expand.grid(runName = "colombia_epi_reports_micro_suspected_inc_confirmed", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
microDatFile <- "~/Documents/Zika/Data/colombia/colombia_reports_microceph_2017.csv"
microDat <- read.csv(microDatFile,stringsAsFactors=FALSE)
incDatFile <- "~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv"
incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

microDat$local <- "colombia"
incDat$local <- "colombia"
jobsI <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                   stateNames="colombia",
                                   microDat=microDat,incDat=incDat,preParTab=NULL,
                                   allowablePars=NULL,incStates=NULL,
                                   mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                   useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                   peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                   microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                   stateWeights=NULL,extra_unfixed=extra_unfixed,
                                   do_call=TRUE,timeout=0)


jobs_colombia <- list(jobsA,jobsB,jobsC, jobsE, jobsF,jobsG,jobsH,jobsI)
