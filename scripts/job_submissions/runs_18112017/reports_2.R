two_inc_states <- c("bahia","riograndedonorte")
two_week_states <- c("bahia","pernambuco")

datFile <- "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile <- "~/Documents/Zika/Data/brazil/brazil_report_zikv_inc_2016.csv"
microDat <- read.csv(datFile,stringsAsFactors=FALSE)
incDat <- read.csv(incFile,stringsAsFactors=FALSE)



peakTimes <- data.frame(local=c("bahia","pernambuco","riograndedonorte"),start=c(855,804,862)-30, end = c(855,804,862)+30)

maxChain <- 5
versions <- c(1,2,3)
priors <- NULL
extra_unfixed <- NULL
combos <- expand.grid(runName = "reports_2_weeks", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsA <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                     stateNames=two_week_states,microDat=microDat,
                                     incDat=NULL,preParTab=NULL,
                                     allowablePars=NULL,incStates=NULL,
                                     mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                     useInc=FALSE,usePeakTimes=TRUE,covMat=NULL,
                                     peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                     microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                     stateWeights=NULL,extra_unfixed=extra_unfixed,do_call=TRUE,timeout=0)


combos <- expand.grid(runName = "reports_2_inc", chainNo=1:maxChain,version=versions,stringsAsFactors=FALSE)
jobsB <- queuer::enqueue_bulk(obj1, combos, "model_fitting_script",
                                         stateNames=two_inc_states,microDat=microDat,incDat=incDat,preParTab=NULL,
                                         allowablePars=NULL,incStates=NULL,,
                                         mcmcPars=mcmcPars,mcmcPars2=mcmcPars2,allPriors=priors,
                                         useInc=TRUE,usePeakTimes=FALSE,covMat=NULL,
                                         peakTimeRange=NULL,peakTime=NULL,sim=FALSE, predict=FALSE,
                                         microChain=NULL,prePeakTimes=peakTimes,normLik=FALSE,
                                         stateWeights=NULL,extra_unfixed=extra_unfixed,do_call=TRUE,timeout=0)

jobs_reports2 <- list(jobsA,jobsB)
