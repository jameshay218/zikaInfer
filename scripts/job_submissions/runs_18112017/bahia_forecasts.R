parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv",stringsAsFactors=FALSE)
maxChain <- 5
datFile <- "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile <- "~/Documents/Zika/Data/brazil/zika_inc_reports.csv"
local <- "bahia"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]
parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365

## Assuming that all forecasting parameters are free
combos <- expand.grid(runName = "bahia_forecast",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_forecasts <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                             microDat=microDat,
                                             incDat=incDat,
                                             parTab=parTab,
                                             mcmcPars1=mcmcPars1,
                                             mcmcPars2=mcmcPars2,
                                             do_call=TRUE,timeout=0)

## Assuming that only incidence reporting could have changed
parTab_incPropn <- parTab
parTab_incPropn[parTab_incPropn$names %in% c("propn","abortion_rate","birth_reduction"),"fixed"] <- 1
parTab_incPropn[parTab_incPropn$names %in% c("abortion_rate","birth_reduction"),"values"] <- 0
combos <- expand.grid(runName = "bahia_incPropn",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_incPropn <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)

## Assuming only microcephaly reporting could have changed
parTab_propn <- parTab
parTab_propn[parTab_propn$names %in% c("incPropn2","abortion_rate","birth_reduction"),"fixed"] <- 1
parTab_propn[parTab_propn$names == "zikv_reporting_change","values"] <- 0
combos <- expand.grid(runName = "bahia_propn",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_propn <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                         microDat=microDat,
                                         incDat=incDat,
                                         parTab=parTab_propn,
                                         mcmcPars1=mcmcPars1,
                                         mcmcPars2=mcmcPars2,
                                         do_call=TRUE,timeout=0)

## Assuming only abortion rate could have changed
parTab_abortion <- parTab
parTab_abortion[parTab_abortion$names %in% c("incPropn2","propn","birth_reduction"),"fixed"] <- 1
parTab_abortion[parTab_abortion$names == "zikv_reporting_change","values"] <- 0
combos <- expand.grid(runName = "bahia_abortion",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_abortion <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab_abortion,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)

## Assuming only birth reductions could have changed
parTab_reduction <- parTab
parTab_reduction[parTab_reduction$names %in% c("incPropn2","abortion_rate","propn"),"fixed"] <- 1
parTab_reduction[parTab_reduction$names == "zikv_reporting_change","values"] <- 0
combos <- expand.grid(runName = "bahia_reduction",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_reduction <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                             microDat=microDat,
                                             incDat=incDat,
                                             parTab=parTab_reduction,
                                             mcmcPars1=mcmcPars1,
                                             mcmcPars2=mcmcPars2,
                                             do_call=TRUE,timeout=0)

jobs_bahia_forecast <- list(jobs_bahia_forecasts, jobs_bahia_incPropn,
                            jobs_bahia_propn,jobs_bahia_abortion,
                            jobs_bahia_reduction)
