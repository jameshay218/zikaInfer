parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv",stringsAsFactors=FALSE)
datFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
incFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
local <- "northeast"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]
parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365
parTab[parTab$names == "c","upper_bounds"] <- 5
parTab[parTab$names == "propn","upper_bounds"] <- 5
parTab[parTab$names == "zikv_reporting_change","values"] <- 0
maxChain <- 15
combos <- expand.grid(runName = "northeast_forecast_noreportingchange",chainNo=1:maxChain,stringsAsFactors=FALSE)

jobs_nejm_forecastsB <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)

parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv",stringsAsFactors=FALSE)
datFile <- "~/Documents/Zika/Data/brazil/deOliveira2017/microCeph_deOliveira2017_clean.csv"
incFile = "~/Documents/Zika/Data/brazil/deOliveira2017/zikv_inc_deOliveira2017_clean.csv"
local <- "northeast"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]
parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365
parTab[parTab$names == "c","upper_bounds"] <- 5
parTab[parTab$names == "propn","upper_bounds"] <- 5
parTab[parTab$names == "zikv_reporting_change","values"] <- 0
maxChain <- 15
combos <- expand.grid(runName = "northeast_forecast_noreportingchange_deOliveira2017",chainNo=1:maxChain,stringsAsFactors=FALSE)

jobs_nejm_forecastsB <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)
