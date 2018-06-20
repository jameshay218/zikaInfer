datFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
incFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
local <- "northeast"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv",stringsAsFactors=FALSE)
maxChain <- 5

incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]

incDat$local <- "bahia"
microDat$local <- "bahia"

parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365
parTab[parTab$names == "N_H","values"] <- as.integer(parTab[parTab$names == "N_H","values"])


combos <- expand.grid(runName = "northeast_forecast",chainNo=1:maxChain,stringsAsFactors=FALSE)

jobs_nejm_forecasts <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)


datFile <- "~/Documents/Zika/Data/brazil/deOliveira2017/microCeph_deOliveira2017_clean.csv"
incFile = "~/Documents/Zika/Data/brazil/deOliveira2017/zikv_inc_deOliveira2017_clean.csv"
local <- "northeast"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)

microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]


incDat$local <- "bahia"
microDat$local <- "bahia"

parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365

combos <- expand.grid(runName = "northeast_forecast_deOliveira",chainNo=1:maxChain,stringsAsFactors=FALSE)

jobs_deOliveira_forecasts <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                                  microDat=microDat,
                                                  incDat=incDat,
                                                  parTab=parTab,
                                                  mcmcPars1=mcmcPars1,
                                                  mcmcPars2=mcmcPars2,
                                                  do_call=TRUE,timeout=0)
