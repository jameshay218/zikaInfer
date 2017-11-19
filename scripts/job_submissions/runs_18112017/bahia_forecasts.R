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
combos <- expand.grid(runName = "bahia_forecast",chainNo=1:maxChain,stringsAsFactors=FALSE)
jobs_bahia_forecasts <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                             microDat=microDat,
                                             incDat=incDat,
                                             parTab=parTab,
                                             mcmcPars1=mcmcPars1,
                                             mcmcPars2=mcmcPars2,
                                             do_call=TRUE,timeout=0)
