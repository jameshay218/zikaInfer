datFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
incFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
local <- "bahia"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]
parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365

combos <- expand.grid(runName = "northeast_forecast",chainNo=1:maxChain,stringsAsFactors=FALSE)

jobs_nejm_forecasts <- queuer::enqueue_bulk(obj1,combos,"model_forecasting",
                                            microDat=microDat,
                                            incDat=incDat,
                                            parTab=parTab,
                                            mcmcPars1=mcmcPars1,
                                            mcmcPars2=mcmcPars2,
                                            do_call=TRUE,timeout=0)
