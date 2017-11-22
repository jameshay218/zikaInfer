##############
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This script is used to create a table of allowable
## mosquito density and epidemic seed time parameters
library(zikaProj)
## Create allowable parameters for all states
data(locationInfo)
states <- unique(locationInfo$rawName)
parTab <- parTabSetup(states, 1, locationInfo,TRUE, NULL, TRUE, FALSE, FALSE, FALSE)


peakTimesB <- data.frame(local=states, start=927 - 60, end=927 + 60)

## Generate peak times
peakTimes <- rep(858,length(states))
peakTimes[which(states == "pernambuco")] <- 804
peakTimes[which(states == "bahia")] <- 855
peakTimes[which(states == "riograndedonorte")] <- 862

peakWidths <- rep(120,length(states))
peakWidths[which(states == "pernambuco")] <- 60
peakWidths[which(states == "bahia")] <- 60
peakWidths[which(states == "riograndedonorte")] <- 60
peakTimesA<- data.frame(local=states,start=peakTimes-peakWidths/2,end=peakTimes+peakWidths/2)



allowablePars_927_60 <- generateAllowableParams(peakTime=NULL, peakTimeRange=NULL,
                                                states,parTab,"",6,peakTimings=peakTimesB)
allowablePars_927_120 <- generateAllowableParams(peakTime=NULL, peakTimeRange=NULL,
                                                 states,parTab,"",6,peakTimings=peakTimesB)
allowablePars_858_60 <- generateAllowableParams(peakTime=NULL, peakTimeRange=NULL,
                                                states,parTab,"",6,peakTimings=peakTimesA)
allowablePars_858_120 <- generateAllowableParams(peakTime=NULL, peakTimeRange=NULL,
                                                 states,parTab,"",6,peakTimings=peakTimesA)

write.csv(allowablePars_927_60,"~/net/home/zika/data/real_data_60/allowablePars_927.csv",row.names=FALSE)
write.csv(allowablePars_927_120,"~/net/home/zika/data/real_data_120/allowablePars_927.csv",row.names=FALSE)
write.csv(allowablePars_858_60,"~/net/home/zika/data/real_data_60/allowablePars_858.csv",row.names=FALSE)
write.csv(allowablePars_858_120,"~/net/home/zika/data/real_data_120/allowablePars_858.csv",row.names=FALSE)
