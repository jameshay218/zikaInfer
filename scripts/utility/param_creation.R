library(zikaProj)
## Create allowable parameters for all states
data(locationInfo)
states <- unique(locationInfo$rawName)
parTab <- parTabSetup(states, 1, locationInfo,TRUE, NULL, TRUE, FALSE, FALSE, FALSE)

allowablePars_927_60 <- generateAllowableParams(peakTime=927, peakTimeRange=60,states,parTab,"",6)
allowablePars_927_120 <- generateAllowableParams(peakTime=927, peakTimeRange=120,states,parTab,"",6)
allowablePars_858_60 <- generateAllowableParams(peakTime=858, peakTimeRange=60,states,parTab,"",6)
allowablePars_858_120 <- generateAllowableParams(peakTime=858, peakTimeRange=120,states,parTab,"",6)

write.csv(allowablePars_927_60,"~/net/home/zika/data/real_data_60/allowablePars_927.csv",row.names=FALSE)
write.csv(allowablePars_927_120,"~/net/home/zika/data/real_data_120/allowablePars_927.csv",row.names=FALSE)
write.csv(allowablePars_858_60,"~/net/home/zika/data/real_data_60/allowablePars_858.csv",row.names=FALSE)
write.csv(allowablePars_858_120,"~/net/home/zika/data/real_data_120/allowablePars_858.csv",row.names=FALSE)
