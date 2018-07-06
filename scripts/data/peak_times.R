devtools::load_all("~/Documents/Zika/zikaInfer")
library(ggplot2)

## Read in Colombia data
colombiaMicro <- read.csv("~/Documents/Zika/Data/colombia/microcephaly_dat_weekly_2017.csv",stringsAsFactors = FALSE)
colombiaMicro <- colombiaMicro[,c("startDay","endDay","microCeph","births","local")]
incCol <- read.csv("~/Documents/Zika/Data/colombia/zikv_inc_2017.csv",stringsAsFactors=FALSE)
incCol_confirmed <- read.csv("~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv",stringsAsFactors=FALSE)

## Read in report data
otherInc <- read.csv("~/Documents/Zika/Data/brazil/zika_inc_reports.csv",stringsAsFactors = FALSE)
stateDat <- read.csv("~/Documents/Zika/Data/brazil/microceph_reports_2016.csv",stringsAsFactors = FALSE)
stateDat <- stateDat[,c("startDay","endDay","microCeph","births","local")]
stateDat_confirmed <- read.csv("~/Documents/Zika/Data/brazil/microceph_reports_2016_confirmed.csv",stringsAsFactors = FALSE)
stateDat_confirmed <- stateDat_confirmed[,c("startDay","endDay","microCeph","births","local")]
#stateDat_confirmed$local <- paste0(stateDat_confirmed$local, " Confirmed")

## Read in data for Northeast Brazil (ie. Lancet data)
#ne_micro <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv",stringsAsFactors = FALSE)
#ne_zikv <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv",stringsAsFactors=FALSE)
## This is ZIKV incidence in pregnant women
ne_micro <- read.csv("~/Documents/Zika/Data/brazil/deOliveira2017/microCeph_deOliveira2017_clean.csv",stringsAsFactors = FALSE)
ne_zikv <- read.csv("~/Documents/Zika/Data/brazil/deOliveira2017/zikv_inc_deOliveira2017_clean.csv",stringsAsFactors=FALSE)
ne_micro$local <- "northeast"
ne_micro <- ne_micro[,c("startDay","endDay","microCeph","births","local")]
ne_zikv$local <- "northeast"
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

## Read in data for Salvador, Brazil
salv_micro <- read.csv("~/Documents/Zika/Data/brazil/salvador_microcephaly.csv",stringsAsFactors = FALSE)
salv_zikv <- read.csv("~/Documents/Zika/Data/brazil/salvador_inc_scaled.csv",stringsAsFactors = FALSE)
salv_micro$local <- "salvador"
salv_zikv$local <- "salvador"
salv_zikv <- salv_zikv[,c("startDay","endDay","inc","N_H","local")]
salv_micro <- salv_micro[,c("startDay","endDay","microCeph","births","local")]

otherDat <- rbind(colombiaMicro,stateDat, salv_micro, ne_micro)
newNames <- c(northeast="Northeast Brazil",colombia="Colombia",bahia="Bahia - reports", 
              pernambuco= "Pernambuco - reports", riograndedonorte="Rio Grande do Norte - reports",
              salvador="Salvador, Brazil")
otherDat$local <- newNames[otherDat$local]
otherDat$meanDay <- (otherDat$startDay + otherDat$endDay)/2
otherDat$local <- factor(otherDat$local, levels=c("Northeast Brazil","Colombia","Pernambuco - reports","Bahia - reports",
                                                  "Rio Grande do Norte - reports","Salvador, Brazil"))
otherDat$confirmed <- "Notified microcephaly"

incCol$local <- "colombia"
incCol <- incCol[,c("startDay","endDay","inc","N_H","local")]
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

otherInc <- otherInc[otherInc$local %in% c("bahia","riograndedonorte","pernambuco"),]
otherInc <- otherInc[,c("startDay","endDay","inc","N_H","local")]
otherInc <- rbind(otherInc,incCol, ne_zikv,salv_zikv)
otherInc$local <- newNames[otherInc$local]
otherInc$confirmed <- "Notified ZIKV"


# Report data confirmed ----------------------------------------------------
otherDat_confirmed <- rbind(stateDat_confirmed,ne_micro)
newNames <- c(northeast="Northeast Brazil",colombia="Colombia",bahia="Bahia - reports", pernambuco= "Pernambuco - reports", riograndedonorte="Rio Grande do Norte - reports")
otherDat_confirmed$local <- newNames[otherDat_confirmed$local]
otherDat_confirmed$meanDay <- (otherDat_confirmed$startDay + otherDat_confirmed$endDay)/2
otherDat_confirmed$local <- factor(otherDat_confirmed$local, levels=c("Northeast Brazil","Colombia","Pernambuco - reports","Bahia - reports","Rio Grande do Norte - reports"))
otherDat_confirmed$confirmed <- "Confirmed microcephaly"
otherInc_confirmed <- incCol_confirmed
otherInc_confirmed$local <- newNames[otherInc_confirmed$local]
otherInc_confirmed <- otherInc_confirmed[,c("startDay","endDay","inc","N_H","local")]
otherInc_confirmed$confirmed <- "Confirmed ZIKV"


otherInc <- rbind(otherInc, otherInc_confirmed)
otherDat <- rbind(otherDat, otherDat_confirmed)
otherDat <- otherDat[complete.cases(otherDat),]
otherInc <- otherInc[complete.cases(otherInc),]
otherInc$confirmed <- as.factor(otherInc$confirmed)
otherDat$confirmed <- as.factor(otherDat$confirmed)

# Peak times --------------------------------------------------------------
###########################################
## 3. PEAK TIME SUMMARIES
###########################################
# Get peak times for microcephaly data
peakTimes <- NULL
locals <- NULL
confirmed_vec <- NULL
peakMicros <- NULL
otherDat$local <- as.character(otherDat$local)
for(local in sort(unique(otherDat$local))){
  print(local)
  tmp <- otherDat[otherDat$local==local,c("startDay","endDay","microCeph","confirmed","births")]
  #if(local != "Colombia") tmp <- tmp[tmp$endDay < 1095,]
  for(confirmed in unique(tmp$confirmed)){
    print(confirmed)
    tmp1 <- tmp[tmp$confirmed == confirmed,]
    times <- rowMeans(tmp1[,c("startDay","endDay")])
    peakRow <- which.max(tmp1[,"microCeph"])
    #print(tmp1[peakRow,])
    peakMicro <- tmp1[peakRow,"microCeph"]/tmp1[peakRow,"births"]
    #peakTime <- mean(as.numeric(times[peakRow,c("startDay","endDay")]))
    peakTime <- times[peakRow]
    peakTimes <- c(peakTimes,peakTime)
    peakMicros <- c(peakMicros, peakMicro)
    locals <- c(locals, local)
    confirmed_vec <- c(confirmed_vec,confirmed)
  }
}
peakTimes_micro <- data.frame(peakTimes,local=locals,confirmed_vec,peakMicros)
#####################
# Get peak times for incidence
peakTimes_inc <- NULL
locals_inc <- NULL
confirmed_vec_inc <- NULL
peakIncs <- NULL
for(local in sort(unique(otherInc$local))){
  tmp <- otherInc[otherInc$local==local,c("startDay","endDay","inc","confirmed","N_H")]
  if(local != "Colombia") tmp <- tmp[tmp$endDay < 950,]
  for(confirmed in unique(tmp$confirmed)){
    tmp1 <- tmp[tmp$confirmed == confirmed,]
    times <- rowMeans(tmp1[,c("startDay","endDay")])
    peakRow <- which.max(tmp1[,"inc"])
    print(tmp1[peakRow,])
    peakInc <- tmp1[peakRow,"inc"]/tmp1[peakRow,"N_H"]
    #peakTime <- mean(as.numeric(times[peakRow,c("startDay","endDay")]))
    peakTime <- times[peakRow]
    peakTimes_inc <- c(peakTimes_inc,peakTime)
    locals_inc <- c(locals_inc, local)
    peakIncs <- c(peakIncs, peakInc)
    confirmed_vec_inc <- c(confirmed_vec_inc,confirmed)
  }
}
peakTimes_inc <- data.frame(peakTimes_inc,local=locals_inc,confirmed_vec_inc,peakIncs)

all_peaks <- merge(peakTimes_micro,peakTimes_inc, id.vars="local",all=TRUE)
all_peaks[all_peaks$local=="Pernambuco - reports","peakTimes_inc"] <- 804
all_peaks$diff <- all_peaks$peakTimes - all_peaks$peakTimes_inc
all_peaks$mean <- (all_peaks$peakTimes+all_peaks$peakTimes_inc)/2
all_peaks_numb <- all_peaks[c(1,2,4,6,8,10),]
all_peaks$peakTimes <- convert_number_to_date(all_peaks$peakTimes)
all_peaks$peakTimes_inc <- convert_number_to_date(all_peaks$peakTimes_inc)
all_peaks$mean <- convert_number_to_date(all_peaks$mean)
all_peaks$local <- as.factor(all_peaks$local)
all_peaks$peakIncs <- all_peaks$peakIncs*10000
all_peaks$peakMicros <- all_peaks$peakMicros*10000

ggplot(all_peaks) + geom_point(aes(y=peakTimes,x=local,group=confirmed_vec),col="blue",position=position_dodge(width=0.5)) + 
  geom_point(aes(y=peakTimes_inc,x=local,group=confirmed_vec_inc),col="red",position=position_dodge(width=0.5))  +
  geom_segment(aes(x=local,xend=local,y=peakTimes_inc,yend=peakTimes))+
  geom_text(aes(x=local,y=mean,label=paste0(round(diff/7,1)," weeks")),vjust=-1.5) +
  coord_flip()

