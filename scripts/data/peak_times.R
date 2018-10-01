##############
## JAMES HAY 09.07.2018 - jameshay218@gmail.com
## Reads in all of the microcephaly and incidence data and produces
## plots and tables showing lag between peak ZIKV and microcephaly
## if you source this file, look for "peakTimes_data" as the final output object

#devtools::load_all("~/Documents/Zika/zikaInfer")
library("zikaInfer")
library(ggplot2)
setwd("~/Documents/Zika/zikaInfer") ## Change to wherever you've saved the git repo

## Read in Colombia data
colombiaMicro <- read.csv("RawData/Colombia/microcephaly_dat_weekly_Dec_2017.csv",stringsAsFactors = FALSE)
colombiaMicro <- colombiaMicro[,c("startDay","endDay","microCeph","births","local")]
incCol <- read.csv("RawData/Colombia/zikv_inc_Dec_2017.csv",stringsAsFactors=FALSE)
incCol_confirmed <- read.csv("RawData/Colombia/zikv_inc_Dec_2017_confirmed.csv",stringsAsFactors=FALSE)

## Read in report data
all_incDat <- read.csv("RawData/Brazil/brazil_report_zikv_inc_2016.csv",stringsAsFactors = FALSE)
stateDat <- read.csv("RawData/Brazil/microceph_reports_2016.csv",stringsAsFactors = FALSE)
stateDat <- stateDat[,c("startDay","endDay","microCeph","births","local")]
stateDat_confirmed <- read.csv("RawData/Brazil/microceph_reports_2016_confirmed.csv",stringsAsFactors = FALSE)
stateDat_confirmed <- stateDat_confirmed[,c("startDay","endDay","microCeph","births","local")]
#stateDat_confirmed$local <- paste0(stateDat_confirmed$local, " Confirmed")

## Read in data for Northeast Brazil (ie. Lancet data)
#ne_micro <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv",stringsAsFactors = FALSE)
#ne_zikv <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv",stringsAsFactors=FALSE)
## This is ZIKV incidence in pregnant women
ne_micro <- read.csv("RawData/Northeast/deOliveira2017/microCeph_deOliveira2017_clean.csv",stringsAsFactors = FALSE)
ne_zikv <- read.csv("RawData/Northeast/deOliveira2017/zikv_inc_deOliveira2017_clean.csv",stringsAsFactors=FALSE)
ne_micro$local <- "northeast"
ne_micro <- ne_micro[,c("startDay","endDay","microCeph","births","local")]
ne_zikv$local <- "northeast"
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

## Read in data for Salvador, Brazil
salv_micro <- read.csv("RawData/Brazil/salvador_microcephaly.csv",stringsAsFactors = FALSE)
salv_zikv <- read.csv("RawData/Brazil/salvador_inc_scaled.csv",stringsAsFactors = FALSE)
salv_micro$local <- "salvador"
salv_zikv$local <- "salvador"
salv_zikv <- salv_zikv[,c("startDay","endDay","inc","N_H","local")]
salv_micro <- salv_micro[,c("startDay","endDay","microCeph","births","local")]

all_microDat <- rbind(colombiaMicro,stateDat, salv_micro, ne_micro)
newNames <- c(northeast="Northeast Brazil",colombia="Colombia",bahia="Bahia", 
              pernambuco= "Pernambuco", riograndedonorte="Rio Grande do Norte",
              salvador="Salvador, Brazil")
all_microDat$local <- newNames[all_microDat$local]
all_microDat$meanDay <- (all_microDat$startDay + all_microDat$endDay)/2
all_microDat$local <- factor(all_microDat$local, levels=c("Northeast Brazil","Colombia","Pernambuco","Bahia",
                                                  "Rio Grande do Norte","Salvador, Brazil"))
all_microDat$confirmed <- "Notified microcephaly"

incCol$local <- "colombia"
incCol <- incCol[,c("startDay","endDay","inc","N_H","local")]
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

all_incDat <- all_incDat[all_incDat$local %in% c("bahia","riograndedonorte","pernambuco"),]
all_incDat <- all_incDat[,c("startDay","endDay","inc","N_H","local")]
all_incDat <- rbind(all_incDat,incCol, ne_zikv,salv_zikv)
all_incDat$local <- newNames[all_incDat$local]
all_incDat$confirmed <- "Notified ZIKV"


# Report data confirmed ----------------------------------------------------
all_microDat_confirmed <- rbind(stateDat_confirmed,ne_micro)
newNames <- c(northeast="Northeast Brazil",colombia="Colombia",bahia="Bahia", pernambuco= "Pernambuco", riograndedonorte="Rio Grande do Norte")
all_microDat_confirmed$local <- newNames[all_microDat_confirmed$local]
all_microDat_confirmed$meanDay <- (all_microDat_confirmed$startDay + all_microDat_confirmed$endDay)/2
all_microDat_confirmed$local <- factor(all_microDat_confirmed$local, levels=c("Northeast Brazil","Colombia","Pernambuco","Bahia","Rio Grande do Norte"))
all_microDat_confirmed$confirmed <- "Confirmed microcephaly"
all_incDat_confirmed <- incCol_confirmed
all_incDat_confirmed$local <- newNames[all_incDat_confirmed$local]
all_incDat_confirmed <- all_incDat_confirmed[,c("startDay","endDay","inc","N_H","local")]
all_incDat_confirmed$confirmed <- "Confirmed ZIKV"


all_incDat <- rbind(all_incDat, all_incDat_confirmed)
all_microDat <- rbind(all_microDat, all_microDat_confirmed)
all_microDat <- all_microDat[complete.cases(all_microDat),]
all_incDat <- all_incDat[complete.cases(all_incDat),]
all_incDat$confirmed <- as.factor(all_incDat$confirmed)
all_microDat$confirmed <- as.factor(all_microDat$confirmed)

# Peak times --------------------------------------------------------------
# Get peak times for microcephaly data
peakTimes <- NULL
locals <- NULL
confirmed_vec <- NULL
peakMicros <- NULL
all_microDat$local <- as.character(all_microDat$local)
## For each location
for(local in sort(unique(all_microDat$local))){
  print(local)
  tmp <- all_microDat[all_microDat$local==local,c("startDay","endDay","microCeph","confirmed","births")]
  
  ## For confirmed and notified cases
  for(confirmed in unique(tmp$confirmed)){
    print(confirmed)
    tmp1 <- tmp[tmp$confirmed == confirmed,]
    
    ## Get mean day (ie. mid point of an epidemiological week)
    times <- rowMeans(tmp1[,c("startDay","endDay")])
    
    ## Find peak microcephaly
    peakRow <- which.max(tmp1[,"microCeph"])
    
    ## Get per birth microcephaly
    peakMicro <- tmp1[peakRow,"microCeph"]/tmp1[peakRow,"births"]
    
    ## Store these results
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

## For each location
for(local in sort(unique(all_incDat$local))){
  tmp <- all_incDat[all_incDat$local==local,c("startDay","endDay","inc","confirmed","N_H")]
  
  ## Make sure that we are only looking at first wave
  ## ie. first half of 2015 - convert_number_to_date(950) = "2015-08-09"
  if(local != "Colombia") tmp <- tmp[tmp$endDay < 950,]
  
  ## For confirmed and notified ZIKV cases
  for(confirmed in unique(tmp$confirmed)){
    tmp1 <- tmp[tmp$confirmed == confirmed,]
    ## Mean of start and end day of EW
    times <- rowMeans(tmp1[,c("startDay","endDay")])
    peakRow <- which.max(tmp1[,"inc"])
    print(tmp1[peakRow,])
    
    ## Get per cap incidence
    peakInc <- tmp1[peakRow,"inc"]/tmp1[peakRow,"N_H"]
    peakTime <- times[peakRow]
    peakTimes_inc <- c(peakTimes_inc,peakTime)
    locals_inc <- c(locals_inc, local)
    peakIncs <- c(peakIncs, peakInc)
    confirmed_vec_inc <- c(confirmed_vec_inc,confirmed)
  }
}
peakTimes_inc <- data.frame(peakTimes_inc,local=locals_inc,confirmed_vec_inc,peakIncs)

## 
all_peaks <- merge(peakTimes_micro,peakTimes_inc, id.vars="local",all=TRUE)

## Peak of ZIKV for Pernambuco is just based on the 3 weeks that reported ZIKV
## in 2015
all_peaks[all_peaks$local=="Pernambuco","peakTimes_inc"] <- 804

## Find lags
all_peaks$diff <- all_peaks$peakTimes - all_peaks$peakTimes_inc
all_peaks$mean <- (all_peaks$peakTimes+all_peaks$peakTimes_inc)/2

## Extract subset of these peaks for data fig 3
all_peaks_numb <- all_peaks[c(1,2,4,6,8,10),]

## Converts integers to Date objects, with 2015-01-01 = 0
all_peaks$peakTimes <- convert_number_to_date(all_peaks$peakTimes)
all_peaks$peakTimes_inc <- convert_number_to_date(all_peaks$peakTimes_inc)
all_peaks$mean <- convert_number_to_date(all_peaks$mean)
all_peaks$local <- as.factor(all_peaks$local)

## Get per 10,000 incidence
all_peaks$peakIncs <- all_peaks$peakIncs*10000
all_peaks$peakMicros <- all_peaks$peakMicros*10000

## Look at peak lags
ggplot(all_peaks) + geom_point(aes(y=peakTimes,x=local,group=confirmed_vec),col="blue",position=position_dodge(width=0.5)) + 
  geom_point(aes(y=peakTimes_inc,x=local,group=confirmed_vec_inc),col="red",position=position_dodge(width=0.5))  +
  geom_segment(aes(x=local,xend=local,y=peakTimes_inc,yend=peakTimes))+
  geom_text(aes(x=local,y=mean,label=paste0(round(diff/7,1)," weeks")),vjust=-1) +
  facet_wrap(~confirmed_vec,nrow=2) + 
  coord_flip()


# Bootstrap peak times ----------------------------------------------------
## Function to bootstrap incidence curves
## dat1 - vector of incidence per unit time for eg. ZIKV infection
## dat2 - vector of incidence per unit time for eg. microcephaly
## weeks1 - vector of reporting times corresponding to dat1 (eg. epidemiological weeks)
## weeks2 - vector of reporting times corresponding to dat2
## N - how many curves to bootstrap
bootstrap_sample <- function(dat1, dat2, weeks1, weeks2, N=1000){
  ## Scale incidence curve to give proportion
  ## of incidence observed at each time
  total_samps_A <- sum(dat1)
  scaled1 <- dat1/total_samps_A
  
  total_samps_B <- sum(dat2)
  scaled2 <- dat2/total_samps_B
  
  ## Empty vectors to store bootstrapped curves
  empty1 <- numeric(length(weeks1))
  names(empty1) <- weeks1
  empty2 <- numeric(length(weeks2))
  names(empty2) <- weeks2
  
  ## Store lags between peaks of bootstrapped incidence curves
  peaks <- numeric(N)
  for(i in 1:N){
    ## Sample times of infection A with replacement
    new1 <- sample(weeks1, size=total_samps_A, prob=scaled1,replace=TRUE)
    new1 <- table(new1)
    empty1[names(new1)] <- as.numeric(new1)
    
    ## Sample times of infection B with replacement
    new2 <- sample(weeks2, size=total_samps_B, prob=scaled2,replace=TRUE)
    new2 <- table(new2)
    empty2[names(new2)] <- as.numeric(new2)
    
    ## Get peaks of this bootstrapped curve
    peak1 <- weeks1[which.max(empty1)]
    peak2 <- weeks2[which.max(empty2)]
    
    ## Find difference between these
    dif <- peak1-peak2
    peaks[i] <- dif
    
  }
  return(peaks)
}

## Alternative function to above - assume poisson distributed observation process
## (as opposed to binomial above)
bootstrap_pois <- function(dat1, dat2, weeks1, weeks2,N=1000){
  peaks <- numeric(N)
  for(i in 1:N){
    new1 <- rpois(length(dat1),lambda=dat1)
    new2 <- rpois(length(dat2),lambda=dat2)
    
    peak1 <- weeks1[which.max(new1)]
    peak2 <- weeks2[which.max(new2)]
    dif <- peak1-peak2
    peaks[i] <- dif
  }
  return(peaks)
}

#bootstrapped <- NULL
all_incDat$local <- as.character(all_incDat$local)
all_incDat$confirmed <- as.character(all_incDat$confirmed)
all_microDat$local <- as.character(all_microDat$local)
all_microDat$confirmed <- as.character(all_microDat$confirmed)

peakTime_data <- NULL

for(local in unique(all_incDat$local)){
  print(local)
  tmpInc <- all_incDat[all_incDat$local == local,]
  tmpMicro <- all_microDat[all_microDat$local == local,]
  for(confirmed_zikv in unique(tmpInc$confirmed)){
    for(confirmed_micro in unique(tmpMicro$confirmed)){
      ## No incidence curve to sample for Pernambuco
      if(local != "Pernambuco"){
        dat1 <- tmpInc[tmpInc$confirmed == confirmed_zikv,]
        ## If not Colombia, need to only look at first wave
        if(local != "Colombia") dat1 <- dat1[dat1$endDay < 950,]
        inc <- dat1$inc
        weeks2 <- rowMeans(dat1[,c("startDay","endDay")])
        
        dat2 <- tmpMicro[tmpMicro$confirmed == confirmed_micro,]
        microCeph <- dat2$microCeph
        weeks1 <- rowMeans(dat2[,c("startDay","endDay")])
        
        #y <- bootstrap_pois(microCeph,inc,weeks1,weeks2,N=50)
        y1 <- bootstrap_sample(microCeph,inc,weeks1,weeks2,N=50)
        
        ## get 95% confidence intervals and look by week rather than day
        #bootstrapped_pois[[local]] <- c(mean(y/7),quantile(y/7,c(0.025,0.975)))
        peakTime_data <- rbind(peakTime_data,c(local, confirmed_zikv, confirmed_micro, mean(y1/7),quantile(y1/7,c(0.025,0.975))))
      }
    }
  }
}
colnames(peakTime_data) <- c("Location","ZIKV","Microcephaly","Mean","2.5%","97.5%")
peakTime_data <- as.data.frame(peakTime_data)
peakTime_data[,"Mean"] <- as.numeric(as.character(peakTime_data[,"Mean"]))
peakTime_data[,"2.5%"] <- as.numeric(as.character(peakTime_data[,"2.5%"]))
peakTime_data[,"97.5%"] <- as.numeric(as.character(peakTime_data[,"97.5%"]))

ggplot(peakTime_data) + 
  geom_pointrange(aes(x=Location,ymin=`2.5%`,ymax=`97.5%`,y=Mean)) + 
  facet_grid(ZIKV~Microcephaly) + 
  scale_y_continuous(limits=c(0,50)) + 
  xlab("Lag in weeks")+
  coord_flip()
