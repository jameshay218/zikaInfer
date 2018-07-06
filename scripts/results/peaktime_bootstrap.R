bootstrap_sample <- function(dat1, dat2, weeks1, weeks2, N=1000){
  total_samps_A <- sum(dat1)
  scaled1 <- dat1/total_samps_A
  
  total_samps_B <- sum(dat2)
  scaled2 <- dat2/total_samps_B
  empty1 <- numeric(length(weeks1))
  names(empty1) <- weeks1
  empty2 <- numeric(length(weeks2))
  names(empty2) <- weeks2
  
  peaks <- numeric(N)
  for(i in 1:N){
     new1 <- sample(weeks1, size=total_samps_A, prob=scaled1,replace=TRUE)
    new1 <- table(new1)
    empty1[names(new1)] <- as.numeric(new1)
    
    new2 <- sample(weeks2, size=total_samps_B, prob=scaled2,replace=TRUE)
    new2 <- table(new2)
    empty2[names(new2)] <- as.numeric(new2)
  
    peak1 <- weeks1[which.max(empty1)]
    peak2 <- weeks2[which.max(empty2)]
    
    dif <- peak1-peak2
    peaks[i] <- dif
    
  }
  return(peaks)
}


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
otherDat$confirmed <- "Suspected microcephaly"

incCol$local <- "colombia"
incCol <- incCol[,c("startDay","endDay","inc","N_H","local")]
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

otherInc <- otherInc[otherInc$local %in% c("bahia","riograndedonorte","pernambuco"),]
otherInc <- otherInc[,c("startDay","endDay","inc","N_H","local")]
otherInc <- rbind(otherInc,incCol, ne_zikv,salv_zikv)
otherInc$local <- newNames[otherInc$local]
otherInc$confirmed <- "Suspected ZIKV"


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

#bootstrapped <- NULL
bootstrapped_pois <- NULL
bootstrapped_norm <- NULL
for(local in unique(otherInc$local)){
  print(local)
  if(local != "Pernambuco - reports"){
    dat1 <- otherInc[otherInc$local == local,]
    if(local != "Colombia") dat1 <- dat1[dat1$endDay < 950,]
    inc <- dat1$inc
    weeks2 <- rowMeans(dat1[,c("startDay","endDay")])

    dat2 <- otherDat[otherDat$local == local,]
    microCeph <- dat2$microCeph
    weeks1 <- rowMeans(dat2[,c("startDay","endDay")])
  
    y <- bootstrap_pois(microCeph,inc,weeks1,weeks2,N=10000)
    y1 <- bootstrap_sample(microCeph,inc,weeks1,weeks2,N=1000)
    bootstrapped_pois[[local]] <- c(mean(y/7),quantile(y/7,c(0.025,0.975)))
    bootstrapped_norm[[local]] <- c(mean(y1/7),quantile(y1/7,c(0.025,0.975)))
  }
}

