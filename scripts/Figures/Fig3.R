##############
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## Reads in all of the microcephaly and incidence data and produces
## incidence plots.
## The "Reading in data" section must be modified to specify full file paths to all of the
## incidence data. The save locations of the plots are specified throughout the script in sections
## labelled "Save plots"
## The end of this script also generates Fig3 of the main text
devtools::load_all("~/Documents/Zika/zikaInfer")
#library(zikaInfer)
library(cowplot)
library(ggplot2)
library(extrafont)
library(zoo)
setwd("~/Documents/Zika/zikaInfer")

use_locations <- c("Colombia","Pernambuco","Bahia","Rio Grande do Norte","Salvador, Brazil","Colombia","Northeast Brazil")

# Reading in data ---------------------------------------------------------
## Read in Brazil data
brazilMicro <- read.csv("RawData/Brazil/brazil_microCeph_28022017.csv",stringsAsFactors=FALSE)
brazilIncFaria <- read.csv("RawData/Brazil/zikv_inc_faria2016.csv",stringsAsFactors=FALSE)
brazilInc2017 <- otherInc <- read.csv("RawData/Brazil/zika_inc_reports.csv",stringsAsFactors = FALSE)

## Read in Colombia data
colombiaMicro <- read.csv("RawData/Colombia/microcephaly_dat_weekly_Dec_2017.csv",stringsAsFactors = FALSE)
colombiaMicro <- colombiaMicro[,c("startDay","endDay","microCeph","births","local")]
incCol <- read.csv("RawData/Colombia/zikv_inc_Dec_2017.csv",stringsAsFactors=FALSE)
incCol_confirmed <- read.csv("RawData/Colombia/zikv_inc_Dec_2017_confirmed.csv",stringsAsFactors=FALSE)


## Read in report data
stateDat <- read.csv("RawData/Brazil/microceph_reports_2016.csv",stringsAsFactors = FALSE)
stateDat <- stateDat[,c("startDay","endDay","microCeph","births","local")]
stateDat_confirmed <- read.csv("RawData/Brazil/microceph_reports_2016_confirmed.csv",stringsAsFactors = FALSE)
stateDat_confirmed <- stateDat_confirmed[,c("startDay","endDay","microCeph","births","local")]
#stateDat_confirmed$local <- paste0(stateDat_confirmed$local, " Confirmed")

## Read in data for Northeast Brazil (ie. Lancet data)
#ne_micro <- read.csv("RawData//brazil/Northeast/northeast_microceph.csv",stringsAsFactors = FALSE)
#ne_zikv <- read.csv("RawData//brazil/Northeast/northeast_zikv.csv",stringsAsFactors=FALSE)
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


use_states <- all_states <- zikaInfer::get_epidemic_locations(brazilMicro,TRUE)$include

## Combine all data into a list, generate mean report days, and get full location names
allDats <- list("Brazil"=brazilMicro,"Brazil Inc"=brazilIncFaria,"Brazil 2017"=brazilInc2017,
                "Colombia"=colombiaMicro,"Reports"=stateDat,"Reports Confirmed"=stateDat_confirmed,
                "Northeast micro"=ne_micro,"Northeast inc"=ne_zikv,"Salvador"=salv_micro)

allDats <- lapply(allDats, function(x){
    x$meanDay <- (x$startDay + x$endDay)/2
    return(x)
    })

allDats <- lapply(allDats, function(x){
  x$local <-   locationInfo$fullName[match(x$local, locationInfo$rawName)]
  return(x)
})


# Generate epidemic peak times by state -----------------------------------
## Generate peak times for states with incidence data

## Store peak times for states with report data
peakTimes <- data.frame(local=locationInfo$fullName, rawName=locationInfo$rawName, peakTimeReport=NA,
                        startReports=858-60,endReports=858+60,stringsAsFactors=FALSE)
peakTimes[peakTimes$rawName %in% c("pernambuco","bahia","riograndedonorte"),"peakTimeReport"] <- c(804,855,862)
peakTimes[peakTimes$rawName %in% c("pernambuco","bahia","riograndedonorte"),"startReports"] <- c(804,855,862) - 60
peakTimes[peakTimes$rawName %in% c("pernambuco","bahia","riograndedonorte"),"endReports"] <- c(804,855,862) + 60

## Store peak times  
brazilIncFaria$meanDay <- (brazilIncFaria$startDay + brazilIncFaria$endDay)/2
fariaPeaks <- plyr::ddply(brazilIncFaria, c("local"), .fun=function(x) x[which.max(x$inc),"meanDay"])
colnames(fariaPeaks) = c("rawName","peakTimes")
fariaPeaks <- data.frame(fariaPeaks,startFaria=fariaPeaks[,2]-30,endFaria=fariaPeaks[,2]+30,stringsAsFactors=FALSE)
peakTimes <- merge(peakTimes, fariaPeaks, by="rawName")

brazilIncFaria$local <- locationInfo$fullName[match(brazilIncFaria$local, locationInfo$rawName)]
brazilIncFaria$meanDay <- (brazilIncFaria$startDay + brazilIncFaria$endDay)/2
brazilMicro$meanDay <- (brazilMicro$startDay + brazilMicro$endDay)/2
brazilMicro$local <- locationInfo$fullName[match(brazilMicro$local, locationInfo$rawName)]

order <- unique(locationInfo$fullName[sort(locationInfo$incidenceOrder,index.return=TRUE)$ix])
brazilMicro$local <- factor(brazilMicro$local, levels=order)
brazilIncFaria$local <- factor(brazilIncFaria$local, levels=order)
peakTimes$local <- factor(peakTimes$local, levels=order)

brazilInc2017$meanDay <- (brazilInc2017$startDay + brazilInc2017$endDay)/2
brazilInc2017$local <- locationInfo$fullName[match(brazilInc2017$local, locationInfo$rawName)]
brazilInc2017$local <- factor(brazilInc2017$local, levels=order)

# Plots -------------------------------------------------------------------
## Generate plot labels
labels <- rep(getDaysPerMonth(3),4)
labels <- c(0,cumsum(labels))
labels_names <- as.yearmon(as.Date(labels,origin="2013-01-01"))

## Monthly microcephaly data from Bruno Zoca data
p1 <- ggplot() + 
  geom_rect(data=peakTimes,aes(xmin=startReports,xmax=endReports,ymin=0,ymax=Inf,group=local),alpha=0.5,fill="red") +
  geom_rect(data=peakTimes,aes(xmin=startFaria,xmax=endFaria,ymin=0,ymax=Inf,group=local),alpha=0.5,fill="orange") +
  geom_area(data=brazilMicro,aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  geom_vline(data=peakTimes,aes(xintercept = peakTimeReport,group=local),col="black",lty="dashed") +
  facet_wrap(~local,scales="fixed",ncol=4) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,0.02)) + 
  theme(axis.text.x=element_text(angle=90,vjust=0.5,size=8,family="Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        panel.grid.minor = element_blank()) +
  ylab("Per birth microcephaly incidence") +
  #xlab("Report Date") +
  xlab("") +
  scale_x_continuous(breaks=labels,labels=labels_names)

labels_inc <- rep(getDaysPerMonth(6),4)
labels_inc <- c(0,cumsum(labels_inc))
labels_names_inc <- as.yearmon(as.Date(labels_inc,origin="2013-01-01"))

## Incidence data from Faria et al
p_inc <- ggplot() +
  geom_area(data=brazilIncFaria,aes(x=meanDay,y=inc/N_H),stat="identity",alpha=0.5,fill="red",col="black") +
  facet_wrap(~local,scales="free_y",ncol=4) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,size=8,family="Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        panel.grid.minor = element_blank()) +
  ylab("Per capita ZIKV incidence (2015)") +
  #xlab("Report Date") +
  xlab("") +
  scale_x_continuous(breaks=labels_inc,labels=labels_names_inc)


## Incidence data from 2017 sources
p_inc2017 <- ggplot() +
  geom_area(data=brazilInc2017,aes(x=meanDay,y=inc/N_H),stat="identity",alpha=0.5,fill="red",col="black") +
  facet_wrap(~local,scales="free_y",ncol=4) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,vjust=0.5,size=8,family="Arial"),
        axis.text.y=element_text(size=8,family="Arial"),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        panel.grid.minor = element_blank()) +
  ylab("Per capita ZIKV incidence (2017)") +
  #xlab("Report Date") +
  xlab("") +
  scale_x_continuous(breaks=labels_inc,labels=labels_names_inc)


# Report data plots -------------------------------------------------------
otherDat <- rbind(colombiaMicro,stateDat, salv_micro)
newNames <- c(northeast="Northeast Brazil",colombia="Colombia",bahia="Bahia - reports", 
              pernambuco= "Pernambuco - reports", riograndedonorte="Rio Grande do Norte - reports",
              salvador="Salvador, Brazil")
otherDat$local <- newNames[otherDat$local]
otherDat$meanDay <- (otherDat$startDay + otherDat$endDay)/2
otherDat$local <- factor(otherDat$local, levels=c("Northeast Brazil","Colombia","Pernambuco - reports","Bahia - reports","Rio Grande do Norte - reports","Salvador, Brazil"))
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


peakTimes2 <- data.frame(local=c("Northeast Brazil","Colombia","Pernambuco - reports","Bahia - reports",
                                 "Rio Grande do Norte - reports"),
                         start=c(NA,NA, 804-60,NA,NA),
                         end=c(NA,NA,804+60,NA,NA))
actual2 <- data.frame(local=c("Northeast Brazil","Colombia","Pernambuco - reports","Bahia - reports",
                              "Rio Grande do Norte - reports"),peakTime=c(NA,NA,804,NA,NA))

peakTimes2$local <- as.character(peakTimes2$local)
peakTimes2[peakTimes2$local == "Pernambuco - reports","local"] <- "Pernambuco, Brazil"
actual2$local <- as.character(actual2$local)
actual2[actual2$local == "Pernambuco - reports","local"] <- "Pernambuco, Brazil"

blank_labels <- data.frame("confirmed"=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))

otherInc$local <- as.character(otherInc$local)
otherInc[otherInc$local == "Bahia - reports","local"] <- "Bahia, Brazil"
otherInc[otherInc$local == "Rio Grande do Norte - reports","local"] <- "Rio Grande do Norte, Brazil"
otherInc[otherInc$local == "Pernambuco - reports","local"] <- "Pernambuco, Brazil"

otherDat$local <- as.character(otherDat$local)
otherDat[otherDat$local == "Bahia - reports","local"] <- "Bahia, Brazil"
otherDat[otherDat$local == "Rio Grande do Norte - reports","local"] <- "Rio Grande do Norte, Brazil"
otherDat[otherDat$local == "Pernambuco - reports","local"] <- "Pernambuco, Brazil"

all_peaks_numb <- read.csv("~/Documents/Zika/zikaInfer/RawData/peak_times_for_plot.csv",stringsAsFactors=FALSE)
all_peaks_numb[all_peaks_numb$local %in% c("Bahia","Pernambuco","Rio Grande do Norte"),"local"] <- paste0(all_peaks_numb[all_peaks_numb$local %in% c("Bahia","Pernambuco","Rio Grande do Norte"),"local"],", Brazil")


# NE Brazil ----------------------------------------------------------------
local <- "Northeast Brazil"
scale2 <- 0.006*10000
y_lab_pos <- 0.005*10000
p2_ne <- ggplot() + 
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),
                                                                    y=y_lab_pos,vjust=-0.5) +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*100*10000/N_H,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +   facet_wrap(~local,scales="free_y",ncol=1) + 
 
  scale_fill_manual(values=c("red","blue","red","blue"),
                    labels=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                    limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))+
  theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/100))) +
  scale_fill_manual(values=c("blue","red","darkblue","darkred")) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)

# Colombia ----------------------------------------------------------------
local <- "Colombia"
scale2 <- 0.012*10000
scale3 <- 70
y_lab_pos <- 0.0105*10000
p2_c <- ggplot() + 
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),y=y_lab_pos,vjust=-0.5) +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*scale3*10000/N_H,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +   facet_wrap(~local,scales="free_y",ncol=1) + 
  theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/scale3))) +
  scale_fill_manual(values=c("red","blue","red","blue"),
    labels=c("Confirmed ZIKV    ","Confirmed microcephaly    ","Notified ZIKV    ","Notified microcephaly"    ),
    limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                     labels=c("Confirmed ZIKV    ","Confirmed microcephaly    ","Notified ZIKV    ","Notified microcephaly"    ))+
  guides(fill=guide_legend(ncol=4))+
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        text=element_text(size=8,family="Arial"),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text=element_text(size=8,family="Arial"),
        legend.title=element_blank(),
        legend.position = "bottom") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)
legend <- cowplot::get_legend(p2_c)
p2_c <- p2_c + theme(legend.position="none")
plot(legend)
# p2_b --------------------------------------------------------------------
local <- "Bahia, Brazil"
scale <- 200
scale2 <- 0.1*10000
y_lab_pos <- 0.085*10000
p2_b <- ggplot() +
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),y=y_lab_pos,vjust=-0.5) +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*200*10000/N_H,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),stat="identity",position=position_stack(reverse=TRUE),col="black") +
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_fill_manual(values=c("red","blue","red","blue"),
                   labels=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                   limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/200))) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        text=element_text(size=8,family="Arial"),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_p --------------------------------------------------------------------
local <- "Pernambuco, Brazil"
scale2 <- 0.11*10000
y_lab_pos <- 0.09*10000
p2_p <- ggplot() +
  geom_rect(data=peakTimes2[peakTimes2$local==local,],aes(xmin=start,xmax=end,ymin=0,ymax=Inf,group=local),alpha=0.5,fill="red") +
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),y=y_lab_pos,vjust=-0.5) +
  geom_vline(data=actual2[actual2$local==local,],aes(xintercept = peakTime,group=local),col="black",lty="dashed") +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*1*10000/N_H,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") +   
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_fill_manual(values=c("red","blue","red","blue"),
                    labels=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                    limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))+
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2)) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        text=element_text(size=8,family="Arial"),
        axis.text.x=element_text(angle=90,hjust=0.5,size=8,family="Arial"), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_r --------------------------------------------------------------------
local <- "Rio Grande do Norte, Brazil"
scale2 <- 0.05*10000
y_lab_pos <- 0.041*10000
p2_r <- ggplot() + 
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),y=y_lab_pos,vjust=-0.5) +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*200*10000/N_H,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") +  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/200))) +
  scale_fill_manual(values=c("red","blue","red","blue"),
                    labels=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                    limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))+
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_s --------------------------------------------------------------------
local <- "Salvador, Brazil"
scale2 <- 0.08*10000
y_lab_pos <- 0.065*10000
p2_s <- ggplot() + 
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*0.5*10000/N_H,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph*10000/births,alpha=confirmed,fill=confirmed),position=position_stack(reverse=TRUE),stat="identity",col="black") + 
  geom_segment(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=peakTimes_inc,xend=peakTimes),y=y_lab_pos,yend=y_lab_pos)+
  geom_text(data=all_peaks_numb[all_peaks_numb$local == local,],aes(x=(peakTimes+peakTimes_inc)/2,label=paste0(signif(diff/7,3)," weeks")),y=y_lab_pos,vjust=-0.5) +
  
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(0.5))) +
  scale_fill_manual(values=c("red","blue","red","blue"),
                    labels=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"),
                    limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly")) +
  scale_alpha_manual(values=c(0.9,0.9,0.4,0.4),
                     limits=c("Confirmed ZIKV","Confirmed microcephaly","Notified ZIKV","Notified microcephaly"))+
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        text=element_text(size=8,family="Arial"),
        axis.text.x=element_text(angle=90,hjust=0.5,size=8,family="Arial"), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none") +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)

textDat <- data.frame(x=labels,labels=labels_names)
text <- ggplot(textDat) + geom_label(aes(x=x,label=))

p2 <- plot_grid(p2_ne,p2_c,p2_b,p2_p,p2_r,ncol=1,align="v",rel_heights=c(1,1,1,1,1.4))
p3 <- plot_grid(p1,p2,rel_widths=c(2,1))

p4 <- plot_grid(p2_ne,p2_c,p2_b,p2_r,p2_s,p2_p,ncol=2,align="v",rel_heights = c(1,1,1.3,1,1,1.3))
p4 <- plot_grid(p4, legend, ncol=1,rel_heights=c(1,0.1))

## Save plot
cairo_ps("Fig3.eps",width=7,height=7,family="Arial")
print(p4)
dev.off()

## Save plot
svg("Fig3.svg",width=7,height=7,family="Arial")
print(p4)
dev.off()


png("data_plot.png",width=7,height=7,family="Arial",units = "in",res=300)
print(p4)
dev.off()

# Save plots --------------------------------------------------------------
###################
## Save plot
cairo_ps("microtableau_data.eps",width=7,height=8,family="Arial")
print(p1)
dev.off()

png("microtableau.png",width=7,height=8,family="Arial",units = "in",res=300)
print(p1)
dev.off()
##################


###################
## Save plot inc
cairo_ps("faria2016_zikv_incidence_data.eps",width=7,height=8,family="Arial")
print(p_inc)
dev.off()

png("faria2016_zikv_incidence_data.png",width=7,height=8,family="Arial",units = "in",res=300)
print(p_inc)
dev.off()
##################

###################
## Save plot inc 2017
cairo_ps("report2017_zikv_incidence_data.eps",width=7,height=8,family="Arial")
print(p_inc2017)
dev.off()

png("report2017_zikv_incidence_data.png",width=7,height=8,family="Arial",units = "in",res=300)
print(p_inc2017)
dev.off()
##################

###################
## Save plots
png("all_data.png",width=10,height=8,units="in",res=300)
print(p3)
dev.off()

cairo_ps("all_data.eps",width=10,height=8,family="Arial")
print(p3)
dev.off()

cairo_ps("report_data.eps",width=3,height=6,family="Arial")
print(p2)
dev.off()
##################

