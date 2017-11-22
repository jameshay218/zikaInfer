##############
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## Reads in all of the microcephaly and incidence data and produces
## incidence plots.
## The "Reading in data" section must be modified to specify full file paths to all of the
## incidence data. The save locations of the plots are specified throughout the script in sections
## labelled "Save plots"
library(zikaProj)
library(cowplot)
library(zoo)
library(ggplot2)
library(extrafont)

# Reading in data ---------------------------------------------------------
## Read in Brazil data
brazilMicro <- read.csv("~/Documents/Zika/Data/brazil/brazil_microCeph_28022017.csv",stringsAsFactors=FALSE)
brazilIncFaria <- read.csv("~/Documents/Zika/Data/brazil/zikv_inc_faria2016.csv",stringsAsFactors=FALSE)
brazilInc2017 <- otherInc <- read.csv("~/Documents/Zika/Data/brazil/zika_inc_reports.csv",stringsAsFactors = FALSE)

## Read in Colombia data
colombiaMicro <- read.csv("~/Documents/Zika/Data/colombia/microcephaly_dat_weekly.csv",stringsAsFactors = FALSE)
colombiaMicro <- colombiaMicro[,c("startDay","endDay","microCeph","births","local")]
incCol <- read.csv("~/Documents/Zika/Data/colombia/zikv_inc_2016.csv",stringsAsFactors=FALSE)

## Read in report data
stateDat <- read.csv("~/Documents/Zika/Data/brazil/microceph_reports_2016.csv",stringsAsFactors = FALSE)
stateDat <- stateDat[,c("startDay","endDay","microCeph","births","local")]

## Read in data for Northeast Brazil (ie. NEJM data)
ne_micro <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv",stringsAsFactors = FALSE)
ne_zikv <- read.csv("~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv",stringsAsFactors=FALSE)
ne_micro$local <- "northeast"
ne_micro <- ne_micro[,c("startDay","endDay","microCeph","births","local")]
ne_zikv$local <- "northeast"
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

use_states <- all_states <- zikaProj::get_epidemic_locations(brazilMicro,TRUE)$include

## Combine all data into a list, generate mean report days, and get full location names
allDats <- list("Brazil"=brazilMicro,"Brazil Inc"=brazilIncFaria,"Brazil 2017"=brazilInc2017,
                "Colombia"=colombiaMicro,"Reports"=stateDat,
                "Northeast micro"=ne_micro,"Northeast inc"=ne_zikv)

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

order <- locationInfo$fullName[sort(locationInfo$incidenceOrder,index.return=TRUE)$ix]
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
otherDat <- rbind(colombiaMicro,stateDat, ne_micro)
newNames <- c(northeast="Northeast Brazil NEJM",colombia="Colombia",bahia="Bahia - reports", pernambuco= "Pernambuco - reports", riograndedonorte="Rio Grande do Norte - reports")
otherDat$local <- newNames[otherDat$local]
otherDat$meanDay <- (otherDat$startDay + otherDat$endDay)/2
otherDat$local <- factor(otherDat$local, levels=c("Northeast Brazil NEJM","Colombia","Pernambuco - reports","Bahia - reports","Rio Grande do Norte - reports"))

incCol$local <- "colombia"
incCol <- incCol[,c("startDay","endDay","inc","N_H","local")]
ne_zikv <- ne_zikv[,c("startDay","endDay","inc","N_H","local")]

otherInc <- otherInc[otherInc$local %in% c("bahia","riograndedonorte","pernambuco"),]
otherInc <- otherInc[,c("startDay","endDay","inc","N_H","local")]
otherInc <- rbind(otherInc,incCol, ne_zikv)
otherInc$local <- newNames[otherInc$local]

peakTimes2 <- data.frame(local=c("Northeast Brazil NEJM","Colombia","Pernambuco - reports","Bahia - reports",
                                 "Rio Grande do Norte - reports"),
                         start=c(NA,NA, 804-30,NA,NA),
                         end=c(NA,NA,804+30,NA,NA))
actual2 <- data.frame(local=c("Northeast Brazil NEJM","Colombia","Pernambuco - reports","Bahia - reports",
                              "Rio Grande do Norte - reports"),peakTime=c(NA,NA,804,NA,NA))

# NE Brazil ----------------------------------------------------------------
local <- "Northeast Brazil NEJM"
scale2 <- 0.01
p2_ne <- ggplot() + 
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*70/N_H),stat="identity",alpha=0.4,fill="red",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/70))) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# Colombia ----------------------------------------------------------------
local <- "Colombia"
scale2 <- 0.01
p2_c <- ggplot() + 
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*70/N_H),stat="identity",alpha=0.4,fill="red",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/70))) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_b --------------------------------------------------------------------
local <- "Bahia - reports"
scale <- 200
scale2 <- 0.1
p2_b <- ggplot() + 
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*200/N_H),stat="identity",alpha=0.4,fill="red",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/200))) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_p --------------------------------------------------------------------
local <- "Pernambuco - reports"
scale2 <- 0.1
p2_p <- ggplot() +
  geom_rect(data=peakTimes2[peakTimes2$local==local,],aes(xmin=start,xmax=end,ymin=0,ymax=Inf,group=local),alpha=0.5,fill="red") +
  geom_vline(data=actual2[actual2$local==local,],aes(xintercept = peakTime,group=local),col="black",lty="dashed") +
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*1/N_H),stat="identity",alpha=0.4,fill="red",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2)) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)


# p2_r --------------------------------------------------------------------
local <- "Rio Grande do Norte - reports"
scale2 <- 0.1
p2_r <- ggplot() + 
  geom_area(data=otherInc[otherInc$local == local,],
            aes(x=startDay,y=inc*200/N_H),stat="identity",alpha=0.4,fill="red",col="black") +
  geom_area(data=otherDat[otherDat$local == local,],
            aes(x=meanDay,y=microCeph/births),stat="identity",alpha=0.5,fill="blue",col="black") + 
  facet_wrap(~local,scales="free_y",ncol=1) + theme_bw() + 
  scale_y_continuous(expand=c(0,0),limits=c(0,scale2),sec.axis=sec_axis(~.*(1/200))) +
  theme(axis.text.y=element_text(size=8,family="Arial"),
        axis.title.x=element_blank(),
        axis.title=element_text(size=10,family="Arial"),
        strip.text=element_text(size=8,family="Arial"),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=90,hjust=0.5,size=8,family="Arial"), 
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits=c(365,max(labels)),breaks=labels,labels=labels_names)

textDat <- data.frame(x=labels,labels=labels_names)
text <- ggplot(textDat) + geom_label(aes(x=x,label=))

p2 <- plot_grid(p2_ne,p2_c,p2_b,p2_p,p2_r,ncol=1,align="v",rel_heights=c(1,1,1,1,1.4))
p3 <- plot_grid(p1,p2,rel_widths=c(2,1))


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

