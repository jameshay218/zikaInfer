###########################################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This is a script to generate statistics
## for the Zika manuscript. It generates statistics
## on the microcephaly and incidence data as follows:
## 1. Reading in data
##    a. Plot all data
## 2. Raw Pearson correlations between state data
##    a. Pearson correlation coefficient
##    b. Median and 95% quantiles
##    c. Boxplot
##    d. Cross-correlations
##    e. Spatio-temporal correlations
##    f. Correlation plotsother
## 3. Summary stats on peak times
## 4. Summary stats on max cases
## 5. Save a .csv file with all of the summary statistics (peak time and max cases)
###########################################
library(zikaInfer)
library(ggplot2)
library(dplyr)
library(plyr)
library(ncf)
library(Hmisc)
library(corrplot)
###########
## NOTE
## I must have moved the below source file at some point - I think all of the functions
## were moved into the main package
###########
#source("plotting_help_functions.R")

brazilIncDat <- read.csv("~/Documents/Zika/Data/brazil/zikv_inc_faria2016.csv",stringsAsFactors = FALSE)
brazilMicroDat <- read.csv("~/Documents/Zika/Data/brazil/brazil_microCeph_28022017.csv",stringsAsFactors = FALSE)

## Flags for script
SUBSET_STATES <- FALSE
EPIDEMIC_STATES <- FALSE
NORMALISE <- FALSE
SAVE_PLOT <- FALSE
LOG_CORRELATION <- FALSE

# Data --------------------------------------------------------------------
###########################################
## 1. DATA FORMATTING
###########################################
## Read in data
if(EPIDEMIC_STATES){
  print("Using only `epidemic` states")
} else {
  print("Using all states")
}

if(NORMALISE){
  print("Normalised")
}


## Get subset of data and then reformat data to give columns as states
brazilMicroDat$per_birth <- brazilMicroDat$microCeph/brazilMicroDat$births
per_birth <- ddply(brazilMicroDat,~local, function(x) max(x$per_birth))
per_birth$V1 <- per_birth$V1*10000
per_birth$V1 <- as.integer(per_birth$V1)
microDat <- brazilMicroDat[,c("startDay","microCeph","local","births")]
incDat <- brazilIncDat
microcephaly <- NULL
totalBirths <- NULL

## MICROCEPHALY
## For each state, get column of microcephaly and births
for(i in sort(unique(microDat$local))){
  tmp_births <- microDat[microDat$local==i,"births"]
  tmp_microceph <- microDat[microDat$local==i,"microCeph"]
  microcephaly <- cbind(microcephaly, tmp_microceph)
  totalBirths <- cbind(totalBirths, tmp_births)
}

## Correct column names for each state
colnames(microcephaly) <- colnames(totalBirths) <- sort(unique(microDat$local))

# Make data frame and get totals
microDat <- as.data.frame(microcephaly)
microDat <- cbind(microDat,"Total"=rowSums(microDat))
totalBirths <- as.data.frame(totalBirths)
totalBirths <- cbind(totalBirths,"Total"=rowSums(totalBirths))

## INCIDENCE
inc <- NULL
for(i in sort(unique(incDat$local))){
  tmp_inc <- incDat[incDat$local==i,"inc"]
  inc <- cbind(inc, tmp_inc)
}
colnames(inc) <- sort(unique(incDat$local))
formattedIncDat <- as.data.frame(inc)
formattedIncDat <- cbind(formattedIncDat,"Total"=rowSums(formattedIncDat))

## Subset data frame for cross correlations
totalMicroSubset <- brazilMicroDat[brazilMicroDat$local == "bahia",c("buckets","startDay","endDay","microCeph")]
totalMicroSubset$microCeph <- microDat$Total
totalIncSubset <- incDat[incDat$local=="bahia",c("buckets","startDay","endDay","inc")]
totalIncSubset$inc <- formattedIncDat$Total

## Convert to daily normalised inc
perDayI <- totalIncSubset$inc/totalIncSubset$buckets
perDayI <- rep(perDayI,totalIncSubset$buckets)
perDayI <- data.frame("time"=seq(min(totalIncSubset$startDay),max(totalIncSubset$endDay)-1,by=1),"inc"=perDayI)
perDayI$inc <- perDayI$inc/sum(perDayI$inc)
weekly_inc <- colSums(matrix(perDayI$inc,nrow=7))
weekly_inc <- data.frame("week"=seq(53,length(weekly_inc)+52,by=1),"inc"=weekly_inc)

perDayM <- totalMicroSubset$microCeph/totalMicroSubset$buckets
perDayM <- rep(perDayM,totalMicroSubset$buckets)
perDayM <- data.frame("time"=seq(min(totalMicroSubset$startDay),max(totalMicroSubset$endDay)-1,by=1),"micro"=perDayM)
perDayM$micro <- perDayM$micro/sum(perDayM$micro)
weekly_microceph <- colSums(matrix(perDayM$micro,nrow=7))
weekly_microceph <- data.frame("week"=seq(1,length(weekly_microceph),by=1),"microceph"=weekly_microceph)
weekly_microceph <- weekly_microceph[weekly_microceph$week >= min(weekly_inc$week),]
#weekly_microceph <- weekly_microceph[(nrow(weekly_microceph)-nrow(weekly_inc)):nrow(weekly_microceph),]

normalised_microDat <- apply(microDat[,colnames(microDat) != "Total"], 2, function(x) x/sum(x))
normalised_incDat <- apply(formattedIncDat[,colnames(formattedIncDat) != "Total"],2,function(x) x/sum(x))

## MAKE PLOTS
######################
## Microcephaly
######################
x <- (unique(brazilMicroDat$startDay) + unique(brazilMicroDat$endDay))/2 - 365
labels <- generate_plot_date_labels(brazilMicroDat)
breaks <- x[seq(1,length(x),by=4)]
labels <- labels[seq(1,length(x),by=4)]
birth_dat <- brazilMicroDat
birth_dat$startDay <- x

## State plots
micro_data_plot <- ggplot() + 
  geom_line(data=birth_dat,aes(x=startDay,y=microCeph,col=local)) + 
  facet_wrap(~local,ncol=4,scales="free_y") + 
  theme_bw() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  xlab("Date") +
  ylab("Number of suspected microcephaly cases per month")

## Total plot
total_micro <- data.frame("time"=unique(birth_dat$startDay),"micro"=microDat[,"Total"])
birthPlot <- ggplot() + 
  geom_line(data=total_micro,aes(x=time,y=micro)) + 
  theme_bw() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  xlab("") +
  ylab("Number of reported\n microcephaly cases per month")

    

######################
## Incidence
######################
months <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
labels <- unlist(lapply(c("2015","2016"), function(x) paste(months,x,sep=" ")))
breaks <- x[seq(1,length(x),by=1)]
labels <- labels[seq(1,length(x),by=1)]
inc_dat <- incDat
inc_dat$startDay <- (unique(inc_dat$startDay) + unique(inc_dat$endDay))/2 - 730

## State plots
inc_data_plot <- ggplot() + 
  geom_line(data=inc_dat,aes(x=startDay,y=inc,col=local)) + 
  facet_wrap(~local,ncol=5,scales="free_y") + 
  theme_bw() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  xlab("Date") +
  ylab("Number of suspected Zika cases per week")

## Total plot
total_inc <- data.frame("time"=unique(inc_dat$startDay),"inc"=formattedIncDat[,"Total"])
total_inc_plot <- ggplot() + 
  geom_line(data=total_inc,aes(x=time,y=inc)) + 
  theme_bw() + 
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  scale_x_continuous(breaks=breaks,labels=labels) +
  xlab("") +
  ylab("Number of suspected ZIKV\n cases per week")


colnames(inc_dat)[7] <- "State"
inc_dat$inc <- inc_dat$inc*100000/inc_dat$N_H
inc_dat[inc_dat$State == "bahia","inc"] <- inc_dat[inc_dat$State == "bahia","inc"]/10

p2 <- ggplot() + 
  geom_line(data=inc_dat,aes(x=startDay,y=inc,col=State)) + 
  geom_point(data=inc_dat,aes(x=startDay,y=inc,col=State),size=0.5) + 
  theme_bw() + 
  ylab("No. suspected ZIKV cases per 100,000*") +
  theme(legend.position="none",axis.text.x=element_text(angle=45, hjust=1)) +
  annotate("rect",xmin=cumsum(rep(getDaysPerMonth(),2)) - rep(getDaysPerMonth(),2),
           xmax=cumsum(rep(getDaysPerMonth(),2)),ymin=-1,ymax=3,alpha=0.4,fill=rep(c("grey","white"),12))+
  scale_x_continuous(breaks=breaks,labels=labels) +
  coord_cartesian(xlim=c(0,342.5),ylim=c(0,2.5))+
  xlab("Date")
  
  
if(SAVE_PLOT){
  png("micro_data_plot.png",height=600,width=600);plot(micro_data_plot);dev.off()
  svg("micro_data_plot.svg"); plot(micro_data_plot);dev.off()
  png("inc_data_plot.png",height=600,width=600); plot(inc_data_plot);dev.off()
  svg("inc_data_plot.svg"); plot(inc_data_plot);dev.off()
  png("micro_data_plot.png",height=600,width=600);plot(micro_data_plot);dev.off()
  svg("micro_data_plot.svg"); plot(micro_data_plot);dev.off()
  png("inc_data_plot.png",height=600,width=600); plot(inc_data_plot);dev.off()
  svg("inc_data_plot.svg"); plot(inc_data_plot);dev.off()
  png("total_micro_plot.png",height=600,width=600);plot(birthPlot);dev.off()
  svg("total_micro_plot.svg"); plot(birthPlot);dev.off()
  png("total_inc_plot.png",height=600,width=600);plot(total_inc_plot);dev.off()
  svg("total_inc_plot.svg"); plot(total_inc_plot);dev.off()
}
stats_table <- NULL

# Correlations ------------------------------------------------------------
###########################################
## 2. DATA CORRELATIONS
###########################################
## a.
## Pearson correlations for all state pairs
## Only take lower triangle
microCor <- cor(microDat, method="pearson")[lower.tri(cor(microDat[,colnames(microDat) != "Total"]))]
incCor <- cor(formattedIncDat, method="pearson")[lower.tri(cor(formattedIncDat[,colnames(formattedIncDat) != "Total"]))]

log_microCor <- cor(log(microDat[,colnames(microDat) != "Total"]+1), method="pearson")[lower.tri(cor(microDat[,colnames(microDat) != "Total"]))]
log_incCor <- cor(log(formattedIncDat[,colnames(formattedIncDat) != "Total"]+1), method="pearson")[lower.tri(cor(formattedIncDat[,colnames(formattedIncDat) != "Total"]))]

normalised_microCor <- cor(normalised_microDat, method="pearson")[lower.tri(cor(normalised_microDat))]
normalised_incCor <- cor(normalised_incDat, method="pearson")[lower.tri(cor(normalised_incDat))]
                                                       
## b.
## Get quantiles and median
microCor_range <- quantile(c(microCor),c(0.025,0.5,0.975))
incCor_range <- quantile(c(incCor),c(0.025,0.5,0.975))

log_microCor_range <- quantile(c(log_microCor),c(0.025,0.5,0.975))
log_incCor_range <- quantile(c(log_incCor),c(0.025,0.5,0.975))

normalised_microCor_range <- quantile(c(normalised_microCor),c(0.025,0.5,0.975))
normalised_incCor_range <- quantile(c(normalised_incCor),c(0.025,0.5,0.975))


print("Microcephaly correlation stats:")
print(c("mean"=mean(c(microCor)),microCor_range))
print("ZIKV incidence stats:")
print(c("mean"=mean(c(incCor)),incCor_range))

print("Log microcephaly correlation stats:")
print(c("mean"=mean(c(log_microCor)),log_microCor_range))
print("Log ZIKV incidence stats:")
print(c("mean"=mean(c(log_incCor)),log_incCor_range))

print("Normalised microcephaly correlation stats:")
print(c("mean"=mean(c(normalised_microCor)),normalised_microCor_range))
print("Normalised ZIKV incidence stats:")
print(c("mean"=mean(c(normalised_incCor)),normalised_incCor_range))

stats <- data.frame("Pearson_Correlation"=c("Microcephaly_cor","ZIKV_cor","logMicrocephaly_cor","logZIKV_cor","normMicrocephaly_cor","normZIKV_cor"),
                    "mean"=c(mean(c(microCor)),mean(c(incCor)),mean(c(log_microCor)),mean(c(log_incCor)),mean(c(normalised_microCor)),mean(c(normalised_incCor))),
                    "lower"=c(microCor_range[1],incCor_range[1],log_microCor_range[1],log_incCor_range[1],normalised_microCor_range[1],normalised_incCor_range[1]),
                    "median"=c(microCor_range[2],incCor_range[2],log_microCor_range[2],log_incCor_range[2],normalised_microCor_range[2],normalised_incCor_range[2]),
                    "upper"=c(microCor_range[3],incCor_range[3],log_microCor_range[3],log_incCor_range[3],normalised_microCor_range[3],normalised_incCor_range[3]))
            
                    

## Convert correlation matrices to vectors
microCor_melt <- c(microCor)
incCor_melt <- c(incCor)

## c.
if(SAVE_PLOT) png("raw_corr_boxplot.png",height=600,width=600)
boxplot(microCor_melt,incCor_melt,names=c("Microcephaly","ZIKV"),ylab="Pearson Correlation Coefficient")
if(SAVE_PLOT) dev.off()
if(SAVE_PLOT) svg("raw_corr_boxplot.svg")
boxplot(microCor_melt,incCor_melt,names=c("Microcephaly","ZIKV"),ylab="Pearson Correlation Coefficient")
if(SAVE_PLOT) dev.off()

# Cross Correlations ------------------------------------------------------


###########################
## d. CROSS CORRELATIONS
###########################
result <- ccf(weekly_microceph$microceph,weekly_inc$inc,lag.max=100)
#result$lag <- result$lag + 20
if(SAVE_PLOT) png("cross_corr.png",height=600,width=600)
plot(result,main="Cross-correlation between microcephaly and ZIKV incidence curves",xlab="Lag (weeks)")
if(SAVE_PLOT) dev.off()
if(SAVE_PLOT) svg("cross_corr.svg")
plot(result,main="Cross-correlation between microcephaly and ZIKV incidence curves",xlab="Lag (weeks)")
if(SAVE_PLOT) dev.off()

# Spatio Correlations -----------------------------------------------------
##################################
## e. SPATIO-TEMPORAL CORRELATIONS
##################################
## Get country data
map_df <- generate_country_boundaries()
brazil1 <- raster::getData("GADM",country="Brazil",level=1)
map <- fortify(brazil1)
map$id <- as.integer(map$id)
dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
#dat$state <- get_states()[as.character(dat$state)]
# Finding centroids of each state as location for labels
centers <- data.frame(rgeos::gCentroid(brazil1, byid = TRUE))
centers$state <- dat$state

## Convert to factors and order
#centers$State <- convert_name_to_state_factor(centers$state)
#order <- get_correct_order()
#centers$state <- factor(centers$state,levels=order)
#centers <- centers[order(centers$state),]
#centers <- centers[centers$state %in% unique(colnames(microDat)),]

##############################
## Microcephaly Data
##############################
transpose_micro <- t(microDat[,colnames(microDat) != "Total"])
if(LOG_CORRELATION) transpose_micro <- t(log(microDat[,colnames(microDat) != "Total"]+1))
normalised_transpose_micro <- t(apply(transpose_micro,1,function(x) x/sum(x)))
filename <- "micro_spatio"
if(NORMALISE){ transpose_micro <- normalised_transpose_micro; filename <- "micro_spatio_normalised"}

micro_result <- Sncf(centers$x,centers$y,z=transpose_micro,w=transpose_micro,latlon=TRUE,quiet=TRUE,df=10)
## Lag by one month
k=1; transpose_micro_lag1 <- t(apply(transpose_micro, 1, function(x) c(rep(NA, k), x)[1 : length(x)]))
micro_result_lag1 <- Sncf(centers$x,centers$y,z=transpose_micro,w=transpose_micro_lag1,latlon=TRUE,na.rm=TRUE,quiet=TRUE,df=10)
## Lag by two months
k=2; transpose_micro_lag2 <- t(apply(transpose_micro, 1, function(x) c(rep(NA, k), x)[1 : length(x)]))
micro_result_lag2 <- Sncf(centers$x,centers$y,z=transpose_micro,w=transpose_micro_lag2,latlon=TRUE,quiet=TRUE,df=10)

if(SAVE_PLOT){
  png(paste(filename,".png",sep=""),height=300,width=700)
  par(mfrow=c(1,3),mar=c(2,2,2,2),cex=1.5)
  plot(micro_result$real$predicted$y~micro_result$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result$boot$boot.summary$predicted$x, rev(micro_result$boot$boot.summary$predicted$x)), 
          c(micro_result$boot$boot.summary$predicted$y[2,],rev(micro_result$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result$real$predicted$y~micro_result$real$predicted$x)
  
  plot(micro_result_lag1$real$predicted$y~micro_result_lag1$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result_lag1$boot$boot.summary$predicted$x, rev(micro_result_lag1$boot$boot.summary$predicted$x)), 
          c(micro_result_lag1$boot$boot.summary$predicted$y[2,],rev(micro_result_lag1$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result_lag1$real$predicted$y~micro_result_lag1$real$predicted$x)
  
  plot(micro_result_lag2$real$predicted$y~micro_result_lag2$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result_lag2$boot$boot.summary$predicted$x, rev(micro_result_lag2$boot$boot.summary$predicted$x)), 
          c(micro_result_lag2$boot$boot.summary$predicted$y[2,],rev(micro_result_lag2$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result_lag2$real$predicted$y~micro_result_lag2$real$predicted$x)
  dev.off()
  
  svg(paste(filename,".svg",sep=""),height=4,width=8)
  par(mfrow=c(1,3),mar=c(2,2,2,2),cex=1.5)
  plot(micro_result$real$predicted$y~micro_result$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result$boot$boot.summary$predicted$x, rev(micro_result$boot$boot.summary$predicted$x)), 
          c(micro_result$boot$boot.summary$predicted$y[2,],rev(micro_result$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result$real$predicted$y~micro_result$real$predicted$x)
  
  plot(micro_result_lag1$real$predicted$y~micro_result_lag1$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result_lag1$boot$boot.summary$predicted$x, rev(micro_result_lag1$boot$boot.summary$predicted$x)), 
          c(micro_result_lag1$boot$boot.summary$predicted$y[2,],rev(micro_result_lag1$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result_lag1$real$predicted$y~micro_result_lag1$real$predicted$x)
  
  plot(micro_result_lag2$real$predicted$y~micro_result_lag2$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(micro_result_lag2$boot$boot.summary$predicted$x, rev(micro_result_lag2$boot$boot.summary$predicted$x)), 
          c(micro_result_lag2$boot$boot.summary$predicted$y[2,],rev(micro_result_lag2$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(micro_result_lag2$real$predicted$y~micro_result_lag2$real$predicted$x)
  dev.off()
}


##############################
## Incidence Data
##############################
## Get country data
map_df <- generate_country_boundaries()
brazil1 <- raster::getData("GADM",country="Brazil",level=1)
map <- fortify(brazil1)
map$id <- as.integer(map$id)
dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
#dat$state <- get_states()[as.character(dat$state)]
# Finding centroids of each state as location for labels
centers <- data.frame(gCentroid(brazil1, byid = TRUE))
centers$state <- dat$state

## Convert to factors and order
#centers$state <- convert_name_to_state_factor(centers$state)
#order <- get_correct_order()
#centers$state <- factor(centers$state,levels=order)
#centers <- centers[order(centers$state),]
#centers <- centers[centers$state %in% unique(colnames(formattedIncDat)),]

transpose_inc <- t(formattedIncDat[,colnames(formattedIncDat) != "Total"])
if(LOG_CORRELATION) transpose_inc <- t(log(formattedIncDat[,colnames(formattedIncDat) != "Total"]+1))
normalised_transpose_inc <- t(apply(transpose_inc,1,function(x) x/sum(x)))
filename <- "inc_spatio"
if(NORMALISE){ transpose_inc <- normalised_transpose_inc; filename <- "micro_spatio_inc"}

inc_result <- Sncf(centers$x,centers$y,w=transpose_inc,z=transpose_inc,latlon=TRUE,quiet=TRUE,df=10)
## Lag by one week
k=1
transpose_inc_lag1 <- t(apply(transpose_inc, 1, function(x) c(rep(NA, k), x)[1 : length(x)]))
inc_result_lag1 <- Sncf(centers$x,centers$y,z=transpose_inc,w=transpose_inc_lag1,latlon=TRUE,na.rm=TRUE,quiet=TRUE,df=10)

## Lag by two week
k=2
transpose_inc_lag2 <- t(apply(transpose_inc, 1, function(x) c(rep(NA, k), x)[1 : length(x)]))
inc_result_lag2 <- Sncf(centers$x,centers$y,z=transpose_inc,w=transpose_inc_lag2,latlon=TRUE,quiet=TRUE,df=10)


if(SAVE_PLOT){
  png(paste(filename,".png",sep=""),height=300,width=700)
  par(mfrow=c(1,3),mar=c(2,2,2,2),cex=1.5)
  plot(inc_result$real$predicted$y~inc_result$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result$boot$boot.summary$predicted$x, rev(inc_result$boot$boot.summary$predicted$x)), 
          c(inc_result$boot$boot.summary$predicted$y[2,],rev(inc_result$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result$real$predicted$y~inc_result$real$predicted$x)
  
  plot(inc_result_lag1$real$predicted$y~inc_result_lag1$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result_lag1$boot$boot.summary$predicted$x, rev(inc_result_lag1$boot$boot.summary$predicted$x)), 
          c(inc_result_lag1$boot$boot.summary$predicted$y[2,],rev(inc_result_lag1$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result_lag1$real$predicted$y~inc_result_lag1$real$predicted$x)
  
  plot(inc_result_lag2$real$predicted$y~inc_result_lag2$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result_lag2$boot$boot.summary$predicted$x, rev(inc_result_lag2$boot$boot.summary$predicted$x)), 
          c(inc_result_lag2$boot$boot.summary$predicted$y[2,],rev(inc_result_lag2$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result_lag2$real$predicted$y~inc_result_lag2$real$predicted$x)
  dev.off()
  
  svg(paste(filename,".svg",sep=""),height=4,width=8)
  par(mfrow=c(1,3),mar=c(2,2,2,2),cex=1.5)
  plot(inc_result$real$predicted$y~inc_result$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result$boot$boot.summary$predicted$x, rev(inc_result$boot$boot.summary$predicted$x)), 
          c(inc_result$boot$boot.summary$predicted$y[2,],rev(inc_result$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result$real$predicted$y~inc_result$real$predicted$x)
  
  plot(inc_result_lag1$real$predicted$y~inc_result_lag1$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result_lag1$boot$boot.summary$predicted$x, rev(inc_result_lag1$boot$boot.summary$predicted$x)), 
          c(inc_result_lag1$boot$boot.summary$predicted$y[2,],rev(inc_result_lag1$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result_lag1$real$predicted$y~inc_result_lag1$real$predicted$x)
  
  plot(inc_result_lag2$real$predicted$y~inc_result_lag2$real$predicted$x,ylim=c(-0.2,1),type='l',ylab="",xlab="")
  polygon(c(inc_result_lag2$boot$boot.summary$predicted$x, rev(inc_result_lag2$boot$boot.summary$predicted$x)), 
          c(inc_result_lag2$boot$boot.summary$predicted$y[2,],rev(inc_result_lag2$boot$boot.summary$predicted$y[10,])),
          col="grey75",border=NA)
  abline(h=0,lty=3)
  lines(inc_result_lag2$real$predicted$y~inc_result_lag2$real$predicted$x)
  dev.off()
}

# Correlation Plots -------------------------------------------------------


##################################
## f. CORRELATION PLOTS
##################################
## Incidence unadjusted
col1 <- colorRampPalette(c("blue","steelblue1","white","yellow","red"))
use_inc <- formattedIncDat
filename <- "inc_corr"
if(NORMALISE){
  filename <- "inc_corr_normalised"
  use_inc <- normalised_incDat
}

inc_corr_matrix <- rcorr(as.matrix(use_inc[,colnames(use_inc) != "Total"]),type="pearson")
## Bonferroni correction
m <- ncol(use_inc) -1
m <- (m^2 - m)/2
inc_corr_matrix$P <- inc_corr_matrix$P * m
png(paste(filename,".png",sep=""),height=600,width=600)
corrplot(inc_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=inc_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()
svg(paste(filename,".svg",sep=""))
corrplot(inc_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=inc_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()

## Inc Log values
inc_corr_matrix <- rcorr(log(as.matrix(formattedIncDat[,colnames(formattedIncDat) != "Total"])+1),type="pearson")
## Bonferroni correction
m <- ncol(formattedIncDat) -1
m <- (m^2 - m)/2
inc_corr_matrix$P <- inc_corr_matrix$P * m
png("inc_corr_log.png",height=600,width=600)
corrplot(inc_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=inc_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()
svg("inc_corr_log.svg")
corrplot(inc_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=inc_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()

##################################
## Microceph unadjusted
use_micro <- microDat
filename <- "micro_corr"
if(NORMALISE){
  filename <- "micro_corr_normalised"
  use_micro <- normalised_microDat
}

micro_corr_matrix <- rcorr(as.matrix(use_micro[,colnames(use_micro) != "Total"]),type="pearson")
## Bonferroni correction
m <- ncol(use_micro) -1
m <- (m^2 - m)/2
micro_corr_matrix$P <- micro_corr_matrix$P * m
png(paste(filename,".png",sep=""),height=600,width=600)
corrplot(micro_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=micro_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()
svg(paste(filename,".svg",sep=""))
corrplot(micro_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=micro_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()

## Microceph Log values
micro_corr_matrix <- rcorr(log(as.matrix(microDat[,colnames(microDat) != "Total"])+1),type="pearson")
## Bonferroni correction
m <- ncol(microDat) -1
m <- (m^2 - m)/2
micro_corr_matrix$P <- micro_corr_matrix$P * m
png("micro_corr_log.png",height=600,width=600)
corrplot(micro_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=micro_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()
svg("micro_corr_log.svg")
corrplot(micro_corr_matrix$r,method="shade",diag = FALSE,tl.col = "black",p.mat=micro_corr_matrix$P,sig.level=0.05,col=col1(20))
dev.off()

# Peak times --------------------------------------------------------------
###########################################
## 3. PEAK TIME SUMMARIES
###########################################
# Get peak times for microcephaly data
peakTimes <- NULL
locals <- NULL
times <- brazilMicroDat[brazilMicroDat$local=="bahia",c("startDay","endDay")]

for(local in sort(unique(brazilMicroDat$local))){
 tmp <- brazilMicroDat[brazilMicroDat$local==local,c("startDay","endDay","microCeph")]
 peakRow <- which.max(tmp[,"microCeph"])
 peakTime <- mean(as.numeric(times[peakRow,c("startDay","endDay")]))
 peakTimes <- c(peakTimes,peakTime)
 locals <- c(locals, local)
}
total_peak <- mean(as.numeric(times[which.max(microDat[,"Total"]),c("startDay","endDay")]))
peakTimes <- c(peakTimes,total_peak)
locals <- c(locals,"Total")
############
# MICRO RESULT
peakTimes_micro <- data.frame("state"=locals,"peakTime"=peakTimes)
############

#####################
# Get peak times for incidence
peakTimes <- NULL
locals <- NULL
times <- incDat[incDat$local=="bahia",c("startDay","endDay")]
for(local in sort(unique(incDat$local))){
  tmp <- incDat[incDat$local==local,c("startDay","endDay","inc")]
  peakRow <- which.max(tmp[,"inc"])
  peakTime <- mean(as.numeric(times[peakRow,c("startDay","endDay")]))
  peakTimes <- c(peakTimes,peakTime)
  locals <- c(locals, local)
}
total_peak <- mean(as.numeric(times[which.max(formattedIncDat[,"Total"]),c("startDay","endDay")]))
peakTimes <- c(peakTimes,total_peak)
locals <- c(locals,"Total")
############
# INC RESULT
peakTimes_inc <- data.frame("state"=locals,"peakTime"=peakTimes)
############

peakTimes <- merge(peakTimes_inc,peakTimes_micro,by="state",all=TRUE)
colnames(peakTimes) <- c("state","ZIKV peak","microcephaly peak")
print("ZIKV peak statistics:")
print(as.Date(c("mean"=mean(na.omit(peakTimes[,"ZIKV peak"])),quantile(na.omit(peakTimes[,"ZIKV peak"]),c(0.025,0.5,0.975))),origin="2013-01-01"))

print("Microcephaly peak statistics:")
print(as.Date(c("mean"=mean(na.omit(peakTimes[,"microcephaly peak"])),quantile(na.omit(peakTimes[,"microcephaly peak"]),c(0.025,0.5,0.975))),origin="2013-01-01"))

# Max cases ---------------------------------------------------------------
###########################################
## 4. MAX CASES
###########################################
# Max cases
max_micro_cases <- as.numeric(colSums(microDat[,colnames(microDat) != "Total"]))
names(max_micro_cases) <- colnames(microDat[,colnames(microDat) != "Total"])
print("Max microcephaly cases:")
print(c("mean"=mean(max_micro_cases),quantile(max_micro_cases,c(0.025,0.5,0.975))))
print(summary(max_micro_cases))

max_inc_cases <- as.numeric(colSums(formattedIncDat[,colnames(formattedIncDat) != "Total"]))
names(max_inc_cases) <- colnames(formattedIncDat[,colnames(formattedIncDat) != "Total"])
print("Max ZIKV cases:")
print(c("mean"=mean(max_inc_cases),quantile(max_inc_cases,c(0.025,0.5,0.975))))
print(summary(max_inc_cases))

max_cases <- t(plyr::rbind.fill.matrix("ZIKV"=t(max_inc_cases),"Microcephaly"=t(max_micro_cases)))
colnames(max_cases) <- c("ZIKV","Microcephaly")
max_cases <- data.frame("state"=rownames(max_cases),max_cases)

# Table -------------------------------------------------------------------
#############################################
## 5. TABLE OF STATS
#############################################
all_stats <- merge(peakTimes,max_cases,by="state",all=TRUE)
colnames(all_stats) <- c("State","Peak Time ZIKV","Peak Time Microcephaly","Total ZIKV","Total Microcephaly")
tmp_stats <- all_stats[all_stats$State != "Total",]

mean_peak_inc <- as.Date(mean(na.omit(tmp_stats[,2])),"2013-01-1")
quantiles_peak_inc <- as.Date(quantile(na.omit(tmp_stats[,2]),c(0.025,0.5,0.975)),"2013-01-1")
mean_peak_micro <- as.Date(mean(na.omit(tmp_stats[,3])),"2013-01-1")
quantiles_peak_micro <- as.Date(quantile(na.omit(tmp_stats[,3]),c(0.025,0.5,0.975)),"2013-01-1")

mean_inc <- mean(na.omit(tmp_stats[,4]))
mean_micro <- mean(na.omit(tmp_stats[,5]))
quantiles_inc <- quantile(na.omit(tmp_stats[,4]),c(0.025,0.5,0.975))
quantiles_micro <- quantile(na.omit(tmp_stats[,5]),c(0.025,0.5,0.975))

all_stats[,2] <- as.Date(all_stats[,2],"2013-01-01")
all_stats[,3] <- as.Date(all_stats[,3],"2013-01-01")

write.table(all_stats,"data_stats.csv",sep=",",row.names=FALSE)

if(SAVE_PLOT){
  png("peaktime_box.png",width=600,height=500); boxplot(all_stats[,c(2,3)]); dev.off()
  svg("peaktime_box.svg"); boxplot(all_stats[,c(2,3)]); dev.off()
}


