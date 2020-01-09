######################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This script produces model trajectories from MCMC fits
## It also produces Figure 2 from the manuscript.
## The user must have available the MCMC chains for
## each of the fits below (see scripts/job_submissions)
## Users should also check plots/google_searches.R
## to produce the top half of figure 2.
library(zikaInfer)
library(cowplot)
library(zoo)
library(ggplot2)
library(extrafont)
library(coda)

## Where are the MCMC chains stored
#topDir <- "~/net/home/zika/outputs/"
topDir <- "/media/james/JH USB/zika/outputs"
topDir <- "H:/all_results_backup/zika/zika_usb/outputs/"


###################################################
## NEJM Northeast Brazil model fit
chainWD <- paste0(topDir, "/northeast_deOliveira2017/northeast/model_1")
parTab <- read_inipars(chainWD)
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,750000,TRUE,FALSE,FALSE)[["chain"]]
datFile = "~/Documents/Zika/Data/brazil/deOliveira2017/microCeph_deOliveira2017_clean.csv"
incFile = "~/Documents/Zika/Data/brazil/deOliveira2017/zikv_inc_deOliveira2017_clean.csv"
local = "northeast"
localName = "Northeast Brazil NEJM"
incScale=0.015
runs=200
ylim=0.01
NE_plot <- indiv_model_fit(datFile,incFile, local, localName,
                           incScale, runs, ylim,xlim=730,bot=FALSE,standalone=FALSE,
                            parTab=parTab, chain=chain, forecast=FALSE, weeks=FALSE)
NE_plot

###################################################
## Pernambuco Brazil model fit
chainWD <- paste0(topDir, "pernambuco_peak_confirmed/pernambuco/model_1")
parTab <- read_inipars(chainWD)
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,750000,TRUE,FALSE,FALSE)[["chain"]]
datFile = "~/Documents/Zika/Data/brazil/microceph_reports_2016_confirmed.csv"
incFile = NULL
local = "pernambuco"
localName = "Pernambuco"
incScale=10
runs=200
ylim=0.02
pern_plot <- indiv_model_fit(datFile,incFile, local, localName,
                           incScale, runs, ylim,xlim=0,bot=FALSE,standalone=FALSE,
                           parTab=parTab, chain=chain, forecast=FALSE, weeks=TRUE)
pern_plot


###################################################
## Bahia reports data fit
chainWD = paste0(topDir,"bahia/bahia/model_1")
datFile = "H:/all_results_backup/zika/Data/brazil/microceph_reports_2016.csv"
incFile = "H:/all_results_backup/zika/Data/brazil/brazil_report_zikv_inc_2016.csv"
parTab <- read_inipars(chainWD)
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,750000,TRUE,FALSE,FALSE)[["chain"]]
local = "bahia"
localName = "Bahia - reports"
incScale=0.015
runs= 100
ylim <- 0.03

bahia_plot <- indiv_model_fit(datFile,incFile, local, localName,
                      incScale, runs, ylim,xlim=730,bot=FALSE,standalone=FALSE,
                      parTab=parTab, chain=chain, forecast=FALSE, weeks=FALSE)
bahia_plot[[1]]

## Generate residual plot
all_dat <- bahia_plot[[2]]

inc_preds <- all_dat$incPredictions
inc_preds$index <- rep(1:200, each=length(unique(inc_preds$time)))
inc_dat <- all_dat$incDat
inc_dat <- inc_dat[,c("meanDay","inc","N_H")]
colnames(inc_dat) <- c("time","obs","N_H")
inc_dat$obs <- inc_dat$obs/inc_dat$N_H
inc_preds <- merge(inc_preds, inc_dat)
inc_preds$residual <- inc_preds$obs - inc_preds$inc

inc_preds_bounds <- plyr::ddply(inc_preds, ~time, function(x) quantile(x$residual,c(0.025,0.5,0.975)))
colnames(inc_preds_bounds) <- c("time","lower","median","upper")
inc_preds_bounds$dat <- "ZIKV infection"

micro_preds <- all_dat$microPredictions
micro_preds$index <- rep(1:200, each=length(unique(micro_preds$time)))
micro_dat <- all_dat$data
micro_dat <- micro_dat[,c("meanDay","microCeph","births")]
colnames(micro_dat) <- c("time","obs","births")
micro_dat$obs <- micro_dat$obs/micro_dat$births
micro_preds <- merge(micro_preds, micro_dat)
micro_preds$residual <- micro_preds$obs - micro_preds$micro

micro_preds_bounds <- plyr::ddply(micro_preds, ~time, function(x) quantile(x$residual,c(0.025,0.5,0.975)))
colnames(micro_preds_bounds) <- c("time","lower","median","upper")
micro_preds_bounds$dat <- "Microcephaly"

all_bounds <- rbind(inc_preds_bounds, micro_preds_bounds)
options(scipen = 999)


labels <- rep(getDaysPerMonth(3),4)
labels <- c(0,cumsum(labels))
labels_names <- as.yearmon(as.Date(labels,origin="2013-01-01"))
x_lower <- 730

p1 <- ggplot(all_bounds) + 
  geom_hline(yintercept=0,linetype="dashed",col="grey70") +
  geom_ribbon(aes(x=time,ymax=upper,ymin=lower,fill=dat),alpha=0.3) + 
  geom_line(aes(x=time,y=median,col=dat)) +
  ylab("Residual per birth/capita incidence\n(observed - predicted)") +
  xlab("Date") +
  scale_fill_manual(values=c("blue","forestgreen")) +
  scale_colour_manual(values=c("blue","forestgreen")) +
  scale_x_continuous(limits=c(x_lower,max(labels)),breaks=labels,labels=labels_names) +
  facet_wrap(~dat,scales="free_y",ncol=1) +
  theme_classic() +
  theme(legend.position="none",
        axis.text=element_text(size=6),
        axis.title = element_text(size=8),
        strip.text=element_text(size=8))

p2 <- ggplot(all_bounds) + 
  geom_histogram(aes(x=median,fill=dat),bins=25) + 
  geom_vline(xintercept=0,linetype="dashed",col="grey70") +
  #ylab("Residual per birth/capita incidence\n(observed - predicted)") +
  ylab("Count") +
  xlab("Residual per birth/capita incidence\n(observed - predicted)") +
  scale_fill_manual(values=c("blue","forestgreen")) +
  scale_colour_manual(values=c("blue","forestgreen")) +
  scale_x_continuous(breaks= scales::pretty_breaks(5)) +
  #scale_x_continuous(limits=c(x_lower,max(labels)),breaks=labels,labels=labels_names) +
  facet_wrap(~dat,scales="free",ncol=1) +
  theme_classic() +
  theme(legend.position="none",
        axis.text=element_text(size=6),
        axis.title = element_text(size=8),
        strip.text=element_text(size=8))

pdf("E:/James/Google Drive/Thesis/Corrections/residual_plot.pdf",width=7.5,height=5)
p1 + p2 + plot_layout(widths=c(1.5,1))
dev.off()
png("E:/James/Google Drive/Thesis/Corrections/residual_plot.pdf",width=7.5,height=5,res=300,units="in")
p1 + p2 + plot_layout(widths=c(1.5,1))
dev.off()


inc_preds <- reshape2::dcast(inc_preds, time~index,value.var="inc")


###################################################
## Rio Grande Do Norte reports data fit
chainWD = paste0(topDir,"riograndedonorte_confirmed/riograndedonorte/model_1")
datFile = "~/Documents/Zika/Data/brazil/microceph_reports_2016_confirmed.csv"
incFile = "~/Documents/Zika/Data/brazil/brazil_report_zikv_inc_2016.csv"
parTab <- read_inipars(chainWD)
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,750000,TRUE,FALSE,FALSE)[["chain"]]
local = "riograndedonorte"
localName = "Rio Grande do Norte - reports"
incScale=0.008
runs=100
ylim <- 0.03

rio_plot <- indiv_model_fit(datFile,incFile, local, localName,
                              incScale, runs, ylim,xlim=730,bot=TRUE,standalone=FALSE,
                              parTab=parTab, chain=chain, forecast=FALSE, weeks=FALSE)
rio_plot
###################################################
## Colombia weekly inc fit
chainWD = paste0(topDir,"colombia_inc_confirmed/colombia/model_1")
datFile = "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly_2017.csv"
incFile = "~/Documents/Zika/Data/colombia/zikv_inc_2017_confirmed.csv"
parTab <- read_inipars(chainWD)
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,750000,TRUE,FALSE,FALSE)[["chain"]]
local = "colombia"
localName = "Colombia - reports"
incScale=0.03
runs=100
ylim <- 0.005

colombia_plot <- indiv_model_fit(datFile,incFile, local, localName,
                            incScale, runs, ylim,xlim=730,bot=FALSE,standalone=FALSE,
                            parTab=parTab, chain=chain, forecast=FALSE, weeks=FALSE)
colombia_plot


comb_plot <- plot_grid(NE_plot,colombia_plot,bahia_plot,rio_plot,ncol=1,align="v",rel_heights=c(1,1,1,1.3))


###################################################
## FORECAST PLOTS
###################################################
###################################################
## Bahia FORECAST reports data fit
topDir <- "/media/james/JH USB/zika/forecast_all/"
chainWD = paste0(topDir,"bahia_forecast")
datFile = "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile <- "~/Documents/Zika/Data/brazil/zika_inc_reports.csv"
parTab <- read.csv("~/Documents/Zika/zikaInfer/scripts/job_submissions/parTab_forecast.csv")
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,50000,TRUE,FALSE,TRUE)[["chain"]]
local = "bahia"
localName = "Bahia - reports"
incScale=0.015
runs= 100
ylim <- 0.03
bahia_forecast <- indiv_model_fit(datFile,incFile, local, localName,
                              incScale, runs, ylim,xlim=730,bot=FALSE,standalone=TRUE,
                              parTab=parTab, chain=chain, forecast=TRUE, forecastPostChange=TRUE, weeks=FALSE) +
    theme(axis.text.x=element_text(angle=0,hjust=0.5),                                                                            
          legend.position="top")

###################################################
## NEJM Northeast Brazil FORECAST fit
#topDir <- "~/net/home/zika/outputs/"
chainWD <- paste0(topDir, "northeast_forecast")
parTab <- read.csv("~/Documents/Zika/zikaInfer/scripts/job_submissions/parTab_forecast.csv")
chain <- lazymcmc::load_mcmc_chains(chainWD,parTab,FALSE,1,50000,TRUE,FALSE,FALSE)[["chain"]]
datFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
incFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
local = "northeast"
localName = "Northeast Brazil NEJM"
incScale=0.015
runs=200
ylim=0.03
NE_forecast <- indiv_model_fit(datFile,incFile, local, localName,
                           incScale, runs, ylim,xlim=730,bot=FALSE,standalone=TRUE,
                           parTab=parTab, chain=chain, forecast=TRUE, forecastPostChange=TRUE, weeks=FALSE)

####################
### SAVE PLOTS
####################
png("modelfits.png",width=4,height=6,units="in",res=300)
print(comb_plot)
dev.off()

cairo_ps("modelfits.eps",width=4,height=6,family="Arial")
print(comb_plot)
dev.off()

png("bahia_forecast.png",width=5,height=3,units="in",res=300)
print(bahia_forecast)
dev.off()

cairo_ps("bahia_forecast.eps",width=5,height=3,family="Arial")
print(bahia_forecast)
dev.off()

png("NE_forecast.png",width=5,height=3,units="in",res=300)
print(NE_forecast)
dev.off()

cairo_ps("NE_forecast.eps",width=5,height=3,family="Arial")
print(NE_forecast)
dev.off()


####################
## ESSAY FIGURE 2
## Need forecast_plots.R too
fig2 <- plot_grid(bahia_forecast,trend_plot,ncol=1,align="hv",rel_heights = c(1.3,1))
png("Fig2.png",width=6,height=6,units="in",res=300)
print(fig2)
dev.off()

svg("Fig2.svg",width=7.5,height=6,family="Arial")
print(fig2)
dev.off()

cairo_ps("Fig2.eps",width=7.5,height=6,family="Arial")
print(fig2)
dev.off()
