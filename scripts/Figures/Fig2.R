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


###################################################
## NEJM Northeast Brazil model fit
chainWD <- paste0(topDir, "northeast_deOliveira2017/northeast/model_1")
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
datFile = "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile = "~/Documents/Zika/Data/brazil/brazil_report_zikv_inc_2016.csv"
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
bahia_plot
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
