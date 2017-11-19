#source("~/Documents/Zika/plots/plottingScripts/plotting_help_functions.R")
#source("~/Documents/Zika/plots/plottingScripts/stat_help.R")
library(zikaProj)
library(cowplot)
library(zoo)
library(ggplot2)
library(extrafont)

chainWD = "~/Documents/Zika/28.02.2017_chains/multi_all1/model_1"
datFile = "~/Documents/Zika/Data/allDat28.02.17.csv"
incDatFile = "~/Documents/Zika/Data/inc_data_120317.csv"
runs=200
datCutOff=1214
incScale=0.1
p <- main_model_fits(chainWD,datFile,incDatFile,runs,datCutOff,incScale)


topDir <- "/media/james/JH USB/outputs/"
chainWD <- paste0(topDir, "northeast_forecast")
datFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_microceph.csv"
incFile = "~/Documents/Zika/Data/brazil/Northeast/northeast_zikv.csv"
local = "northeast"
localName = "Northeast Brazil NEJM"
incScale=0.01
runs=200
ylim=0.02
p1 <- indiv_model_fit(chainWD,datFile,incFile,local,
                      localName,incScale,runs,ylim,TRUE,xlim=730, parTab=parTab, chain=chain)
p1


chainWD = "~/Documents/Zika/28.02.2017_chains/bahia/bahia/model_1"
datFile = "~/Documents/Zika/Data/microceph_bahia_2016.csv"
incFile = "~/Documents/Zika/Data/inc_2016.csv"
local = "bahia"
localName = "Bahia - reports"
incScale=0.01
runs= 100
ylim <- 0.05

p1 <- zikaProj::indiv_model_fit(chainWD,datFile,incFile,local,localName,incScale,runs,ylim,TRUE,xlim=730, parTab=parTab)
p1



#chainWD = "~/Documents/Zika/28.02.2017_chains/colombia_inc/bahia/model_1"

datFile <- "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile <- "~/Documents/Zika/Data/brazil/zika_inc_reports.csv"

chainWD = "~/net/home/zika/outputs/bahia_early_all/"
#datFile = "~/Documents/Zika/Data/colombia/microcephaly_dat_weekly.csv.csv"
#incFile = "~/Documents/Zika/Data/colombia/zikv_inc_2016.csv"
local = "bahia"
localName = "bahia"
incScale=0.01
runs=100
ylim <- 0.02

p3 <- zikaProj::indiv_model_fit(chainWD,datFile,incFile,local,localName,incScale,runs,ylim,FALSE,xlim=730,forecast=TRUE)
p3

chainWD = "~/Documents/Zika/28.02.2017_chains/riograndedonorte/riograndedonorte/model_1"
datFile = "~/Documents/Zika/Data/microceph_bahia_2016.csv"
incFile = "~/Documents/Zika/Data/inc_2016.csv"
local = "riograndedonorte"
localName = "Rio Grande do Norte - reports"
incScale=0.01
runs=100
ylim <- 0.03

p4 <- indiv_model_fit(chainWD,datFile,incFile,local,localName,incScale,runs,ylim,bot=TRUE)
p4

chainWD = "~/Documents/Zika/28.02.2017_chains/pernambuco_peak/pernambuco/model_1"
datFile = "~/Documents/Zika/Data/microceph_bahia_2016.csv"
incFile = NULL
local = "pernambuco"
localName = "Pernambuco - reports"
incScale= 1
runs=100
ylim <- 0.1

p5 <- indiv_model_fit(chainWD,datFile,incFile,local,localName,incScale,runs,ylim)
p5

pwow <- plot_grid(p,plot_grid(p1,p2,p3,p5,p4,ncol=1,align="v",rel_heights=c(1,1,1,1,1.4)),rel_widths=c(2,1))
png("~/Documents/Zika/plots/manuscript/Fig5_modelfits.png",width=10,height=8,units="in",res=300)
pwow
dev.off()

cairo_ps("~/Documents/Zika/plots/manuscript/Fig5_modelfits.eps",width=10,height=8,family="Arial")
pwow
dev.off()
cairo_ps("~/Documents/Zika/plots/manuscript/Fig5_modelfits.pdf",width=10,height=8,family="Arial")
pwow
dev.off()


pwow1 <- plot_grid(p1,p3,p2,p4,ncol=2,align="v",rel_heights=c(1,1.4,1,1.4))
png("~/Documents/Zika/plots/manuscript/Fig5_modelfits_b.png",width=8,height=5,units="in",res=300)
pwow1
dev.off()

####################
## ESSAY FIGURE 2
## Need forecast_plots.R too
big_plot <- plot_grid(plot_grid(p1,p2,rel_heights=c(1,1.25),ncol=1,align="v"),plot_grid(pwow,pwow1,ncol=1),ncol=2,rel_widths=c(1.2,1))
png("~/Documents/Zika/plots/manuscript/Fig2_essay.png",width=8,height=5,units="in",res=300)
big_plot
dev.off()

cairo_ps("~/Documents/Zika/plots/manuscript/Fig2_essay.eps",width=8,height=5,family="Arial")
big_plot
dev.off()

cairo_ps("~/Documents/Zika/plots/manuscript/Fig2_essay.pdf",width=8,height=5,family="Arial")
big_plot
dev.off()
