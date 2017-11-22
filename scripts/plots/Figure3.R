######################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This script produces Figure 3 from the manuscript.
## The user must have available the MCMC chains for
## each of the fits below (see scripts/job_submissions)
library(cowplot)
library(zoo)
library(ggplot2)
library(coda)
library(extrafont)
library(hexbin)
library(RColorBrewer)
library(zikaProj)

chainWD <- "~/net/home/zika/outputs/"
parTab <- read.csv("~/net/home/zika/inputs/parTab_fixed_switch_early.csv",stringsAsFactors=FALSE)


heatmap_theme <- theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
                       plot.title=element_text(size=8,family="Arial",color="black"),
                       axis.text.y=element_text(size=8,family="Arial",color="black"),
                       axis.title.x=element_text(size=10,family="Arial",color="black"),
                       axis.title.y=element_text(size=10,family="Arial",color="black"),
                       legend.text=element_text(size=8,family="Arial",color="black"),
                       legend.title=element_text(size=8,family="Arial",color="black"),
                       strip.text=element_text(size=8,family="Arial")) + theme_bw()

# 2d densities ------------------------------------------------------------
bahia_chain <- as.data.frame(lazymcmc::load_mcmc_chains(paste0(chainWD,"bahia_forecast"),parTab,FALSE,100,50000,
                                          TRUE,FALSE,FALSE)[["chain"]])
bahia_chain$propn_increase <- bahia_chain$incPropn2/bahia_chain$incPropn

bahia_chain_noreportingchange <- as.data.frame(lazymcmc::load_mcmc_chains(
                                                             paste0(chainWD,"bahia_forecast_noreportingchange"),
                                                             parTab,FALSE,100,50000,
                                                             TRUE,FALSE,FALSE)[["chain"]])

datFile = "~/Documents/Zika/Data/brazil/microceph_reports_2016.csv"
incFile = "~/Documents/Zika/Data/brazil//zika_inc_reports.csv"
local <- "bahia"
microDat <- read.csv(datFile,
                     stringsAsFactors=FALSE)
incDat <- read.csv(incFile,
                   stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == local,]
incDat <- incDat[incDat$local == local,]

parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365
f <- create_forecast_function(parTab, microDat, incDat=incDat, ts=seq(0,3003,by=1), FALSE)
aborted <- NULL

for(i in 1:nrow(bahia_chain_noreportingchange)){
    pars <- get_index_pars(bahia_chain_noreportingchange,i)
    pars["baselineProb"] <- exp(pars["baselineProb"])
    aborted[i] <- sum(f(pars,TRUE)$aborted$aborted)
}
bahia_chain_noreportingchange$abortions <- aborted


heatmap_theme <- heatmap_theme + theme(
  axis.text.x=element_text(size=8,family="Arial",color="black"),
  axis.text.y=element_text(size=8,family="Arial",color="black"),
  legend.text=element_text(size=6,family="Arial",color="black")
)

p1 <- ggplot() + geom_hex(data=bahia_chain,aes(x=birth_reduction,y=propn_increase,fill=..density..))+ 
    scale_fill_gradient2(low="darkturquoise",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
ylab("Relative increase in\n ZIKV reporting rate") +
    xlab("Proportion decrease in\nZIKV affected births") +
    labs(fill="Density") + heatmap_theme
p2 <- ggplot() + geom_hex(data=bahia_chain,aes(x=propn_increase,y=propn,fill=..density..))+ 
    scale_fill_gradient2(low="darkturquoise",mid="#FAFDB8",high="#9E0142",midpoint= 0.005) +
ylab("Microcephaly reporting\nrate in first wave") +
    xlab("Relative increase in\n ZIKV reporting rate") +
    labs(fill="Density") + heatmap_theme
p3 <- ggplot() + geom_hex(data=bahia_chain_noreportingchange,aes(x=birth_reduction,y=abortion_rate,fill=..density..),bins=30)+ 
    scale_fill_gradient2(low="turquoise4",mid="#FAFDB8",high="#9E0142",midpoint= 0.0075) +
scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=0.25)) +
    ylab("Proportion of microcephaly\naffected births aborted") +
    xlab("Relative reduction in infection\nrisk in pregnant women") +
    labs(fill="Density") + heatmap_theme
p4 <- ggplot() + geom_hex(data=bahia_chain_noreportingchange,aes(x=abortion_rate,y=abortions,fill=..density..),bins=200)+ 
    scale_fill_gradient2(low="turquoise4",mid="#FAFDB8",high="#9E0142",midpoint= 0.002) +
ylab("Total number of microcephaly\naffected births aborted") +
    xlab("Proportion of microcephaly\naffected births aborted") +
    labs(fill="Density") + heatmap_theme

p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
p3 <- p3 + theme(legend.position="none")
p4 <- p4 + theme(legend.position="none")



                                        # Main plot ---------------------------------------------------------------
fig3 <- plot_grid(p1,p2,p3,p4,ncol=2)

png("Fig3.png",width=7,height=6,units="in",res=300)
fig3
dev.off()

cairo_ps("Fig3.eps",width=7,height=6,family="Arial")
fig3
dev.off()


find_mode <- function(x){
  y <- density(x)
  return(y$x[which.max(y$y)])
}

chain <- bahia_chain_noreportingchange
for(par in c("abortion_rate","birth_reduction","propn","propn_increase","abortions")){
  print(paste0(par, " mode: ", find_mode(chain[,par])))
  print(paste0(par, " quantiles: ", quantile(chain[,par],c(0.025,0.5,0.975))))
}
