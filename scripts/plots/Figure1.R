######################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This script produces Figure 1 from the manuscript.
## The user must have available the MCMC chains for
## each of the fits below (see scripts/job_submissions)
## Users should also provide the melted, combined MCMC chain
## see scripts/utility/combine_all_chains.R
library(zikaProj)
library(ggplot2)
library(plyr)


meltedChain <- data.table::fread("~/Documents/Zika/results_20112017/combinedChains.csv",data.table=FALSE,stringsAsFactors=FALSE)
meltedChain <- meltedChain[meltedChain$version == "model_1",]
chainwd <- "~/net/home/zika/outputs/"


## Read in individual MCMC chains
bahia_chain <- zikaProj::load_mcmc_chains(location = paste0(chainwd, "bahia/bahia/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
colombia_chain <- zikaProj::load_mcmc_chains(location = paste0(chainwd, "colombia_inc/colombia/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
northeast_chain <- zikaProj::load_mcmc_chains(location = paste0(chainwd, "northeast/northeast/model_1"),asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
rio_chain <- zikaProj::load_mcmc_chains(location = paste0(chainwd, "riograndedonorte/riograndedonorte/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)


## Format these chains for plotting
colombia_dat<- make_micro_dat(colombia_chain,10000)
colombia_dat$Model <- "Colombia"
bahia_dat<- make_micro_dat(bahia_chain,10000)
bahia_dat$Model <- "Bahia, Brazil (reports)"
nejm_dat<- make_micro_dat(northeast_chain, 10000)
nejm_dat$Model <- "Northeast Brazil"
rio_dat<- make_micro_dat(rio_chain,10000)
rio_dat$Model <- "Rio Grande do Norte, Brazil (reports)"

dat <- rbind(colombia_dat, bahia_dat, nejm_dat, rio_dat)
dat$Model <- as.factor(dat$Model)
order <- c(1,2,3,4)
dat$Model <- factor(dat$Model,levels(dat$Model)[order])

## Plot posterior densities
chain <- meltedChain
chain[chain$variable %in% c("lower","upper","range","maxWeek"),"value"] <- chain[chain$variable %in% c("lower","upper","range","maxWeek"),"value"]/7
chain <- chain[chain$variable %in% c("lower","upper","range","maxWeek"),]
varNames <- c("lower"="First risk week",
              "upper"="Last risk week",
              "maxWeek"="Peak risk week",
              "range"="Number of risk weeks")
chain$variable <- varNames[chain$variable]
chain <- chain[chain$version %in% "model_1",]

densities <- ggplot(chain) +
  geom_violin(aes(x=runName,y=value,fill=runName),draw_quantiles=c(0.5),scale="width") +
  facet_wrap(~variable,ncol=1) +
  coord_cartesian(ylim=c(0,40)) + theme_bw() + 
  ylab("")+
  geom_hline(yintercept=c(14,28),linetype="dashed",col="red",alpha=0.5)+
  theme(legend.text=element_text(size=14),strip.text=element_text(size=12),axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=12,angle=45,hjust=1),
        axis.text.y=element_text(size=14),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),panel.grid.minor=element_blank(),
        legend.position="none",,plot.title=element_text(hjust=0.5,size=14))


## Chain for trimester risk bounds
chain <- meltedChain
chain$runName <- as.character(chain$runName)
chain <- chain[chain$runName %in% c("Rio Grande do Norte, Brazil (reports)", "Bahia, Brazil (reports)",
                                    "Colombia","Northeast Brazil") &
                 chain$variable %in% c("tr1","tr2","tr3"),]
trimesterRisks <- plyr::ddply(chain,c("runName", "variable"), function(x) quantile(x$value, c(0.025,0.5,0.975)))
colnames(trimesterRisks) <- c("Model","variable","low","mid","high")

## Plot microcephaly risk curves
p <- ggplot(data=dat) +
    geom_vline(xintercept=c(14*7, 28*7),col="grey40",lty="dashed",size=1) +
    geom_ribbon(aes(x=weeks,ymax=upper,ymin=lower,fill=Model,group=Model),alpha=0.6) +
    geom_line(aes(x=weeks,y=means,group=Model),col="black") +
    
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=0,xend=14*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /2),ymin=low,ymax=high),col="red") +
    
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7,xend=28*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /2),ymin=low,ymax=high),col="red") +
    
    
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7,xend=40*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /2),ymin=low,ymax=high),col="red") +

    facet_wrap(~Model,scales="free_y") +
    ylab("Probability of developing\nmicrocephaly given ZIKV infection")+
    xlab("Gestational age (in weeks) at time of infection") +
    theme_bw() +
    theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=8,family="Arial",color="black"),
          axis.text.y=element_text(size=8,family="Arial",color="black"),
          axis.title.x=element_text(size=10,family="Arial",color="black"),
          axis.title.y=element_text(size=10,family="Arial",color="black"),
          legend.position="none",
          legend.title=element_blank(),
          legend.text.align = 0.5,
          legend.text=element_text(size=8,family="Arial",color="black"),
          legend.direction="horizontal") +
    scale_x_continuous(expand=c(0,0),limits=c(0,40*7),breaks=seq(0,280,by=7*4), labels=seq(0,40,by=4))


## Save plots
png("Fig1.png",width=7,height=4,units="in",res=300)
p
dev.off()

cairo_ps("Fig1.eps",width=7,height=4,family="Arial")
p
dev.off()

png("densities.png",width=7,height=4,units="in",res=300)
densities
dev.off()

cairo_ps("densities.eps",width=7,height=4,family="Arial")
densities
dev.off()
