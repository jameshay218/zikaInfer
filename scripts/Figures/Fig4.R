######################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## This script produces Figure 4 from the manuscript.
## The user must have available the MCMC chains for
## each of the fits below (see scripts/job_submissions)
## Users should also provide the melted, combined MCMC chain
## see scripts/utility/combine_all_chains.R
library(zikaInfer)
library(ggplot2)
library(plyr)


meltedChain <- data.table::fread("~/Documents/Zika/results_20112017/all_chains.csv",data.table=FALSE,stringsAsFactors=FALSE)
meltedChain <- meltedChain[meltedChain$version == "model_1",]
chainwd <- "/media/james/JH USB/zika/outputs"


## Read in individual MCMC chains
bahia_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainwd, "bahia/bahia/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
colombia_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainwd, "colombia_inc/colombia/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
northeast_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainwd, "northeast/northeast/model_1"),asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
rio_chain <- zikaInfer::load_mcmc_chains(location = paste0(chainwd, "riograndedonorte/riograndedonorte/model_1"), asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)


## Format these chains for plotting
colombia_dat<- make_micro_dat(colombia_chain,10000)
colombia_dat$Model <- "Colombia"
bahia_dat<- make_micro_dat(bahia_chain,10000)
bahia_dat$Model <- "Bahia, Brazil"
nejm_dat<- make_micro_dat(northeast_chain, 10000)
nejm_dat$Model <- "Northeast Brazil"
rio_dat<- make_micro_dat(rio_chain,10000)
rio_dat$Model <- "Rio Grande do Norte, Brazil"

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
        legend.position="none",plot.title=element_text(hjust=0.5,size=14))


## Chain for trimester risk bounds
chain <- meltedChain
chain$runName <- as.character(chain$runName)
#chain <- chain[chain$runName %in% c("Rio Grande do Norte, Brazil (reports)", "Bahia, Brazil (reports)",
#                                    "Colombia","Northeast Brazil") 
chain <- chain[chain$runName %in% c("Rio Grande do Norte, Brazil (reports)", "Bahia, Brazil (reports)",
                                    "Colombia (suspected)","Colombia (confirmed)","Northeast Brazil, Lancet",
                                    "Pernambuco (confirmed, reports)", "Salvador, Bahia",
                                    "Rio Grande do Norte, Brazil (confirmed)") &
                 chain$variable %in% c("tr1","tr2","tr3"),]
trimesterRisks <- plyr::ddply(chain,c("runName", "variable"), function(x) quantile(x$value, c(0.025,0.5,0.975)))
colnames(trimesterRisks) <- c("Model","variable","low","mid","high")

## Regenerate microcephaly curve parameters from melted chain
meltedChain <- meltedChain[meltedChain$runName %in% c("Rio Grande do Norte, Brazil (reports)", "Bahia, Brazil (reports)",
                                                      "Colombia (suspected)","Colombia (confirmed)","Northeast Brazil, Lancet",
                                                      "Pernambuco (confirmed, reports)", "Salvador, Bahia",
                                                      "Rio Grande do Norte, Brazil (confirmed)"),]
microChain <- meltedChain[meltedChain$variable %in% c("mean","var","c") & meltedChain$version == "model_1",c("runName","variable","value","chain")]
allTmp <- NULL
for(ver in unique(microChain$runName)){
  allTmp[[ver]] <- NULL
  for(var in unique(microChain$variable)){
    tmp <- microChain[microChain$runName == ver & microChain$variable == var,]
    for(chain in unique(tmp$chain)){
      x <- 1:nrow(microChain[microChain$runName == ver & microChain$variable == var & microChain$chain == chain,])
      tmp1 <- microChain[microChain$runName == ver & microChain$variable == var & microChain$chain == chain,]
      tmp1 <- cbind(tmp1, sampno=x)
      allTmp[[ver]] <- rbind(allTmp[[ver]], tmp1)
    }
  }
  allTmp[[ver]]$seq <- 1:nrow(allTmp[[ver]])
}
dat <- NULL
for(i in 1:length(allTmp)){
  allTmp[[i]] <- reshape2::dcast(allTmp[[i]], runName + sampno+chain~variable,value.var="value")
  tmpDat <- make_micro_dat(allTmp[[i]],1000)
  tmpDat$Model <- unique(allTmp[[i]]$runName)
  dat <- rbind(dat,tmpDat)
}
dat$confirmed <- "notified"
dat[dat$Model %in% c("Colombia (confirmed)","Rio Grande do Norte, Brazil (confirmed)"),"confirmed"] <- "confirmed"
dat1 <- dat
dat1[dat1$Model == "Colombia (confirmed)","Model"] <- "Colombia"
dat1[dat1$Model == "Colombia (confirmed)","Model"] <- "Colombia"
dat1[dat1$Model == "Bahia, Brazil (reports)","Model"] <- "Bahia, Brazil"
dat1[dat1$Model == "Rio Grande do Norte, Brazil (confirmed)","Model"] <- "Rio Grande do Norte, Brazil"
dat1[dat1$Model == "Rio Grande do Norte, Brazil (reports)","Model"] <- "Rio Grande do Norte, Brazil"
dat1[dat1$Model == "Pernambuco (confirmed, reports)","Model"] <- "Pernambuco, Brazil (confirmed)*"

trimesterRisks$confirmed <- "notified"
trimesterRisks[trimesterRisks$Model %in% c("Colombia (confirmed)","Rio Grande do Norte, Brazil (confirmed)"),"confirmed"] <- "confirmed"
trimesterRisks1 <- trimesterRisks
trimesterRisks1[trimesterRisks1$Model == "Colombia (confirmed)","Model"] <- "Colombia"
trimesterRisks1[trimesterRisks1$Model == "Colombia (suspected)","Model"] <- "Colombia"
trimesterRisks1[trimesterRisks1$Model == "Bahia, Brazil (reports)","Model"] <- "Bahia, Brazil"
trimesterRisks1[trimesterRisks1$Model == "Rio Grande do Norte, Brazil (confirmed)","Model"] <- "Rio Grande do Norte, Brazil"
trimesterRisks1[trimesterRisks1$Model == "Rio Grande do Norte, Brazil (reports)","Model"] <- "Rio Grande do Norte, Brazil"
trimesterRisks1[trimesterRisks1$Model == "Pernambuco (confirmed, reports)","Model"] <- "Pernambuco, Brazil (confirmed)*"
trimesterRisks <- trimesterRisks1

dat1$Model <- as.character(dat1$Model)
dat1[dat1$Model == "Northeast Brazil, Lancet","Model"] <- "Northeast Brazil"
dat1[dat1$Model == "Colombia (suspected)","Model"] <- "Colombia"
dat1[dat1$Model == "Colombia (confirmed)","Model"] <- "Colombia"
dat1$Model <- as.factor(dat1$Model)
order <- c(3,2,1,5,6,4)
dat1$Model <- factor(dat1$Model,levels(dat1$Model)[order])

trimesterRisks$Model <- as.character(trimesterRisks$Model)
trimesterRisks[trimesterRisks$Model == "Northeast Brazil, Lancet","Model"] <- "Northeast Brazil"
trimesterRisks[trimesterRisks$Model == "Colombia (notified)","Model"] <- "Colombia"
trimesterRisks$Model <- as.factor(trimesterRisks$Model)
order <- c(3,2,1,5,6,4)
trimesterRisks$Model <- factor(trimesterRisks$Model,levels(trimesterRisks$Model)[order])

## Plot microcephaly risk curves

p <- ggplot(data=dat1) +
    geom_vline(xintercept=c(14*7, 28*7),col="grey40",lty="dashed",size=1) +
    geom_ribbon(aes(x=weeks,ymax=upper,ymin=lower,fill=Model, alpha=confirmed)) +
    geom_line(aes(x=weeks,y=means,group=confirmed, alpha=confirmed),col="black") + 
    geom_line(data=dat1[dat1$Model == "Pernambuco, Brazil (confirmed)*",],aes(x=weeks,y=means,group=confirmed,alpha=confirmed),col="purple",size=1.5) +

    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=0,xend=14*7,y=mid,yend=mid, alpha=confirmed),col="red",linetype="dashed",size=0.5) +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=low,yend=low, alpha=confirmed),col="red",size=0.5) +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=high,yend=high, alpha=confirmed),col="red",size=0.5) +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr1",],aes(x=(14*7 /2),ymin=low,ymax=high, alpha=confirmed),col="red") +
    
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7,xend=28*7,y=mid,yend=mid, alpha=confirmed),col="red",linetype="dashed") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=low,yend=low, alpha=confirmed),col="red",size=0.5) +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=high,yend=high, alpha=confirmed),col="red",size=0.5) +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr2",],aes(x=14*7 + (14*7 /2),ymin=low,ymax=high, alpha=confirmed),col="red",size=0.5) +
    
    
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7,xend=40*7,y=mid,yend=mid, alpha=confirmed),col="red",linetype="dashed") +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=low,yend=low, alpha=confirmed),col="red",size=0.5) +
    geom_segment(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=high,yend=high, alpha=confirmed),col="red",size=0.5) +
    geom_linerange(data=trimesterRisks[trimesterRisks$variable=="tr3",],aes(x=28*7 + (14*7 /2),ymin=low,ymax=high, alpha=confirmed),col="red",size=0.5) +
    scale_alpha_manual(values=c(0.2,0.8)) +
    facet_wrap(~Model,scales="free_y",ncol=2) +
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
png("Fig4.png",width=7,height=6,units="in",res=300)
p
dev.off()

cairo_ps("Fig4.eps",width=7,height=6,family="Arial",fallback_resolution = 600)
p
dev.off()

png("densities.png",width=7,height=4,units="in",res=300)
densities
dev.off()

cairo_ps("densities.eps",width=7,height=4,family="Arial")
densities
dev.off()
