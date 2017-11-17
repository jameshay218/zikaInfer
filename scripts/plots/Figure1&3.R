library(zikaProj)
library(ggplot2)
library(plyr)

make_micro_dat <- function(chain, samp_no,scale=FALSE){
  samples <- sample(nrow(chain), samp_no)
  microCurves <- matrix(nrow=samp_no,ncol=40*7)
  trimesterCurves <- matrix(nrow=samp_no,ncol=40*7)
  tm1 <- tm2 <- tm3 <- numeric(samp_no)
  i <- 1
  for(samp in samples){
    pars <- chain[samp,]
    pars <- as.numeric(chain[samp,])
    names(pars) <- colnames(chain)
    probs <- generate_micro_curve(pars)
    probs <- probs
    tm1[i] <- mean(probs[1:14*7])
    tm2[i] <- mean(probs[(14*7 + 1):(28*7)])
    tm3[i] <- mean(probs[(28*7 + 1):length(probs)])
    microCurves[i,] <- probs
    i <- i + 1
  }
  #microCurves <- microCurves/max(microCurves)
   means <- colMeans(microCurves)
  lower <- apply(microCurves,2,function(x) quantile(x, 0.025))
  upper <- apply(microCurves,2,function(x) quantile(x,0.975))
  if(scale){
    means <- means/max(upper)
    lower <- lower/max(upper)
    upper <- upper/max(upper)
  }
  dat <- data.frame(weeks=1:(40*7),means,upper,lower)
  return(dat)
}

brazil_chain <- zikaProj::load_mcmc_chains(location = "~/Documents/Zika/28.02.2017_chains/multi_all1/model_1",asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
colombia_chain <- zikaProj::load_mcmc_chains(location = "~/Documents/Zika/28.02.2017_chains/colombia_inc/bahia/model_1", asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
northeast_chain <- zikaProj::load_mcmc_chains(location = "~/Documents/Zika/28.02.2017_chains/northeast/model_1/",asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)
rio_chain <- zikaProj::load_mcmc_chains(location = "~/Documents/Zika/28.02.2017_chains/riograndedonorte/riograndedonorte/model_1", asList = FALSE,convertMCMC = FALSE,unfixed = FALSE,thin = 10,burnin = 750000)

colombia_dat<- make_micro_dat(colombia_chain,10000)
colombia_dat$Model <- "Colombia"
bahia_dat<- make_micro_dat(bahia_chain,10000)
bahia_dat$Model <- "Bahia, Brazil - reports"
nejm_dat<- make_micro_dat(northeast_chain, 10000)
nejm_dat$Model <- "Northeast Brazil"
rio_dat<- make_micro_dat(rio_chain,10000)
rio_dat$Model <- "Rio Grande do Norte, Brazil - reports"

dat <- rbind(colombia_dat, bahia_dat, nejm_dat, rio_dat)
dat$Model <- as.factor(dat$Model)
order <- c(1,2,3,4)
dat$Model <- factor(dat$Model,levels(dat$Model)[order])

p <- ggplot(data=dat) +
  geom_vline(xintercept=c(14*7, 28*7),col="grey40",lty="dashed",size=1) +
  geom_ribbon(aes(x=weeks,ymax=upper,ymin=lower,fill=Model,group=Model),alpha=0.6) +
    geom_line(aes(x=weeks,y=means,group=Model),col="black",lty="dashed") +
  geom_segment(data=wow[wow$variable=="tr1",],aes(x=0,xend=14*7,y=mid,yend=mid)) +
    facet_wrap(~Model,scales="free_y") +
  ylab("Probability of developing\nmicrocephaly given ZIKV infection")+
  xlab("Week of gestation at time of infection") +
  theme_bw() +
  #scale_fill_manual(values = c("#F8766D","#C49A00","#00B6EB","#FB61D7","blue") ) +
  theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
        #plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
        panel.grid.minor=element_blank(),
        plot.title=element_text(size=8,family="Arial",color="black"),
        axis.text.y=element_text(size=8,family="Arial",color="black"),
        axis.title.x=element_text(size=10,family="Arial",color="black"),
        axis.title.y=element_text(size=10,family="Arial",color="black"),
        #legend.position=c(0.99,0.99),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text=element_text(size=8,family="Arial",color="black"),
        legend.justification=c(0,0)) +
  #scale_y_continuous(expand=c(0,0),limits=c(0,1)) +
  scale_x_continuous(expand=c(0,0),limits=c(0,40*7),breaks=seq(0,280,by=7*4), labels=seq(0,40,by=4))

cairo_ps("~/Documents/Zika/plots/manuscript/Fig3_riskcurves.eps",width=8,height=6,family="Arial")
p
dev.off()



chain <- data.table::fread("~/Documents/Zika/test/combinedChains2.csv",data.table=FALSE,stringsAsFactors=FALSE)
runNameConvert <- c("multi_all1" = "Brazil Model",
                    "colombia_inc/bahia" = "Colombia",
                    "bahia/bahia" = "Bahia - Reports",
                    "riograndedonorte/riograndedonorte" = "Rio Grande do Norte - Reports",
                    "pernambuco_peak/pernambuco" = "Pernambuco - Reports",
                    "reports_3_inc"="Bahia + Pernambuco\n+ Rio Grande do Norte",
                    "reports_2_inc"="Bahia + Rio Grande do Norte") 
chain$runName <- as.character(chain$runName)
chain$runName <- runNameConvert[chain$runName]
chain$runName <- factor(chain$runName)
chain$runName <- factor(chain$runName, levels(chain$runName)[c(4,5,1,3,2,6,7)])

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

cowplot::plot_grid(p, densities,rel_widths = c(10,5),rel_heights = c(5,10))


chain <- data.table::fread("~/Documents/Zika/test/combinedChains_m1_final.csv",data.table=FALSE)
chain <- chain[chain$runName %in% c("riograndedonorte/riograndedonorte", "bahia/bahia",
                                    "colombia_inc/bahia","northeast") &
                 chain$variable %in% c("tr1","tr2","tr3"),]
change_names <- c("riograndedonorte/riograndedonorte"="Rio Grande do Norte, Brazil - reports",
                  "bahia/bahia" = "Bahia, Brazil - reports",
                  "colombia_inc/bahia"="Colombia",
                  "northeast"="Northeast Brazil")
wow <- ddply(chain,c("runName", "variable"), function(x) quantile(x$value, c(0.025,0.5,0.975)))
wow$runName <- change_names[wow$runName]
colnames(wow) <- c("Model","variable","low","mid","high")
p <- ggplot(data=dat) +
    geom_vline(xintercept=c(14*7, 28*7),col="grey40",lty="dashed",size=1) +
    geom_ribbon(aes(x=weeks,ymax=upper,ymin=lower,fill=Model,group=Model),alpha=0.6) +
    geom_line(aes(x=weeks,y=means,group=Model),col="black") +
    
    geom_segment(data=wow[wow$variable=="tr1",],aes(x=0,xend=14*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=wow[wow$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=wow[wow$variable=="tr1",],aes(x=(14*7 /3),xend=(14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=wow[wow$variable=="tr1",],aes(x=(14*7 /2),ymin=low,ymax=high),col="red") +
    
    geom_segment(data=wow[wow$variable=="tr2",],aes(x=14*7,xend=28*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=wow[wow$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=wow[wow$variable=="tr2",],aes(x=14*7 + (14*7 /3),xend=14*7 + (14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=wow[wow$variable=="tr2",],aes(x=14*7 + (14*7 /2),ymin=low,ymax=high),col="red") +
    
    
    geom_segment(data=wow[wow$variable=="tr3",],aes(x=28*7,xend=40*7,y=mid,yend=mid),col="red",linetype="dashed") +
    geom_segment(data=wow[wow$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=low,yend=low),col="red") +
    geom_segment(data=wow[wow$variable=="tr3",],aes(x=28*7 + (14*7 /3),xend=28*7 + (14*7*2 /3),y=high,yend=high),col="red") +
    geom_linerange(data=wow[wow$variable=="tr3",],aes(x=28*7 + (14*7 /2),ymin=low,ymax=high),col="red") +
    
    
    facet_wrap(~Model,scales="free_y") +
    ylab("Probability of developing\nmicrocephaly given ZIKV infection")+
    xlab("Gestational age (in weeks) at time of infection") +
    theme_bw() +
                                        #scale_fill_manual(values = c("#F8766D","#C49A00","#00B6EB","#FB61D7","blue") ) +
    theme(axis.text.x=element_text(size=8,family="Arial",color="black"),
                                        #plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
          panel.grid.minor=element_blank(),
          plot.title=element_text(size=8,family="Arial",color="black"),
          axis.text.y=element_text(size=8,family="Arial",color="black"),
          axis.title.x=element_text(size=10,family="Arial",color="black"),
          axis.title.y=element_text(size=10,family="Arial",color="black"),
                                        #legend.position=c(0.99,0.99),
          legend.position="none",
          legend.title=element_blank(),
          legend.text.align = 0.5,
          legend.text=element_text(size=8,family="Arial",color="black"),
          legend.direction="horizontal") +
                                        #scale_y_continuous(expand=c(0,0),limits=c(0,1)) +
    scale_x_continuous(expand=c(0,0),limits=c(0,40*7),breaks=seq(0,280,by=7*4), labels=seq(0,40,by=4))
png("~/Documents/Zika/plots/manuscript/to_use/Fig1_essay.png",width=7,height=4,units="in",res=300)
p
dev.off()

cairo_ps("~/Documents/Zika/plots/manuscript/to_use/Fig1_essay.eps",width=7,height=4,family="Arial")
p
dev.off()

cairo_ps("~/Documents/Zika/plots/manuscript/to_use/Fig1_essay.pdf",width=7,height=4,family="Arial")
p
dev.off()
