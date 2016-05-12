# Plots posteriors as volcano plots
library(ggplot2)
library(scales)
library(data.table)
trueMu <- 28
trueB <- 50
truesd <- 1.2
burnin <- 10000

allDatb <- NULL
allDatmu <- NULL
allDatsd <- NULL
filenames <- c("high_dat_chain.csv","half_dat_chain.csv","low_dat_chain.csv","lowish_dat_chain.csv","verylow_dat_chain.csv")
#filenames <- c("high_count_chain.csv","half_count_chain.csv","low_count_chain.csv","lowish_count_chain.csv","verylow_count_chain.csv")
samples <- c("90%","50%","10%","5%","1%")
i <- 1
for(file in filenames){
  tmp <- fread(file,data.table=FALSE)
  tmp <- tmp[10000:nrow(tmp),]
  tmp <- as.data.frame(tmp[,c(4,5,17)])
  tmpmu <- data.frame("value"=tmp$muI,"dat"=samples[i])
  allDatmu <- rbind(allDatmu,tmpmu)
  tmpsd <- data.frame("value"=tmp$sdI,"dat"=samples[i])
  allDatsd <- rbind(allDatsd,tmpsd)
  tmpb <- data.frame("value"=tmp$b,"dat"=samples[i])
  allDatb <- rbind(allDatb,tmpb)
  i <- i + 1
}
allDatmu <- as.data.frame(allDatmu[sample(seq(1,nrow(allDatmu),by=1),1000,replace=FALSE),])
colnames(allDatmu) <- c("value","dat")
allDatb <- as.data.frame(allDatb[sample(seq(1,nrow(allDatb),by=1),1000,replace=FALSE),])
colnames(allDatb) <- c("value","dat")
allDatsd <- as.data.frame(allDatsd[sample(seq(1,nrow(allDatsd),by=1),1000,replace=FALSE),])
colnames(allDatsd) <- c("value","dat")

pB1 <- ggplot(allDatb,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_continuous(limits=c(30,60))+
  facet_grid(.~dat)+
  coord_flip()+ 
  ylab("") +
  ggtitle("Bite Rate, b")+
  xlab("")+
  geom_vline(aes(xintercept=trueB),colour="black",linetype="longdash")+
  theme(
    panel.grid.major = element_blank(),
    text=element_text(size=16,colour="gray20"),
    axis.line=element_line(colour="gray20"),
    axis.line.x = element_line(colour = "gray20"),
    axis.line.y=element_line(colour="gray20"),
    legend.position="none"
  )
pmu1 <- ggplot(allDatmu,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  facet_grid(.~dat)+
  scale_y_continuous(breaks=pretty_breaks(n=3))+
  coord_flip()+ 
  ggtitle(expression(mu[I]))+
  ylab("") +
  xlab("Value")+
  geom_vline(aes(xintercept=trueMu),colour="black",linetype="longdash")+
  theme(
    panel.grid.major = element_blank(),
    text=element_text(size=16,colour="gray20"),
    axis.line=element_line(colour="gray20"),
    axis.line.x = element_line(colour = "gray20"),
    axis.line.y=element_line(colour="gray20"),
    legend.position="none"
  )

psd1 <- ggplot(allDatsd,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
 facet_grid(.~dat)+
  scale_x_continuous(limits=c(0,10))+
  scale_y_continuous(breaks=pretty_breaks(n=3))+
  coord_flip()+ 
  ylab("Density") +
  ggtitle(expression(sigma[I])) +
  xlab("")+
  geom_vline(aes(xintercept=truesd),colour="black",linetype="longdash")+
  theme(
    panel.grid.major = element_blank(),
    text=element_text(size=16,colour="gray20"),
    axis.line=element_line(colour="gray20"),
    axis.line.x = element_line(colour = "gray20"),
    axis.line.y=element_line(colour="gray20"),
    legend.position="none"
  )


plots <- grid.arrange(pB1,pmu1,psd1,nrow=3)

ggsave("posterior_dat.png",plots,height=9,width=12)

#
#b1<-ggplot(allDatb, aes(dat, value)) + 
#  geom_boxplot(aes(fill = dat)) +
#  theme(legend.position = "none") + geom_hline(aes(yintercept=trueB))