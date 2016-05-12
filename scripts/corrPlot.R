library(ggplot2)
library(RColorBrewer)
library(data.table)
library(plyr)
source("~/Documents/zikaProj/doc/scripts/declareVariables.R")

countries <- c("pernambuco",
               "bahia",
               "saopaulo",
               "paraiba",
               "maranhao",
               "ceara",
               "sergipe",
               "riodejaneiro",
               "piaui",
               "riograndedonorte",
               "minasgerais",
               "matogrosso",
               "alagoas",
               "para",
               "acre",
               "goias",
               "espiritosanto",
               "tocantins")
country_names <- c(
  "pernambuco"="Pernambuco",
  "bahia"="Bahia",
  "saopaulo"="Sao Paulo",
  "paraiba"="Paraiba",
  "maranhao"="Maranhao",
  "ceara"="Ceara",
  "sergipe"="Sergipe",
  "riodejaneiro"="Rio de Janeiro",
  "piaui"="Piaui",
  "riograndedonorte"="Rio Grande Norte",
  "minasgerais"="Minas Gerais",
  "matogrosso"="Mato Grosso",
  "alagoas"="Algoas",
  "para"="Para",
  "acre"="Acre",
  "espiritosanto"="Espirito Santo",
  "goias"="Goias",
  "tocantins"="Tocantins"
)

allDat <- NULL
for(j in 1:length(countries)){
  place <- countries[j]
  chains <- NULL
  for(i in 1:3){
    filename <- paste(place,"_",i,"_chain.csv",sep="")
    greb <- fread(filename,data.table=FALSE)
    greb <- as.data.frame(greb[30000:nrow(greb),paramCols])
    chains <- rbind(chains, greb)
    
  }
  chains <- chains[sample(nrow(chains),1000),]
  chains <- as.data.frame(chains)
  chains$State <- country_names[as.character(place)]
  
  allDat <- rbind(allDat,chains)
}

getLevel <- function(x,y,prob=0.95) {
  kk <- MASS::kde2d(x,y)
  dx <- diff(kk$x[1:2])
  dy <- diff(kk$y[1:2])
  sz <- sort(kk$z)
  c1 <- cumsum(sz) * dx * dy
  approx(c1, sz, xout = 1 - prob)$y
}

p1 <- ggplot(allDat,aes(x=r0,y=epiStart,fill=State,colour=State)) + 
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+
  #geom_point(alpha=1/5,size=0.75) + 
  scale_alpha_continuous(range=c(0.1,0.5))+
  #stat_ellipse(aes(fill=State,colour=State),level=0.95)+
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  scale_x_continuous(limits=c(0,7.5),breaks=seq(1,7.5,by=0.5),expand=c(0,0))+
  scale_y_continuous(limits=c(-1,2),expand=c(0,0))+
  coord_cartesian(ylim=c(0,1.6),xlim=c(1,7.5))+
  xlab("R0")+
  ylab("epiStart")+
   overallTheme+ 
  mycolours +
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))

p2 <- ggplot(allDat,aes(x=r0,y=baselineProb,fill=State,col=State)) +
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+
  scale_alpha_continuous(range=c(0.1,0.5))+
#geom_point(alpha=1/5,size=0.75) + 
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  scale_x_continuous(limits=c(1,7.5),breaks=seq(1,7.5,by=0.5),expand=c(0,0))+
    xlab("R0")+
  ylab("bp")+
  #scale_x_continuous(limits=c(0,0.0005)) + 
  scale_y_continuous(limits=c(0,3e-4),expand=c(0,0),labels=seq(0,3e-4,by=1e-4),breaks=seq(0,3e-4,by=1e-4)) + 
  overallTheme+
  mycolours+
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))

p3 <- ggplot(allDat,aes(x=r0,y=probMicro,fill=State,col=State)) + 
  scale_y_continuous(limits=c(-1,2),expand=c(0,0))+
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+
  scale_alpha_continuous(range=c(0.1,0.5))+
 #geom_point(alpha=1/5,size=0.75) + 
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  scale_x_continuous(limits=c(1,7.5),breaks=seq(1,7.5,by=0.5),expand=c(0,0))+
  coord_cartesian(ylim=c(0,0.5),xlim=c(1,7.5))+
  #s
  xlab("R0")+
  ylab("pm")+
  #scale_y_log10()+
  #scale_x_continuous(limits=c(0,0.0005)) + 
  #scale_y_continuous(limits=c(0,1)) + 
  overallTheme+
  mycolours+
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))

p4 <- ggplot(allDat,aes(x=epiStart,y=baselineProb,fill=State,col=State)) + 
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+
  scale_alpha_continuous(range=c(0.1,0.5))+
 #geom_point(alpha=1/5,size=0.75) + 
 # guides(colour=FALSE)+
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  xlab("epiStart")+
  ylab("bp")+
  #scale_x_continuous(limits=c(0,0.0005)) + 
  scale_y_continuous(limits=c(0,3e-4),expand=c(0,0),labels=seq(0,3e-4,by=1e-4),breaks=seq(0,3e-4,by=1e-4)) + 
  scale_x_continuous(limits=c(0,1.6),expand=c(0,0))+
  overallTheme+
  mycolours+
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))

p5 <- ggplot(allDat,aes(x=epiStart,y=probMicro,fill=State,col=State)) + 
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+ 
  #scale_alpha_continuous(range=c(0.1,0.5))+
  #geom_point(alpha=1/5,size=0.75) + 
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  xlab("epiStart")+
  ylab("pm")+
  scale_y_continuous(limits=c(-1,1.5),expand=c(0,0))+
  scale_x_continuous(limits=c(-1,2),expand=c(0,0))+
  coord_cartesian(ylim=c(0,0.5),xlim=c(0,1.6))+
  #scale_x_continuous(limits=c(0,0.0005)) + 
  #scale_y_continuous(limits=c(0,1)) + 
  overallTheme+
  mycolours+
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))


p6 <- ggplot(allDat,aes(x=probMicro,y=baselineProb,fill=State,col=State)) + 
  stat_density2d(aes(alpha=..level..),geom="polygon",lwd=0,bins=5)+
  geom_density_2d(aes(colour=State),lwd=0.1,bins=5)+ 
  scale_alpha_continuous(range=c(0.1,0.5))+
  #geom_point(alpha=1/5,size=0.75) + 
 # guides(colour = guide_legend(override.aes= list(alpha = 1, size=2)))+
  scale_x_continuous(limits=c(-1,1.5),expand=c(0,0))+
  scale_y_continuous(limits=c(-1e-4,3e-4),expand=c(0,0),labels=seq(0,3e-4,by=1e-4),breaks=seq(0,3e-4,by=1e-4)) + 
  coord_cartesian(xlim=c(0,0.5),ylim=c(0,3e-4))+
  xlab("pm")+
  ylab("bp")+
  overallTheme+
  mycolours+
  mycolours2 +
  theme(plot.margin=unit(c(0.2,0.5,0.1,0.1),"cm"),text=element_text(size=18),legend.text=element_text(size=18),axis.title=element_text(size=18))
  
allP <- grid_arrange_shared_legend(p1,p2,p3,p4,p5,p6)

allCor <- melt(cor(allDat[,c("r0","probMicro","baselineProb","epiStart")],method="pearson"))
stateCor <- NULL
for(i in 1:length(unique(allDat$State))){
  tmp <- allDat[allDat$State==unique(allDat$State)[i],c("r0","probMicro","baselineProb","epiStart")]
  tmpCor <- melt(cor(tmp,method="pearson"))
  tmpCor <- tmpCor[c(5,9,13,10,14,15),]
  tmpCor$State <- unique(allDat$State)[i]
  stateCor <- rbind(stateCor,tmpCor)
}

r0base <- stateCor[stateCor$Var1=="r0" & stateCor$Var2=="baselineProb",c("value","State")]
r0base$pair <- "R0/bp"
r0prob<- stateCor[stateCor$Var1=="r0" & stateCor$Var2=="probMicro",c("value","State")]
r0prob$pair <- "R0/pm"
r0epi<- stateCor[stateCor$Var1=="r0" & stateCor$Var2=="epiStart",c("value","State")]
r0epi$pair <- "R0/start"
probbase<- stateCor[stateCor$Var1=="probMicro" & stateCor$Var2=="baselineProb",c("value","State")]
probbase$pair <- "pm/bp"
probepi<- stateCor[stateCor$Var1=="probMicro" & stateCor$Var2=="epiStart",c("value","State")]
probepi$pair <- "pm/start"
baseepi<- stateCor[stateCor$Var1=="baselineProb" & stateCor$Var2=="epiStart",c("value","State")]
baseepi$pair <- "bp/start"

pairNames <- c("R0.pm","R0.bp","R0.start","pm.bp","pm.start","bp.start")

allCorDat <- rbind(r0base,r0prob,r0epi,probbase,probepi,baseepi)
tmpMeans <- aggregate(allCorDat[,c("value")],list(allCorDat$pair),mean)
colnames(tmpMeans) <- c("pair","mean")

allCorDat <- join(allCorDat,tmpMeans,by="pair")

corComparisonPlot <- ggplot() + 
  geom_point(dat=allCorDat,aes(x=State,y=value,fill=pair,colour=pair),size=3) +
  scale_colour_brewer(palette="Set1")+
  geom_hline(dat=allCorDat,aes(yintercept=mean,colour=pair))+
  geom_hline(yintercept=0,linetype="dashed")+
  scale_y_continuous(limits=c(-1,1),expand=c(0,0))+
  ylab("Pearson Correlation Coefficient")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle=45,hjust=1)
  )

allP
corComparisonPlot

#tmpDat1 <- allDat[allDat$State=="Bahia",]
#y <- tmpDat1[,1]
#x <- tmpDat1[,4]
#g <- gam(y~s(x))
#plot(g)
