# Plots posteriors as volcano plots
library(ggplot2)
library(scales)
library(data.table)
library(gridExtra)
library(grid)
library(gtable)
library(scales)
library(nleqslv)

setwd("~/tmpZika6")
burnin <-30000
paramCols <- c("probMicro","baselineProb","epiStart","r0")
actualDat <- read.csv("~/Dropbox/Zika/Data/allDat07.04.16.csv")
#correctOrder <- sort(by(actualDat[,"microCeph"],actualDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
#correctOrder <- as.factor(names(correctOrder))
#correctOrder <- places
places <- unique(actualDat[,"local"])
correctOrder <- sort(places)
files <- list.files(pattern="\\.csv$")


allR0 <- NULL
allEpi <- NULL
allProb <- NULL
allBaseline <- NULL
allAttack <- NULL

#skip <- c("bahia_2","espiritosanto_2","goias_1","para_1","riodejaneiro_3","saopaulo_2")

for(i in 1:length(correctOrder)){
    for(j in 1:3){
      tmp <- paste(correctOrder[i],"_",j,sep="")
     # if(tmp %in% skip) next
      
      filename <- paste(correctOrder[i],"_",j,"_chain.csv",sep="")
      print(filename)
      tmp <- fread(filename,data.table=FALSE)
      tmp <- as.data.frame(tmp[burnin:nrow(tmp),paramCols])

      tmp <- tmp[sample(seq(1,nrow(tmp),by=1),1000,replace=FALSE),]
      
      name <- as.character(correctOrder[i])
      print(name)
    
      tmpR0 <- data.frame("value"=tmp$r0,"dat"=name,"chain"=j)
      tmpAttack <- data.frame("value"=sapply(tmp$r0, function(x) nleqslv(0.8, zikaProj::simeq, R0=x)$x),"dat"=name,"chain"=j)
      tmpEpi <- data.frame("value"=tmp$epiStart,"dat"=name,"chain"=j)
      tmpProb <- data.frame("value"=tmp$probMicro,"dat"=name,"chain"=j)
      tmpBase <- data.frame("value"=tmp$baselineProb,"dat"=name,"chain"=j)

      allR0 <- rbind(allR0, tmpR0)
      allAttack <- rbind(allAttack, tmpAttack)
      allEpi <- rbind(allEpi, tmpEpi)
      allProb <- rbind(allProb, tmpProb)
      allBaseline <- rbind(allBaseline, tmpBase)
    }
}

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


allR0$dat <- factor(allR0$dat,correctOrder)
allR0 <- allR0[allR0$dat %in% countries,]

allEpi$dat <- factor(allEpi$dat,correctOrder)
allEpi <- allEpi[allEpi$dat %in% countries,]

allProb$dat <- factor(allProb$dat,correctOrder)
allProb <- allProb[allProb$dat %in% countries,]

allBaseline$dat <- factor(allBaseline$dat,correctOrder)
allBaseline <- allBaseline[allBaseline$dat %in% countries,]

allAttack$dat <- factor(allAttack$dat,correctOrder)
allAttack <- allAttack[allAttack$dat %in% countries,]

oldTheme <- theme(
  strip.text.x=element_text(size=8,angle=90),
  #panel.grid.major = element_blank(),
  panel.grid.minor=element_blank(),
  panel.grid.major.x = element_blank(),
  text=element_text(size=12,colour="gray20"),
  axis.line=element_line(colour="gray20"),
  axis.line.x = element_blank(),
  axis.line.y=element_line(colour="gray20"),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.x=element_blank(),
  legend.position="none"
)

mytheme <- theme(
  #panel.grid.major = element_blank(),
  panel.grid.minor=element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(colour="gray20",linetype="dashed",size=0.25),
  text=element_text(size=16,colour="gray20"),
  axis.line=element_line(colour="gray20"),
  #axis.line.x = element_blank(),
  axis.line.y=element_line(colour="gray20"),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  #axis.title.y =element_blank(),
  strip.text.y=element_text(size=18,angle=0),
  strip.text = element_blank(),
  legend.position="none",
  panel.margin=unit(0,"lines"),
  plot.margin=unit(c(0.1,0.5,0.1,0.1),"cm"),
  panel.grid = element_line(colour="gray20"),
  panel.background=element_rect(colour="gray40"),
  strip.background=element_rect(colour="black")
)

R0plot <- ggplot(allR0,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_continuous(limits=c(0,5),breaks=seq(0,5,by=0.5))+
  facet_grid(dat~.,labeller=as_labeller(country_names))+
  #coord_flip()+ 
  ylab("") +
  ggtitle("Basic Reproductive Number, R0")+
  xlab("Value")+
  mytheme +
  theme(
    plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
  )



attackplot <- ggplot(allAttack,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_continuous(limits=c(0.5,1),breaks=seq(0.5,1,by=0.1),labels=c("50%","60%","70%","80%","90%","100%"))+
  facet_grid(dat~.)+
  #coord_flip()+ 
  ylab("") +
  ggtitle("Attack Rate")+
  xlab("Value")+
  mytheme +
  theme(
    plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
  )+
  coord_cartesian(ylim=c(-50,50))



epiplot <- ggplot(allEpi,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_continuous(limits=c(0,2),breaks=seq(0,2,by=3/12),labels=c("07/2013","10/2013","01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015"))+
  facet_grid(dat~.)+
  ylab("") +
  ggtitle("Epidemic Start Time")+
  xlab("Value")+
  mytheme



baseplot <- ggplot(allBaseline,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_log10(limits=c(1e-5,1e-3),breaks=c(1e-5,1e-4,1e-3),labels=c("0.001%","0.01%","0.1%"))+
  facet_grid(dat~.)+
  ylab("") +
  ggtitle("Baseline Probability of Microcephaly")+
  xlab("Value (log scale)")+
  mytheme+
  theme(
    plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
  )



probplot <- ggplot(allProb,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_log10(limits=c(0.01,0.5),breaks=c(0.01,0.05,0.1,0.2,0.5),labels=c("1%","5%","10%","20%","50%"))+
  facet_grid(dat~.)+
  ylab("") +
  ggtitle("Probability of Microcephaly Given Infection")+
  xlab("Value (log scale)")+
  mytheme+
  theme(
    plot.margin=unit(c(0.1,0.5,0.1,0),"cm")
  )  
  
strips <- R0plot + theme(strip.text=element_text())
tmp <- ggplot_gtable(ggplot_build(strips))
tmp1 <- gtable_filter(tmp,"strip-right")

# g1 <- grid.arrange(arrangeGrob(textGrob(label="State",gp=gpar(fontsize=14,colour="grey20")),
#                                tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.35,9,0.53)),
#                    R0plot,attackplot,ncol=3,widths=c(3,10,10))
# g2 <- grid.arrange(arrangeGrob(textGrob(label="State",gp=gpar(fontsize=14,colour="grey20")),
#                               tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.35,9,0.53)),
#                   probplot,baseplot,ncol=3,widths=c(3,10,10))

g3 <- grid.arrange(arrangeGrob(textGrob(label="State",gp=gpar(fontsize=16,colour="grey20")),
                               tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.55,7,0.8)),
                   R0plot,attackplot,arrangeGrob(textGrob(label="State",gp=gpar(fontsize=16,colour="grey20")),
                                                 tmp1,rectGrob(gp=gpar(col="white")),heights=c(0.55,7,0.8)),
                   probplot,baseplot,ncol=3,widths=c(3,10,10))

p <- grid.arrange(R0plot, epiplot,baseplot,probplot,ncol=1)

# png("r0.png")
# plot(R0plot)
# dev.off()
# 
# png("baseline.png")
# plot(baseplot)
# dev.off()
# 
# png("probMicro.png")
# plot(probplot)
# dev.off()
# 
# png("epiStart.png")
# plot(epiplot)
# dev.off()
# 
# png("attackRate.png")
# plot(attackplot)
# dev.off()

plots <- list(attackplot, R0plot, epiplot, probplot,baseplot)
newPlots <- NULL
index <- 1
for(p1 in plots){
  # Convert the plot to a grob
  gt <- ggplotGrob(p1)
  
  # Get the positions of the panels in the layout: t = top, l = left, ...
  panels <-c(subset(gt$layout, name == "panel", select = t:r))

  # Add a row below the x-axis tick mark labels,
  # the same height as the strip
  gt = gtable_add_rows(gt, gt$height[min(panels$t)-1], max(panels$b) + 2)
 
  # Get the strip grob
  stripGrob = gtable_filter(gt, "strip-top")

  
  # Insert the strip grob into the new row
  gt = gtable_add_grob(gt, stripGrob, t = max(panels$b) + 3, l = min(panels$l), r = max(panels$r))

  # remove the old strip
  
 
  gt = gt[-(min(panels$t)-1), ]

  newPlots[[index]] <- gt
  index <- index + 1
}
all <- grid.arrange(newPlots[[2]],newPlots[[4]],newPlots[[3]],newPlots[[5]])

filenames <- c("attackPlot.svg","R0plot.svg","epiPlot.svg","probPlot.svg","baselinePlot.svg")
index <- 1
library(RSvgDevice)
# for(p1 in newPlots){
#   devSVG(filenames[index])
#   plot(p1)
#   dev.off()
#   index <- index + 1
# }
