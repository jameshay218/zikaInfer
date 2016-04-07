# Plots posteriors as volcano plots
library(ggplot2)
library(scales)
library(data.table)
library(gridExtra)
library(gtable)

setwd("~/tmpZika5")
burnin <-30000
paramCols <- c("probMicro","baselineProb","epiStart","r0")
actualDat <- read.csv("~/Dropbox/Zika/Data/allDat07.04.16.csv")
correctOrder <- sort(by(actualDat[,"microCeph"],actualDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
correctOrder <- as.factor(names(correctOrder))
files <- list.files(pattern="\\.csv$")


allR0 <- NULL
allEpi <- NULL
allProb <- NULL
allBaseline <- NULL
allAttack <- NULL

for(i in 1:length(correctOrder)){
    for(j in 1:3){
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

allR0$dat <- factor(allR0$dat,correctOrder)
allEpi$dat <- factor(allEpi$dat,correctOrder)
allProb$dat <- factor(allProb$dat,correctOrder)
allBaseline$dat <- factor(allBaseline$dat,correctOrder)
allAttack$dat <- factor(allAttack$dat,correctOrder)

R0plot <- ggplot(allR0,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
 # scale_x_continuous(limits=c(0,5))+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("R0")+
  xlab("Value")+
  theme(
      strip.text.x=element_text(size=8,angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(size=16,colour="gray20"),
      axis.line=element_line(colour="gray20"),
      axis.line.x = element_blank(),
      axis.line.y=element_line(colour="gray20"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
  )


attackplot <- ggplot(allAttack,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  # scale_x_continuous(limits=c(0,5))+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("Attack Rate")+
  xlab("Value")+
  theme(
    strip.text.x=element_text(size=8,angle=90),
    panel.grid.major = element_blank(),
    panel.grid.minor=element_blank(),
    text=element_text(size=16,colour="gray20"),
    axis.line=element_line(colour="gray20"),
    axis.line.x = element_blank(),
    axis.line.y=element_line(colour="gray20"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x=element_blank(),
    legend.position="none"
  )


epiplot <- ggplot(allEpi,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
 # scale_x_continuous(limits=c(0,5))+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("Epidemic Start Time")+
  xlab("Value")+
  theme(
      strip.text.x=element_text(size=8,angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(size=16,colour="gray20"),
      axis.line=element_line(colour="gray20"),
      axis.line.x = element_blank(),
      axis.line.y=element_line(colour="gray20"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
  )



baseplot <- ggplot(allBaseline,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_log10()+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("Baseline Probability of Microcephaly")+
  xlab("Value")+
  theme(
      strip.text.x=element_text(size=8,angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(size=16,colour="gray20"),
      axis.line=element_line(colour="gray20"),
      axis.line.x = element_blank(),
      axis.line.y=element_line(colour="gray20"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
  )



probplot <- ggplot(allProb,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
  scale_x_log10(breaks=c(0.00001,0.0001,0.001,0.01,0.1,0.2,0.5,1))+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("Probability of Microcephaly Given Infection")+
  xlab("Log Probability")+
  theme(
      strip.text.x=element_text(size=8,angle=90),
      panel.grid.major = element_blank(),
      panel.grid.minor=element_blank(),
      text=element_text(size=16,colour="gray20"),
      axis.line=element_line(colour="gray20"),
      axis.line.x = element_blank(),
      axis.line.y=element_line(colour="gray20"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x=element_blank(),
      legend.position="none"
  )

p <- grid.arrange(R0plot, epiplot,baseplot,probplot,ncol=1)
png("r0.png")
plot(R0plot)
dev.off()

png("baseline.png")
plot(baseplot)
dev.off()

png("probMicro.png")
plot(probplot)
dev.off()

png("epiStart.png")
plot(epiplot)
dev.off()

png("attackRate.png")
plot(attackplot)
dev.off()

plots <- list(attackplot, R0plot, epiplot, probplot)
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

filenames <- c("attackPlot.png","R0plot.png","epiPlot.png","probPlot.png")
index <- 1
for(p1 in newPlots){
  png(filenames[index])
  plot(p1)
  dev.off()
  index <- index + 1
}
