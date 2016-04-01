# Plots posteriors as volcano plots
library(ggplot2)
library(scales)
library(data.table)

setwd("~/tmpZika3")
burnin <- 10000
paramCols <- c("probMicro","baselineProb","epiStart","r0")

files <- list.files(pattern="\\.csv$")

allR0 <- NULL
allEpi <- NULL
allProb <- NULL
allBaseline <- NULL

for(file in files){
    print(file)
    tmp <- fread(file,data.table=FALSE)
    tmp <- as.data.frame(tmp[burnin:nrow(tmp),paramCols])

    tmp <- tmp[sample(seq(1,nrow(tmp),by=1),1000,replace=FALSE),]
    name <- strsplit(file,"_")
    name <- name[[1]][1]
    print(name)
    
    tmpR0 <- data.frame("value"=tmp$r0,"dat"=name)
    tmpEpi <- data.frame("value"=tmp$epiStart,"dat"=name)
    tmpProb <- data.frame("value"=tmp$probMicro,"dat"=name)
    tmpBase <- data.frame("value"=tmp$baselineProb,"dat"=name)

    allR0 <- rbind(allR0, tmpR0)
    allEpi <- rbind(allEpiStart, tmpEpi)
    allProb <- rbind(allProb, tmpProb)
    allBaseline <- rbind(allBaseline, tmpBase)
}

R0plot <- ggplot(allR0,aes(x=value))+
  stat_density(aes(ymax=..density..,ymin=-..density..,fill=dat,color=dat),trim=TRUE,geom="ribbon",position="identity")+
 # scale_x_continuous(limits=c(0,5))+
  facet_grid(~dat)+
  coord_flip()+ 
  ylab("Density") +
  ggtitle("Bite Rate, b")+
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
      legend.position="none"
  )
