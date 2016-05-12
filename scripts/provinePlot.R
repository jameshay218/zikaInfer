library(rworldmap)
library(raster)
library(ggplot2)
library(dplyr)
library(rgeos)
source("~/Documents/zikaProj/doc/scripts/declareVariables.R")

# Get Brazilian state data to plot
brazil1 <- raster::getData("GADM",country="Brazil",level=1)
map <- fortify(brazil1)

map$id <- as.integer(map$id)
dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
map_df <- inner_join(map,dat,by="id")

levels(map_df$state) = states

# Read in results from model fitting ie. parameter estimates
greb <- read.csv("~/Dropbox/Zika/results5/finalstats2.csv")
greb[grep("r0",greb$X),"X"] <- "r0"
greb[grep("probMicro",greb$X),"X"] <- "probMicro"
greb[grep("baselineProb",greb$X),"X"] <- "baselineProb"
greb[grep("epiStart",greb$X),"X"] <- "epiStart"
colnames(greb)[ncol(greb)] <- "state"

# Reformat these estimates to be joined to the map plot data frame. This is ugly, but works...
plotLabels <- NULL
for(i in 1:length(unique(map_df$state))){
  place <- as.character(unique(map_df$state)[i])
  # If this state doesn't have good parameter estimates (ie. no epidemic), then we'll fill it in with NAs
  if(!(place %in% states_with_data)){
    omg <- matrix(nrow=1,ncol=nrow(tmp)+1)
    omg[1,] <-rep(NA,7)
  } else {
    tmp <- greb[as.character(greb$state)==place,c("X","Mean")]
    omg <- matrix(nrow=1,ncol=nrow(tmp)+1)
    R0 <- tmp[5,2]
    ## Also get attack rate estimate
    attackR <- nleqslv(0.8,zikaProj::simeq,R0=R0)$x
    omg[1,] <- c(tmp[,2],attackR)
  }
  colnames(omg) <- c("probMicro","baselineProb","epiStart","b","r0","lnlike","attackRate")
  omg <- as.data.frame(omg)
  omg$state <- as.character(unique(map_df$state)[i])
  plotLabels <- rbind(plotLabels,omg)
}

# Add parameter estimates to map plot data
map_df <- join(map_df,plotLabels,by="state")

# Finding centroids of each state as location for labels
centers <- data.frame(gCentroid(brazil1, byid = TRUE))
centers$state <- dat$state
levels(centers$state) = states
centers <- join(centers,plotLabels,by="state")

######
## Plots at state level for parameter estimates
######
r0Plot <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=r0),col="black")+
  scale_fill_gradient2(limits=c(1,8),mid="#132B43",low="#56B1F7",high="red",midpoint=5,breaks=seq(1,8,by=1)) +
  ggtitle("Basic Reproductive Number, R0")+
  #geom_text(data = centers, aes(label = r0, x = x, y = y), size = 3)+
  mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))

attackPlot <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=attackRate),col="black")+
  scale_fill_gradient(high="darkgreen",low="#e5f5e0") +
  ggtitle("Attack Rate")+
 # geom_text(data = centers, aes(label = attackRate, x = x, y = y), size = 3)+
  mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))

microPlot <-  ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=probMicro),col="black")+
  scale_fill_gradient(limits=c(0,0.5),low="white",high="#de2d26") +
 #geom_text(data = centers, aes(label = probMicro, x = x, y = y), size = 3)+
  ggtitle("Probability of Zika Associated Microcephaly")+
  mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))

basePlot <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=baselineProb),col="black")+
  scale_fill_gradient(low="white",high="#756bb1") +
 # geom_text(data = centers, aes(label = baselineProb, x = x, y = y), size = 3)+
  ggtitle("Baseline Probability of Microcephaly")+
  mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))

allRegions <- grid.arrange(r0Plot,attackPlot,microPlot,basePlot,ncol=2)


transposedCorrs <- NULL
colnames(stateCor)[4] <- "state"
for(i in 1:length(unique(stateCor$state))){
  place <- as.character(unique(stateCor$state)[i])
  tmp <- stateCor[as.character(stateCor$state)==place,c("value")]
  omg <- matrix(nrow=1,ncol=length(tmp))
  omg[1,] <- tmp
  colnames(omg) <- pairNames
  omg <- as.data.frame(omg)
  omg$state <- place
  transposedCorrs <- rbind(transposedCorrs,omg)
}
ohwow <- NULL
for(i in 1:length(transposedCorrs$state)){
  ohwow[i] <- names(which(transposedCorrs$state[i]==country_names))
}
transposedCorrs$state <- ohwow

map_df <- join(map_df,transposedCorrs,by="state")


####
## Make all correlation plots
pair1P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=R0.bp),col="black")+
  scale_fill_gradient(low="white",high=c25[1],limits=c(0,1)) +
  ggtitle("R0/bp Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme
pair2P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=R0.pm),col="black")+
  scale_fill_gradient(high="white",low=c25[2],limits=c(-1,0)) +  ggtitle("R0/pm Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme
pair3P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=R0.start),col="black")+
  scale_fill_gradient(low="white",high=c25[3],limits=c(0,1)) +  ggtitle("R0/epiStart Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme
pair4P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=pm.bp),col="black")+
  scale_fill_gradient(high="white",low=c25[4],limits=c(-1,0)) +  ggtitle("pm/bp Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme
pair5P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=pm.start),col="black")+
  scale_fill_gradient(high="white",low=c25[5],limits=c(-1,0)) +  ggtitle("pm/epiStart Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme
pair6P <- ggplot() + geom_map(data=map_df,map=map,aes(x=long,y=lat,map_id=id,group=group,fill=bp.start),col="black")+
  scale_fill_gradient(low="white",high=c25[6],limits=c(0,1)) +  ggtitle("bp/epiStart Correlation")+
  theme(plot.title=element_text(size=14))+
  mytheme

allRegionsCor <- grid.arrange(pair2P,pair1P,pair3P,pair4P,pair5P,pair6P,ncol=3)

