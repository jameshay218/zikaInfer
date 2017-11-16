## Loads all MCMC chains from a top level directory (go through each sub directory)
load_all_chains <- function(thin=1,burnin=100000,version=1,microceph_limit=0.001){
  directories <- list.dirs(full.names=FALSE,recursive=FALSE)
  topDir <- getwd()
  allChain <- NULL
  for(directory in directories[1:length(directories)]){
    print(directory)
    setwd(topDir)
    setwd(directory)
    tmpChains <- load_mcmc_chains(TRUE,FALSE,FALSE,thin,NULL,burnin)
    for(i in 1:length(tmpChains)){
      r0 <- r0.vector(tmpChains[[i]])
      attackRate <- calculate_AR(r0)
      if(version == 1) { 
        mode <- calculate_gamma_mode(tmpChains[[i]]$mean,tmpChains[[i]]$var)
        maxWeek <- get_max_micro_week(tmpChains[[i]])
        maxP <- get_max_micro(tmpChains[[i]],directory)
      } else{
        mode <- 0
        maxWeek <- 0
        maxP <- 0
      } 
      lower <- get_lower_micro_bound(tmpChains[[i]],directory,microceph_limit)
      upper <- get_upper_micro_bound(tmpChains[[i]],directory,microceph_limit)
      range <- upper - lower
      tmpChains[[i]] <- cbind(tmpChains[[i]],R0=r0,AR=attackRate,Mode=mode,maxWeek=maxWeek,lower=lower,upper=upper,range=range,max=maxP,chain=i,state=directory)
          }
    tmpChain <- do.call("rbind",tmpChains)
    allChain <- rbind(allChain,tmpChain)
  }
  setwd(topDir)
  allChain
}

density_plot <- function(meltedChain, variableName, xlimits=NULL, title="test",breaks=NULL,labels=NULL){
  tmpDat <- meltedChain[meltedChain$variable==variableName,]
  plot <- ggplot(tmpDat,aes(x=value))+
      stat_density(aes(ymax=..density..,ymin=-..density..,fill=state,color=state),trim=TRUE,geom="ribbon",position="identity")+
      facet_grid(state~.)+

      ylab("") +
      ggtitle(title)+
      xlab("Value")+
      theme_bw() +
      theme(
        plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),
        legend.position="none",
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()
      )
    if(!is.null(xlimits)){
      if(!is.null(breaks) & !is.null(labels)) plot <- plot + scale_x_continuous(limits=xlimits,breaks=breaks,labels=labels)+
          theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),
                legend.position="none",
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                axis.text.x=element_text(angle=45,hjust=1))
      else plot <- plot +
          scale_x_continuous(limits=xlimits)
    }
       return(plot)
}

plot_best_trajectory_single1 <- function(local, chain=NULL, realDat=NULL, 
                                         parTab=NULL, t_pars=NULL, runs=100,incDat=NULL, 
                                         ylabel=TRUE,xlabel=TRUE,
                                         mcmcPars=c("burnin"=50000,"adaptive_period"=100000,"thin"=50),ylimM=NULL, ylimI=NULL, title=NULL){
  allDat <- plot_setup_data(chain,realDat, parTab, t_pars,local,runs)
  bestMicro <- allDat$bestMicro
  bestInc <- allDat$bestInc
  incBounds <- allDat$incBounds
  microBounds <- allDat$microBounds
  dat <- allDat$data
  
  xlim <- c(min(dat[,"startDay"]),max(dat[,"endDay"]))
  xlim <- c(0,1500)
  
  quantiles <- unique(microBounds[,"quantile"])
  botM <- microBounds[microBounds[,"quantile"]==quantiles[1],c("time","micro")]
  topM <- microBounds[microBounds[,"quantile"]==quantiles[2],c("time","micro")]
  
  botI <- incBounds[incBounds[,"quantile"]==quantiles[1],c("time","inc")]
  topI <-incBounds[incBounds[,"quantile"]==quantiles[2],c("time","inc")]
  polygonM <- create_polygons(botM, topM)
  polygonI <- create_polygons(botI, topI)
  
  xlabs <- generate_x_labels(xlim[1],xlim[2])
  myPlot <- microceph_plot(dat,microBounds,bestMicro,polygonM,local,xlim,ylimM,xlabs) + ggtitle(title) +
    theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.y=element_blank(), axis.title.x=element_blank())
  incPlot <- inc_plot(incBounds,bestInc,polygonI,ylimI,xlim,incDat) + ggtitle(title) +
    theme(axis.text.x=element_text(size=12),axis.text.y=element_text(size=12),axis.title.y=element_blank(), axis.title.x=element_blank())
  
  if(!ylabel){
    myPlot <- myPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,0,0),"cm"))
    incPlot <- incPlot + ylab("") + theme(plot.margin=unit(c(0,0.4,0,0),"cm"))
  }
  if(!xlabel) myPlot <- myPlot + xlab("")
  
  combined <- overlapPlots(myPlot,incPlot,FALSE)
  
  return(combined)
}

generate_country_boundaries <- function(country="Brazil"){
  library(rworldmap)
  library(raster)
  brazil1 <- raster::getData("GADM",country=country,level=1)
  map <- fortify(brazil1)
  
  map$id <- as.integer(map$id)
  dat <- data.frame(id=1:(length(brazil1@data$NAME_1)),state=brazil1@data$NAME_1)
  map_df <- inner_join(map,dat,by="id")
  return(map_df)
}



isolate_state <- function(local, i, chains, parTab,microceph_limit=0.001,version=1){
  pars <- as.numeric(chains[100,])
  names(pars) <- colnames(chains)
  
  state_pars <- parTab[parTab$local==local,"names"]
  if(i >= 1) state_pars <- paste(state_pars,".",i,sep="")
  state_pars <- pars[state_pars]
  stateNames <- c(names(state_pars),parTab[parTab$local=="all","names"])
  
  blankNames <- parTab[parTab$local==local,"names"]
  blankNames <- c(blankNames,parTab[parTab$local=="all","names"])
  
  tmpChain <- chains[,stateNames]
  colnames(tmpChain) <- blankNames
  
  r0 <- r0.vector(tmpChain)
  attackRate <- calculate_AR(r0)
  if(version == 1){
    mode <- calculate_gamma_mode(tmpChain$mean,tmpChain$var)
    lower <- get_lower_micro_bound(tmpChain,local,microceph_limit)
    tmpChain <- cbind(chains, mode=mode,lower=lower,R0=r0,AR=attackRate,state=local)
  } else {
   tmpChain <- cbind(chains[,colnames(chains) %in% parNames],R0=r0,AR=attackRate,state=local)
  }
  
  return(tmpChain)
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

generate_country_plot <- function(var,map_df,centers, title){
  brazil1 <- raster::getData("GADM",country="Brazil",level=1)
  map <- fortify(brazil1)
  
  map$id <- as.integer(map$id)
 # centers[,var] <- signif(centers[,var],3)
  plot <- ggplot() + 
    geom_map(data=map_df,map=map,aes_string(x="long",y="lat",map_id="id",group="group",fill=var),col="black")+
    ggtitle(title)+
    geom_text(data = centers, aes_string(label = var, x = "x", y = "y"), size = 3)+
    mytheme + theme(plot.title=element_text(size=14),legend.text=element_text(size=16))
}


generate_mcmc_summaries <- function(chain,states, version="single"){
  summaries <- NULL
  if(version == "single"){
    for(state in states){
      tmpChain <- chain[chain$state==state,]
      tmpChain[is.na(tmpChain)] <- 0
      print(state)
      tmpSum <- summary(as.mcmc(tmpChain[,1:(ncol(tmpChain)-1)]),na.rm=TRUE)
      tmpSum <- data.frame(names=rownames(tmpSum$statistics),tmpSum$statistics,tmpSum$quantiles,state=state,stringsAsFactors = FALSE)
      summaries <- rbind(summaries,tmpSum)
    }
  } else {
    indices <- seq(1,length(states),by=1)
    for(i in seq(0,length(indices)-1,by=1)){
      print(states[i+1])
      tmpChain <- isolate_state(states[i+1],i,chain,parTab)
      tmpChain[is.na(tmpChain)] <- 0
      tmpSum <- summary(as.mcmc(tmpChain[,1:(ncol(tmpChain)-1)]),na.rm=TRUE)
      tmpSum <- data.frame(names=rownames(tmpSum$statistics),tmpSum$statistics,tmpSum$quantiles,state=states[i+1],stringsAsFactors = FALSE)
      summaries <- rbind(summaries,tmpSum)
    }
  }
  return(summaries)
}


generate_country_plot_labels <- function(plot_pars,map_df,states, summaries){
  plotLabels <- NULL
  npar <- length(plot_pars)
  
  places <- as.character(unique(map_df$state))
  places <- places[!(places %in% states)]
  for(place in states){
   
    tmp1 <- summaries[summaries$state==place,c("names","Mean","X2.5.","X97.5.")]
    tmp <- matrix(nrow=1,ncol=npar)
    colNames <- tmp1[tmp1$names %in% plot_pars,"names"]
  
    means <- signif(tmp1[tmp1$names %in% plot_pars,"Mean"],3)
    lowerCI <- signif(tmp1[tmp1$names %in% plot_pars,"X2.5."],3)
    upperCI <- signif(tmp1[tmp1$names %in% plot_pars,"X97.5."],3)
    
    #labels <- NULL
    #for(i in 1:length(means)){
    #  labels[i] <- paste(means[i]," (",lowerCI[i],"-",upperCI[i],")",sep="")
    #}
    #tmp[1,] <- labels
    tmp[1,] <- means
    
    colnames(tmp) <- colNames
    tmp <- as.data.frame(tmp)
    tmp$state <- place
    plotLabels <- rbind(plotLabels,tmp)
  }
  for(place in places){
  
    tmp <- matrix(nrow=1,ncol=npar)
      tmp[1,] <- rep(NA,npar)
        colnames(tmp) <- colNames
    tmp <- as.data.frame(tmp)
    tmp$state <- place
    plotLabels <- rbind(plotLabels,tmp)
  }
  return(plotLabels)

}


remake_multi_plot <- function(filename, mcmcPars,stateNames,incDatFile=NULL){
  t_pars <- c("dur"=3000,"step"=1)
  if(!is.null(incDatFile)){
    incDat <- read.csv(incDatFile)
    incDat <- incDat[incDat$local %in% stateNames,]
  } else {
    incDat <- NULL
  }
  realDat <- read.csv("~/net/home/zika/data/microcephaly_data.csv")
  realDat <- realDat[realDat$local %in% stateNames,]
  parTab <- setupParTable(1,realDat,sharedProb=TRUE)
  result <- read.csv(filename)
  result <- result[result$sampno >= (mcmcPars["adaptive_period"]+mcmcPars["burnin"]),]
  pdf("remade_multi_plot.pdf",width=20,height=20)
  plot(plot_best_trajectory_multi(result,realDat,parTab,t_pars,10,incDat,mcmcPars))
  dev.off()
}


make_micro_dat <- function(chain, samp_no){
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
  means <- means/max(upper)
  lower <- lower/max(upper)
  upper <- upper/max(upper)
  dat <- data.frame(weeks=1:(40*7),means,upper,lower)
  return(dat)
}
