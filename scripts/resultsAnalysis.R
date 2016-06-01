setupStuff <- function(){
  runName <- "real_3_all"
  allDatFile <- "allDat20.04.16.csv"
  realDat <- read.csv(allDatFile)
  parTab <- setupParTable(version=3,realDat=realDat)
  pars <- setupListPars(duration = 1200,N_H = 9000000,N_M=27000000,version=3)
  
  filename2 <- paste(sprintf("%s_%d", runName, chain),"B",sep="_")
  opt_freq <- 5000
  states <- c("all",
              "pernambuco",
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
  parTab <- parTab[parTab$local %in% states,]
  parTab[parTab$names=="propn","fixed"] <- 0
  peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
  testDat <- generate_multiple_data(pars[[1]],parTab,NULL)
  
  plotStates <- unique(parTab$local)[2:length(unique(parTab$local))]
  ps <- NULL
  for(i in 1:length(plotStates)){
    ps[[i]] <- plot_best_trajectories(chain, realDat,parTab,pars,plotStates[i],i-1,100)[[1]]
  }
  allPlots <- do.call("grid.arrange",c(ps,ncol=4))
}

plot_random_microceph_curves <- function(chain, runs){
  tmpChain <- chain[,c("mean","var","scale")]
  samples <- sample(nrow(tmpChain),runs)
  index <- 1
  myPlot <- ggplot()
  for(i in samples){
    tmpPars <- as.numeric(tmpChain[i,])
    gammaMean <- tmpPars[1]
    gammaVar <- tmpPars[2]
    
    rate <- gammaMean/gammaVar
    shape <- gammaMean*rate
    #probs <- dgamma(0:39, pars["shape"],pars["rate"])*pars["scale"]
    
    probs <- dgamma(0:39,shape,rate)*tmpPars[3]
    
    #probs <- c(rep(pars["shape"],13),rep(pars["rate"],13),rep(pars["scale"],14))
    probs[probs > 1] <- 1
    probs <- data.frame(probs=probs,week=seq(0,39,by=1))
    myPlot <- myPlot + geom_line(dat=probs, aes(x=week,y=probs))
    index <- index+1
  }
  return(myPlot)
}

plot_best_trajectories <- function(chain, realDat, parTab, pars, state, number, runs){
  tmp <- plot_best_pars(chain,realDat,parTab,pars[[1]],state,number,runs)
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
  actualDat <- tmp[[2]]
  predictedInc <- tmp[[3]]
  probs <- tmp[[4]]
  predictedMicroceph <- tmp[[1]]
  r0s <- tmp[[5]]
  
  xlabels <- c("01/2014","04/2014","07/2014","10/2014","01/2015","04/2015","07/2015","10/2015","01/2016")
 
   predictedInc$times <- as.integer(predictedInc$times)
  
  probBounds <- sapply(unique(probs$weeks), function(x) quantile(probs[probs$weeks==x,"prob"], c(0.025,0.5,0.975)))
  microBounds <- sapply(unique(predictedMicroceph$day), function(x) quantile(predictedMicroceph[predictedMicroceph$day==x,"number"], c(0.025,0.5,0.975)))
  microMean <- data.frame(value=microBounds[2,],Var1=actualDat[,"meanDay"])
  microBounds <- melt(microBounds[c(1,3),])
  
  incBounds <- sapply(unique(predictedInc$times), function(x) quantile(predictedInc[predictedInc$times==x,"I_H"], c(0.025,0.5,0.975)))
  incMean <- data.frame(value=incBounds[2,],Var1=seq(1,ncol(incBounds),by=1))
  incBounds <- melt(incBounds)
  
  actualDat[,"microCeph"] <- actualDat[,"microCeph"]/actualDat[,"births"]
  
  microBounds[,"Var2"] <- actualDat[,"meanDay"][microBounds[,"Var2"]]
  
  xlabBreaks <- actualDat[,"meanDay"][seq(1,length(actualDat[,"meanDay"]),by=3)]

   myPlot <- ggplot() + geom_point(dat=actualDat,aes(y=microCeph,x=meanDay),col="black",size=2) +
      geom_line(dat=microBounds,aes(y=value,x=Var2,group=Var1),lwd=0.5,linetype=2,col="red") +
      geom_line(dat=microMean,aes(y=value,x=Var1),col="red",lwd=0.5) +
        scale_y_continuous(limits=c(0,0.05))+
       scale_x_continuous(labels=xlabels,breaks=xlabBreaks)+
       theme_bw() +
       ylab("Per capita zika (green)/microcephaly (red)") +
       xlab("Date") +
       ggtitle(country_names[state]) + 
       theme(axis.text.x=element_text(hjust=1,angle=45,size=6),
             panel.grid.minor=element_blank(),
             plot.title=element_text(size=10),
             axis.text.x=element_text(size=6),
             axis.text.y=element_text(size=6),
             axis.title.x=element_text(size=8),
             axis.title.y=element_text(size=8))
   return(myPlot)
   
   incPlot <- ggplot() +
   geom_line(dat=incBounds,aes(y=value,x=Var2,group=Var1),linetype=2,col="green",size=0.5) +
     geom_line(dat=incMean,aes(y=value,x=Var1),col="green",lwd=0.5)+
     xlab("") +
     ylab("") +
     theme(
       element_rect(fill=NA),
     panel.grid=element_blank(),
     axis.text.y=element_text(colour="gray20",size=14),
     axis.line.y=element_line(colour="gray20"),
     axis.title.y=element_text(colour="gray20",size=16,angle=-90),
     axis.line.y=element_blank(),
     axis.line.x=element_blank(),
     axis.ticks.y=element_blank(),
     axis.ticks.x=element_blank(),
     axis.text.x=element_blank(),
     axis.text.y=element_blank()
     )
   
   return(list(myPlot, incPlot,density(r0s)))
}
   
  

