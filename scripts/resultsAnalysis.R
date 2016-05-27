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
  
  actualDat <- tmp[[2]]
  predictedInc <- tmp[[3]]
  probs <- tmp[[4]]
  predictedMicroceph <- tmp[[1]]
  r0s <- tmp[[5]]
  
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
 
  
   myPlot <- ggplot() + geom_point(dat=actualDat,aes(y=microCeph,x=meanDay),col="black") +
      geom_line(dat=microBounds,aes(y=value,x=Var2,group=Var1),linetype=2,col="red") +
      geom_line(dat=microMean,aes(y=value,x=Var1),col="red",lwd=1) +
    geom_line(dat=incBounds,aes(y=value,x=Var2,group=Var1),linetype=2,col="green") +
    geom_line(dat=incMean,aes(y=value,x=Var1),col="green",lwd=1)+
     theme_bw() +
     ylab("Per capita incidence (green)/microcephaly incidence (red)") +
     xlab("Day") +
     ggtitle(state)
   return(list(myPlot,density(r0s)))
}
   
  

