
pdf("traces.pdf")

plot(as.mcmc(chain[,c(indices+1,ncol(chain))]))
dev.off()

tmpChain <- chain[,c("mean","var","scale")]
samples <- sample(nrow(tmpChain),20)
index <- 1
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
  if(index ==1) plot(probs,type='l',ylim=c(0,0.1))
  else lines(probs)
  index <- index+1
}


tmp <- plot_best_pars(chain,realDat,parTab,pars[[1]],"bahia",1,1000)
tmp1 <- tmp[[2]]
tmp2 <- tmp[[3]]
probs <- tmp[[4]]
tmp <- tmp[[1]]

probBounds <- sapply(unique(probs$weeks), function(x) quantile(probs[probs$weeks==x,"prob"], c(0.025,0.5,0.975)))
plot(main="Microcephaly Risk Curve",probBounds[3,],type='l',lty=2,col="blue",ylim=c(0,0.1),ylab="Probability of microcephaly given infection",xlab="Week of Gestation")

lines(probBounds[1,],lty=2,col="blue")
abline(v=16,lty=2)
lines(probBounds[2,],lwd=2,col="blue")

ohdear <- sapply(unique(tmp$day), function(x) quantile(tmp[tmp$day==x,"number"], c(0.025,0.5,0.975)))

tmp2$times <- as.integer(tmp2$times)
ohdear1 <- sapply(unique(tmp2$times), function(x) quantile(tmp2[tmp2$times==x,"I_H"], c(0.025,0.5,0.975)))
plot(main="Bahia", ylab="Per capita incidence (green)/ microceph incidence (red)", 
     xlab="Day",tmp1[,"microCeph"]/tmp1[,"births"] ~ unique(tmp$day),
     ylim=c(0,0.08), pch=21,col="black",bg="black")

lines(ohdear[3,] ~ unique(tmp$day), col="red", lty=2)
lines(ohdear[1,] ~ unique(tmp$day),col="red",lty=2)
lines(ohdear[2,] ~ unique(tmp$day),col="red",lwd=2)
lines(ohdear1[1,],col="green",lty=2)
lines(ohdear1[3,],col="green",lty=2)
lines(ohdear1[2,],col="green",lwd=2)
