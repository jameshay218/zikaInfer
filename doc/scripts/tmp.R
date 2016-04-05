
chains <- NULL
for(i in 1:3){
  filename <- paste("pernambuco",i,"_chain.csv",sep="")
  print(filename)
  tmp <- fread(filename,data.table=FALSE)
  tmp <- tmp[30000:nrow(tmp),c("probMicro","r0","baselineProb","epiStart")]
  print(nrow(tmp))
  tmp <- tmp[seq(1,nrow(tmp),by=10),]
  chains[[i]] <- as.mcmc(tmp)
}
plot(as.mcmc.list(chains))

nchains <- 2

setwd("~/tmpZika5")
allDat <- read.csv("~/Dropbox/Zika/Data/allDat.csv")
places <- unique(allDat$local)
allStats <- NULL
for(i in 1:length(places)){
  chain <- NULL
  for(j in 1:nchains){
    setwd("~/tmpZika5")
    y <- paste(places[i],j,"_chain.csv",sep="")
    filename <- places[i]
    tmpChain <- read.csv(y)
    tmpChain <- tmpChain[(burnin+adaptive):nrow(tmpChain),c("probMicro","baselineProb","epiStart","b","r0","lnlike")]
    chain[[j]] <- as.mcmc(tmpChain)
  }
  
  setwd("~/Dropbox/Zika/results4")
  pdf(paste(filename,".pdf",sep=""))
  plot(mcmc.list(chain))
  dev.off()
  
  stats <- summary(mcmc.list(chain))
  toSave <- cbind(stats$statistics, stats$quantiles)
  toSave <- as.data.frame(toSave)
  toSave <- cbind(toSave,"local"=rep(as.character(places[i]),nrow(toSave)))
  
  allStats[[i]] <- toSave
}
setwd("~/Dropbox/Zika/results4")
finalStats <- do.call("rbind",allStats)
write.csv(finalStats,"finalstats2.csv")
