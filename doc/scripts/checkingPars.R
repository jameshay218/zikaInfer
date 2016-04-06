# Which place to check
i <- 4

allDat <- read.csv("~/Dropbox/Zika/Data/allDat.csv")
print(places[i])
tmpDat <- allDat[allDat$local == places[i],]
filename <- as.character(places[i])
chains <- NULL
for(j in 1:3){
  tmp <- read.csv(paste(filename,j,"_chain.csv",sep=""))
  chains <- rbind(tmp,chains)
}


correctOrder <- sort(by(allDat[,"microCeph"],allDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
correctOrder <- as.factor(names(correctOrder))
places <- unique(allDat$local)
places <- correctOrder
adaptive <- 30000


xlabels <- NULL
xlab <- c("01","02","03","04","05","06","07","08","09","10","11","12")
for(i in xlab){
  xlabels <- c(xlabels,paste(i,"/2014",sep=""))
}
for(i in xlab){
  xlabels <- c(xlabels,paste(i,"/2015",sep=""))
}
xlabels <- c(xlabels,"01/2016","02/2016")


start <- seq(tmpDat[1,"start"],(nrow(tmpDat)/12) - 1/12,by=1/12)
end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
buckets <- cbind(start,end)
tmpDat <- tmpDat[,c("microCeph","births")]
tmpDat[,2] <- tmpDat[,2] - tmpDat[,1]

bestPars <- chains[which.max(chains$lnlike),2:(ncol(chains)-2)]
tmpName <- names(bestPars)
bestPars <- as.numeric(bestPars)
names(bestPars) <- names(pars[[3]])

tmpLn <- NULL
for(i in 1:200){
  tmpPar <- bestPars
  tmpPar["b"] <- i
  tmpLn[i] <- posterior(pars[[1]],pars[[2]],tmpPar,tmpDat,threshold=NULL,times=buckets)
}
plot(tmpLn)
print(which.max(tmpLn))
