# Which place to check
index <- 4

allDat <- read.csv("~/Dropbox/Zika/Data/allDat07/04.16.csv")
places <- unique(allDat$local)
correctOrder <- sort(by(actualDat[,"microCeph"],actualDat[,"local"],sum,simplify=TRUE),decreasing=TRUE)
correctOrder <- as.factor(names(correctOrder))
tmpDat <- allDat[allDat$local == correctOrder[index],]

N_H <- tmpDat[1,"NH"]
L_H <- tmpDat[1,"LH"]

pars <- setupListPars(duration=end[length(end)], N_H = N_H, N_M = N_H*densityM)  

pars <- setupListPars(duration=end[length(end)], N_H = N_H, N_M = N_H*densityM)  
pars[[1]][1] <- 4
pars[[3]]["sampFreq"] <- 30
pars[[3]]["probMicro"] <- 0.2
pars[[3]]["baselineProb"] <- 0.01
pars[[3]]["epiStart"] <- 1.2
pars[[3]]["L_H"] <- L_H
pars[[3]]["constSeed"] <- 0

filename <- as.character(correctOrder[index])
chains <- NULL
for(j in 1:3){
  tmp <- fread(paste(filename,"_",j,"_chain.csv",sep=""),data.table=FALSE)
  chains <- rbind(tmp,chains)
}


start <- seq(0,(nrow(tmpDat)/12) - 1/12,by=1/12)
end <- seq(1/12,nrow(tmpDat)/12,by=1/12)
buckets <- cbind(start,end)

tmpDat <- tmpDat[,c("microCeph","births")]

bestPars <- chains[which.max(chains$lnlike),2:(ncol(chains)-2)]
tmpName <- names(bestPars)
bestPars <- as.numeric(bestPars)
names(bestPars) <- names(pars[[3]])

tmpLn <- NULL
for(i in 1:200){
  tmpPar <- bestPars
  tmpPar["epiStart"] <- i/100
  tmpLn[i] <- posterior(pars[[1]],pars[[2]],tmpPar,tmpDat,threshold=NULL,times=buckets)
}
plot(tmpLn)
print(which.max(tmpLn))
tmpPar <- bestPars
tmpPar["epiStart"] <- which.max(tmpLn)/100
p <- plot_line_dat(list(pars[[1]],pars[[2]],bestPars),tmpDat)
