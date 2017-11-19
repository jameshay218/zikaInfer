## Reads all MCMC chains in the vector of directories. 
## Adds states, run name, model version and chain number.
#############################################
library(data.table)
library(coda)
library(zikaProj)
library(nleqslv)

saveFile <- "~/Documents/Zika/combinedChains.csv"
topDir <- "/media/james/JH USB/outputs/"
thin <- 100
burnin <- 750000

runNames <- c("bahia/bahia","colombia_inc/colombia","northeast/northeast")

correctedRunNames <- c("bahia/bahia"="Brazil, Bahia (reports)",
                       "colombia_inc/colombia"="Colombia",
                       "northeast/northeast"="Northeast Brazil")


models <- c("model_1")
extraMicroParNames <-c("chain","mode","maxWeek","lower","upper","range","max","tr1","tr2","tr3")
extraStateParNames <- c("R0","AR","chain","state")

allChain <- NULL
for(run in runNames){
    print(run)
    for(model in models){
        print(model)
        dir <- paste(topDir,run,"/",model,sep="")
        
        if(!file.exists(dir)){
            print("Directory not found")
            next
        }

        ## Load MCMC chain from this run
        tmpChain <- zikaProj::load_mcmc_chains(dir,FALSE,FALSE,FALSE,thin,burnin,TRUE)
        parTab <- read_inipars(dir)
        scale <- parTab[which(parTab$names=="propn" & parTab$fixed == 1),"values"]
        print(paste0("Fixed reporting proportion used: ", scale))
        tmpChain <- cbind(tmpChain, extra_microceph_calculations(tmpChain,1,0.001,scale))

        ## Find unique states
        states <- unique(parTab[parTab$local != "all","local"])

        ## If this was a run with multiple states, we need to analyse parameters for each in isolation
        if(length(states) > 1){
            ## Isolate index for this state
            indices <- 1:(length(states)-1)
            parNames <- unique(parTab[parTab$local != "all","names"])
            stateParNames <- c(sapply(indices,function(x) paste(parNames,".",x,sep="")))
            parNames <- c(parNames, stateParNames,extraStateParNames)
            
            indices <- 0:(length(states)-1)
            stateChain <- NULL

            ## Get parameters related to all locations
            microParChain <- tmpChain[,colnames(tmpChain) %in% c(parTab[parTab$local == "all","names"],
                                                                 extraMicroParNames)]
            microParChain$state <- "all"

            ## For each state/location, isolate its parameters
            for(i in indices){
                print(states[i+1])
                tmpChain1 <- isolate_location(states[i+1],i,tmpChain,parTab,"chain")
                tmpChain1 <- tmpChain1[,colnames(tmpChain1) %in% parNames]
                tmpChain1$R0 <- r0.vector(tmpChain1)
                tmpChain1$AR <- calculate_AR(tmpChain1$R0)
                stateChain <- rbind(stateChain,tmpChain1)
            }
            ## If only one state/location, just deal with this on its own
        } else {
            tmpChain$state <- states
            microParChain <- tmpChain[,colnames(tmpChain) %in% c(parTab[parTab$local == "all","names"],
                                                                 extraMicroParNames)]
            microParChain$state <- states

            r0 <- r0.vector(tmpChain)
            attackRate <- calculate_AR(r0)
            tmpChain$R0 <- r0
            tmpChain$AR <- attackRate
            stateChain <- tmpChain
        }
            
        meltMicroChain <- reshape2::melt(microParChain,id.vars=c("state","chain"))
        meltChain <- reshape2::melt(stateChain,id.vars=c("state","chain"))
        meltChain <- rbind(meltMicroChain,meltChain)
        meltChain <- cbind(runName=run,version=model,meltChain)
        allChain <- rbind(allChain, meltChain)
    }
}
allChain$runName <- correctedRunNames[allChain$runName]
allChain$chain <- as.factor(allChain$chain)
allChain$state <- as.factor(allChain$state)
allChain$runName <- as.factor(allChain$runName)
allChain$version <- as.factor(allChain$version)

write.table(allChain,saveFile,sep=",",row.names=FALSE)
