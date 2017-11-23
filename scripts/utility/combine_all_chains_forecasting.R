#######################################
## JAMES HAY 20.11.2017 - jameshay218@gmail.com
## Reads all MCMC chains in the vector of directories. You should specify the "topDir", which contains all of the subfolders.
## "runNames" is then the full folder path (within this topDir) that contains each of the MCMC chains. Note that
## the script expects folders within these indicating which version of the model was run ie. model_1, model_2
## "correctRunNames" is simply used to reformat the run names into a neater form
## Adds states, run name, model version and chain number.
#############################################
library(data.table)
library(coda)
library(zikaProj)
library(nleqslv)

saveFile <- "~/Documents/Zika/combinedChains_forecasts.csv" ## Where to save the results
#topDir <- "~/net/home/zika/outputs/" ## Where are the MCMC chains saved?
topDir <- "/media/james/JH USB/Forecast_runs/"

microDat <- read.csv("~/Documents/Zika/Data/brazil/microceph_reports_2016.csv",stringsAsFactors=FALSE)
incDat <- read.csv("~/Documents/Zika/Data/brazil//zika_inc_reports.csv",stringsAsFactors=FALSE)
microDat <- microDat[microDat$local == "bahia",]
incDat <- incDat[incDat$local == "bahia",]

thin <- 10
## Thin the read in MCMC chains - this speeds up the code massively and prevents
## the creation of a massive final csv file
burnin <- 500000 ## Iterations to discard

runNames <- c("bahia_forecast_noreportingchange","bahia_abortion","bahia_forecast",
              "bahia_incPropn","bahia_propn","bahia_reduction")

correctedRunNames <- c("bahia_forecast_noreportingchange"="No ZIKV reporting change",
                       "bahia_abortion"="Abortions only",
                       "bahia_forecast"="Forecast analysis",
                       "bahia_incPropn"="ZIKV reporting change only",
                       "bahia_propn"="Microcephaly reporting change only",
                       "bahia_reduction"="Avoided births only")

extraMicroParNames <-c("chain","mode","maxWeek","lower","upper","range","max","tr1","tr2","tr3")
extraStateParNames <- c("R0","AR","chain","state", "zikv_report_rate_change")
allChain <- NULL

for(run in runNames){
    print(run)
    dir <- paste(topDir,run,sep="")
    
    if(!file.exists(dir)){
        print("Directory not found")
        next
    }
    ## Load MCMC chain from this run
    tmpChain <- zikaProj::load_mcmc_chains(dir,FALSE,FALSE,FALSE,thin,burnin,TRUE)
    parTab <- read.csv("~/net/home/zika/inputs/parTab_forecast.csv",stringsAsFactors=FALSE)

    parTab[parTab$names %in% c("L_H","N_H"),"values"] <- as.numeric(incDat[1,c("L_H","N_H")])
    parTab[parTab$names == "L_H","values"] <- parTab[parTab$names == "L_H","values"]*365
    f <- create_forecast_function(parTab, microDat, incDat=incDat, ts=seq(0,3003,by=1), FALSE)
    aborted <- NULL

    print("Calculating number of abortions")
    ## Calculating abortions
    for(i in 1:nrow(tmpChain)){
        pars <- get_index_pars(tmpChain,i)
        pars["baselineProb"] <- exp(pars["baselineProb"])
        aborted[i] <- sum(f(pars,TRUE)$aborted$aborted)
    }
    tmpChain$abortions <- aborted
    
    scale <- parTab[which(parTab$names=="propn" & parTab$fixed == 1),"values"]
  
    tmpChain$state <- "bahia"
    microParChain <- tmpChain[,colnames(tmpChain) %in% c(parTab[parTab$local == "all","names"],
                                                         extraMicroParNames)]
    microParChain$state <- "bahia"

    r0 <- r0.vector(tmpChain)
    attackRate <- calculate_AR(r0)
    tmpChain$R0 <- r0
    tmpChain$AR <- attackRate
    tmpChain$zikv_report_rate_change <- tmpChain$incPropn2/tmpChain$incPropn
    stateChain <- tmpChain
    
    meltMicroChain <- reshape2::melt(microParChain,id.vars=c("state","chain"))
    meltChain <- reshape2::melt(stateChain,id.vars=c("state","chain"))
    meltChain <- rbind(meltMicroChain,meltChain)
    meltChain <- cbind(runName=run,version="forecast",meltChain)
    allChain <- rbind(allChain, meltChain)
}

allChain <- allChain[allChain$runName %in% names(correctedRunNames),]
allChain$runName <- as.character(allChain$runName)
allChain$runName <- correctedRunNames[allChain$runName]
allChain$chain <- as.factor(allChain$chain)
allChain$state <- as.factor(allChain$state)
allChain$runName <- as.factor(allChain$runName)

write.table(allChain,saveFile,sep=",",row.names=FALSE)
