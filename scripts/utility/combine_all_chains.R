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
library(zikaInfer)
library(nleqslv)

saveFile <- "~/Documents/Zika/combinedChains_21June2018.csv" ## Where to save the results
topDir <- "~/net/home/zika/outputs/" ## Where are the MCMC chains saved?
thin <- 100
## Thin the read in MCMC chains - this speeds up the code massively and prevents
## the creation of a massive final csv file
burnin <- 750000 ## Iterations to discard

runNames <- c("bahia/bahia","colombia_inc/colombia","northeast/northeast",
              "colombia_inc_month/colombia","colombia_peak/colombia","multi_3","multi_3_varied",
              "reports_2_inc","reports_2_weeks","reports_3","reports_3_inc",
              "reports_3_varied_inc",
              "riograndedonorte/riograndedonorte",
              "colombia_inc_suspected/colombia",
              "colombia_inc_confirmed/colombia",
              "northeast_deOliveira2017/northeast",
              "pernambuco_peak_confirmed/pernambuco",
              "riograndedonorte_confirmed/riograndedonorte")

correctedRunNames <- c("bahia/bahia"="Bahia, Brazil (reports)",
                       "colombia_inc/colombia"="Colombia",
                       "northeast/northeast"="Northeast Brazil",
                       "colombia_inc_month/colombia"="Colombia (monthly)",
                       "colombia_peak/colombia"="Colombia (peak only)",
                       "multi_3"="Brazil Model 3 states",
                       "reports_3" = "Reports 3 states, no incidence",
                       "reports_3_inc"="Reports 3 states",
                       "riograndedonorte/riograndedonorte"="Rio Grande do Norte, Brazil (reports)",
                       "colombia_inc_suspected/colombia"="Colombia (suspected)",
                       "colombia_inc_confirmed/colombia"="Colombia (confirmed)",
                       "northeast_deOliveira2017/northeast"="Northeast Brazil, Lancet",
                       "pernambuco_peak_confirmed/pernambuco"="Pernambuco (confirmed, reports)",
                       "riograndedonorte_confirmed/riograndedonorte"="Rio Grande do Norte, Brazil (confirmed)")

models <- c("model_1","model_2")
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
        tmpChain <- zikaInfer::load_mcmc_chains(dir,FALSE,FALSE,FALSE,thin,burnin,TRUE)
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
                tmpChain2 <- tmpChain1[,colnames(tmpChain1) %in% parNames]
                tmpChain2$R0 <- r0.vector(tmpChain1)
                tmpChain2$AR <- calculate_AR(tmpChain1$R0)
                tmpChain2$state <- states[i+1]
                stateChain <- rbind(stateChain,tmpChain2)
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

## Format and save
allChain <- allChain[allChain$runName %in% names(correctedRunNames),]
allChain$runName <- as.character(allChain$runName)
allChain$runName <- correctedRunNames[allChain$runName]
allChain$chain <- as.factor(allChain$chain)
allChain$state <- as.factor(allChain$state)
allChain$runName <- as.factor(allChain$runName)
allChain$version <- as.factor(allChain$version)

write.table(allChain,saveFile,sep=",",row.names=FALSE)
