model_fitting_script <- function(runName,
                                 chainNo,
                                 version=1,
                                 stateNames,
                                 microDat,
                                 incDat,
                                 preParTab=NULL,
                                 allowablePars=NULL,
                                 incStates=NULL,
                                 mcmcPars=NULL,
                                 mcmcPars2=NULL,
                                 allPriors=NULL,
                                 useInc = TRUE,
                                 usePeakTimes=FALSE,
                                 covMat=NULL,
                                 peakTimeRange=60,
                                 peakTime=930,
                                 sim=FALSE,
                                 predict=FALSE,
                                 microChain=NULL,
                                 prePeakTimes=NULL,
                                 normLik=NULL,
                                 stateWeights=NULL,
                                 extra_unfixed=NULL,
                                 sharedProb=TRUE,
                                 mosquito_lifespan=5
                                 ){
    ## Create directory for results if doesn't exist and go into this
    wdRunName <- runName
    if(length(stateNames)==1) wdRunName <- paste(runName,stateNames[1],sep="/")
    
    newModelDir <- paste(getwd(),"/outputs/",wdRunName,"/model_",version,sep="")
    if(!dir.exists(newModelDir)) dir.create(newModelDir,recursive=TRUE)
    
    ## Create filenames to save chains
    arbit_filename <- sprintf("%s/outputs/%s/model_%d/%s_%d", getwd(),wdRunName, version,runName,chainNo)
    filename1 <- sprintf("%s/outputs/%s/model_%d/%s_%d", getwd(),wdRunName, version,runName,chainNo)
    filename2 <- sprintf("%s/outputs/%s/model_%d/%s_%d", getwd(),wdRunName, version,runName, chainNo)
    
##########################################
    ## Data set up
##########################################
    ts <- seq(0,3003,by=1)
    microDat <- microDat[microDat$local %in% stateNames,] # Get subset of microcephaly data

    ## If trying to specify subset of states for incidence data, do this here to subset incidence data
    use_incStates <- stateNames
    if(!is.null(incDat) && !is.null(incStates)) use_incStates <- incStates
    if(!is.null(incDat)) incDat <- incDat[incDat$local %in% use_incStates,]
    
    ## If no parameter table provided, set one up. Otherwise, skip this and use the precomputed parameter
    ## table as the starting point
    parTab <- preParTab
    if(is.null(preParTab)){
        ## If state weightings have been specified, need to flag this for parameter table setup
        useStateWeights <- FALSE
        if(!is.null(stateWeights)) useStateWeights <- TRUE
        data(locationInfo)
        ## Set up parameter table
        ##parTab <- partab_setup(stateNames=stateNames,version=version,
        ##                       realDat=microDat,useInc=useInc,
        ##                       allowablePars=NULL,sharedProb=sharedProb,
        ##                       normLik=normLik,stateWeights=useStateWeights)
        parTab <- parTabSetup(locationNames=stateNames, version=version, locationInfo=locationInfo,
                              useInc=useInc,allowablePars=NULL, sharedProb=sharedProb,
                              normLik=normLik, locationWeights=useStateWeights, forecast=FALSE)

        ## If a vector of state weights have been provided, need to set these explicitly in the parameter table
        if(!is.logical(stateWeights)){
            for(place in stateNames){
                parTab[parTab$local == place & parTab$names == "location_weight","values"] <- stateWeights[place]
            }
        }

        ## If we're excluding some states from incidence data, fix their incidence specific parameters
        if(!is.null(incStates)){
            parTab[!(parTab$local %in% use_incStates) & parTab$names == "incPropn","values"] <- 1
            parTab[!(parTab$local %in% use_incStates) & parTab$names == "incPropn","fixed"] <- 1
            parTab[!(parTab$local %in% use_incStates) & parTab$names == "baselineInc","fixed"] <- 1
            parTab[!(parTab$local %in% use_incStates) & parTab$names == "baselineInc","values"] <- 0
        }
    }
    ## Get parameter table subset of used states
    parTab <- parTab[parTab$local %in% c("all",stateNames),]
    parTab[parTab$names == "propn","upper_bounds"] <- 2

    ## If not provided, generate a data frame of allowable peak times.
    peakTimes <- NULL
    if(usePeakTimes){
        if(is.null(prePeakTimes)){
            peakTimes <- data.frame("start"=peakTime-peakTimeRange/2,"end"=peakTime+peakTimeRange/2,"local"=as.character(stateNames),stringsAsFactors=FALSE)
        } else {
            peakTimes <- prePeakTimes
        }
    }

        
##################
    ## STARTING PARAMETERS
##################
    ## Generate random starting points
    ## If some extra unfixed parameters have been specified, unfix these in the parameter table
    if(!is.null(extra_unfixed)){
        parTab[parTab$names %in% extra_unfixed,"fixed"] <- 0
    }

    ## Mosquito lifespan
    parTab[parTab$names == "L_M","values"] <- mosquito_lifespan
    
    ## Starting parameters from parTab
    startTab <- parTab
    useIncDat <- NULL
    if(useInc) useIncDat <- incDat

    ## If not pre-provided, generate starting points
    if(is.null(preParTab)){
        ## If not allowable starting parameters provided, generate a table of these
        ## This is important if parameters contributing to R0 have been varied
        if(is.null(allowablePars)){
            message("Generating allowable starting parameters...")
            ## Max R0 is 7
            allowablePars <- generateAllowableParams(peakTime=peakTime,peakTimeRange=peakTimeRange,stateNames=stateNames,startTab,"allowablePars.csv", R0max=7, peakTimings=prePeakTimes)
        }
        write.table(allowablePars, paste0(arbit_filename,"_allowablePars.csv"),sep=",")
        
        startTab <- generateStartingParTab(parTab, peakTimes, TRUE, allowableParsFile=allowableParsFile, allowablePars)
        ## If using incidence data, starting point for incidence proportion can be based on the data itself
        if(useInc) {
            useIncDat <- incDat
            for(local in unique(incDat$local)){
                startTab[startTab$local == local & startTab$names=="incPropn","values"] <- max(incDat[incDat$local == local,"inc"])/incDat[incDat$local==local,"N_H"][1]
                startTab[startTab$local == local & startTab$names=="baselineInc","values"] <- 0.0001
            }
        } 
    }
    
    iniParFile <- paste(arbit_filename, "_inipars.csv",sep="")
    write.table(startTab, iniParFile, row.names=FALSE,sep=",")
    truePars <- NULL
    if(sim){
        ## Write starting data and parameters used to generate data
        write.table(microDat,paste(arbit_filename,"_microDat.csv",sep=""),sep=",",row.names=FALSE)
        write.table(incDat,paste(arbit_filename,"_incDat.csv",sep=""),sep=",",row.names=FALSE)
        write.table(parTab,paste(arbit_filename,"_simpars.csv",sep=""),sep=",",row.names=FALSE)
        truePars <- parTab$values
        names(truePars) <- parTab$names
    }

###############################################################
    ## UNIVARIATE SAMPLER - CONFIGURATION CHAIN
###############################################################
    realDat <- microDat
    if(predict) realDat <- NULL
    if(is.null(covMat)){
        message("Starting MCMC chain")
        print(mcmcPars)
        result <- lazymcmc::run_MCMC(parTab=startTab,data=realDat, mcmcPars=mcmcPars, filename=filename1, 
                                      CREATE_POSTERIOR_FUNC=create_posterior, mvrPars=NULL, PRIOR_FUNC=NULL,
                                      OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=peakTimes)
        message("Univariate chain done!")
        ## Use results from first chain as starting conditions and widths
        chain <- read.csv(result$file)
        chain <- chain[chain$sampno >= mcmcPars["adaptive_period"],]
        covMat <- cov(chain[,2:(ncol(chain)-1)])
        startTab$values <- get_best_pars(chain)
        
        
        if(is.null(covMat)) message("Covmat broken...")
    }


###############################################################
    ## MULTIVARIATE SAMPLER - FINAL RESULTS
###############################################################
    print(filename2)
    mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.8)
    result <- lazymcmc::run_MCMC(parTab=startTab,data=realDat, mcmcPars=mcmcPars2, filename=filename2, 
                                 CREATE_POSTERIOR_FUNC=create_posterior, mvrPars=mvrPars, PRIOR_FUNC=NULL,
                                 OPT_TUNING=0.2, ts=ts,incDat=incDat,peakTimes=peakTimes)
    result <- read.csv(result$file)
    result <- result[result$sampno >= mcmcPars2["adaptive_period"],]

###############################################################
    ## PLOT RESULTS
###############################################################
    ## Plot MCMC chains
    plot_MCMC(list(result), parTab,paste(arbit_filename,"_mcmc.pdf",sep=""))
    
     ## Plot estimated microcepahly curves
    to.png(plot(plot_random_microceph_curves(result, 50)),paste(arbit_filename,"_microceph.png",sep=""))
    
     ## Plot trajectories
    weeks <- round(max(ts)/7)
    ## Plot model fit/trajectory
    #to.png(plot_best_trajectory_multi(result,microDat,parTab,ts,200,incDat,mcmcPars,startDay=0,months=NULL,weeks=weeks),paste(filename2,"_plot.png",sep=""))

    if(predict){
        microceph_births <- forecast_microceph(result, microChain,parTab, ts,200,"2013-01-01")
        write.table(microceph_births, paste(arbit_filename,"_predict.csv",sep=""),row.names=FALSE)
    }
    
    return(TRUE)
}


