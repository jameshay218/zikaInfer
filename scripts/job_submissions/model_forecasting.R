model_forecasting <- function(runName, chainNo, microDat, incDat, parTab, mcmcPars1, mcmcPars2){
    ts <- seq(0,3003,by=1)
    ## Create directory for results if doesn't exist and go into this
    newModelDir <- paste(getwd(), "/outputs/",runName,sep="")
    if(!dir.exists(newModelDir)) dir.create(newModelDir,recursive=TRUE)
    
    ## Create filenames to save chains
    arbit_filename <- paste(sprintf("%s/outputs/%s/%s_%d", getwd(),runName,runName,chainNo))
    filename1 <- paste(sprintf("%s/outputs/%s/%s_%d", getwd(),runName,runName,chainNo),"univariate",sep="_")
    filename2 <- paste(sprintf("%s/outputs/%s/%s_%d", getwd(),runName, runName,chainNo),"slice",sep="_")
    
    create_f <- function(parTab,data,PRIOR_FUNC,...){
        microDat <- data
        incDat <- incDat

        startDays <- microDat$startDay
        endDays <- microDat$endDay
        buckets <- microDat$buckets
        births <- microDat$births
        microCeph <- microDat$microCeph
        
        zikv <- incDat$inc
        nh <- incDat$N_H
        inc_buckets <- incDat$buckets
        inc_start <- incDat$startDay
        inc_end <- incDat$endDay
        
        names_pars <- parTab$names
        f <- function(values){
            names(values) <- names_pars
            lik <- posterior_known_inc_seir(values, startDays, endDays,
                                            buckets, microCeph, births,
                                            zikv, nh, inc_buckets,
                                            inc_start, inc_end)
            if(!is.null(PRIOR_FUNC)) lik <- lik + PRIOR_FUNC(values)
            return(lik)
        }
    }
    wow <- lazymcmc::run_MCMC(parTab,data=microDat,mcmcPars1,filename=filename1,
                              CREATE_POSTERIOR_FUNC = create_f,mvrPars=NULL,PRIOR_FUNC=NULL,
                              OPT_TUNING=0.2,incDat=incDat,peakTimes=NULL,ts=ts)
    chain <- read.csv(wow$file)
    startTab <- parTab
    startTab$values <- get_best_pars(chain)
    chain <- chain[chain$sampno > mcmcPars1["adaptive_period"],2:(ncol(chain)-1)]

    covMat <- cov(chain)
    mvrPars <- list(covMat,2.38/sqrt(nrow(parTab[parTab$fixed==0,])),w=0.9)
    wow1 <- lazymcmc::run_MCMC(startTab,data=microDat,mcmcPars2,filename=filename2,
                               CREATE_POSTERIOR_FUNC = create_f,mvrPars=mvrPars,PRIOR_FUNC=NULL,
                               OPT_TUNING=0.2,incDat=incDat,peakTimes=NULL,ts=ts)
    chain1 <- read.csv(wow1$file)
    pdf(paste0(filename1,".pdf"))
    plot(coda::as.mcmc(chain1))
    dev.off()
    
    return(TRUE)
}
