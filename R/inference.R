#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size. Random walk may be on a linear or log scale
#' @param values a vector of the parameters to be explored
#' @param lower_bounds a vector of the low allowable bounds for the proposal
#' @param upper_bounds a vector of the upper allowable bounds for the proposal
#' @param steps a vector of step sizes for the proposal
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
#' @useDynLib zikaProj
proposalfunction <- function(values, lower_bounds, upper_bounds,steps, index){
    mn <- lower_bounds[index]
    mx <- upper_bounds[index]

    rtn <- values
    
    x <- toUnitScale(values[index],mn,mx)

    ## 5th index is step size
    stp <- steps[index]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

    ## Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

    ## Cyclical boundary conditions
    ##if (x < 0) x <- 1 + x	
    ##if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}

#' Multivariate proposal function
#'
#' Given the current parameters and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param values the vector of current parameter values
#' @param fixed set of flags corresponding to the parameter vector indicating which parameters are fixed
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
#' @useDynLib zikaProj
mvr_proposal <- function(values, fixed, covMat){
    proposed <- values
    proposed[fixed] <- proposed[fixed] + mvrnorm(n=1,mu=rep(0,length(proposed[fixed])),covMat[fixed,fixed])
    return(proposed)
}




#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{proposalfunction}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param data the data against which the likelihood is calculated
#' @param ts vector of times to solve the ODE model over
#' @param y0s starting conditions for the ODE model
#' @param param_table a table of parameter data used for information such as bounds and prior function pointers.
#' @param mcmcPars a named vector with parameters for the MCMC procedure. Iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param mvrPars a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param incDat optional data frame of incidence data if including incidence data in the likelihood function
#' @param peakTimes optional parameter - data frame of peak times for Zika incidence for each state
#' @param allPriors flag of whether or not to use priors
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) the last used covarianec matrix; 3) the last used scale size
#' @export
#' @seealso \code{\link{posterior}}, \code{\link{proposalfunction}}
#' @useDynLib zikaProj
run_metropolis_MCMC <- function(data,
                                t_pars,
                                param_table,
                                mcmcPars=c("iterations"=1000,"popt"=0.44,"opt_freq"=50,"thin"=1,"burnin"=100,"adaptive_period"=100,"save_block"=500),
                                filename,
                                mvrPars=NULL,
                                incDat=NULL,
                                peakTimes=NULL,
                                allPriors=NULL,
                                ...
                                ){
    iterations <- mcmcPars["iterations"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    burnin <- mcmcPars["burnin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    save_block <- mcmcPars["save_block"]

    TUNING_ERROR<- 0.1

    ## Extract parameter table into individual vectors for speed
    non_fixed_params <- which(param_table$fixed==0)
    par_names <- param_table$names
    current_params <- start_values <- param_table$values
    par_labels <- param_table$local
    lower_bounds <- param_table$lower_bounds
    upper_bounds <- param_table$upper_bounds
    steps <- param_table$steps
    fixed <- param_table$fixed
    names(current_params) <- names(start_values) <- par_names

    ## Extract data into vectors
    startDays <- data[,"startDay"]
    endDays <- data[,"endDay"]
    buckets <- data[,"buckets"]
    microCeph <- data[,"microCeph"]
    births <- data[,"births"]
    data_locals <- data[,"local"]
    
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- nrow(param_table)

    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",param_table$names,"lnlike")
    ## Arrays to store acceptance rates
    if(is.null(mvrPars)) tempaccepted <- tempiter <- integer(all_param_length)
    else {
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]]
        w <- mvrPars[[2]]
        scale <- mvrPars[[3]]
        scaledCovMat <- covMat*scale
    }
    reset <- integer(all_param_length)
    reset[] <- 0

    ## Create empty chain to store every iteration for the adaptive period
    empty_chain <- chain <- matrix(nrow=opt_freq,ncol=all_param_length+2)
    chain_index <- 1
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=all_param_length+2)
    probab <- posterior(t_pars, current_params, par_names, par_labels, startDays, endDays, buckets, microCeph, births, data_locals, incDat, peakTimes, allPriors, ...)
    ## Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,as.numeric(current_params),probab)
    colnames(tmp_table) <- chain_colnames
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    no_recorded <- 1
    sampno <- 2
    index2 <- 1
    ## Go through chain
    for (i in 1:(iterations+adaptive_period+burnin)){
        if(is.null(mvrPars)) {
            ## For each parameter (Gibbs)
            proposal <- current_params
            j <- non_fixed_params[index2]
            index2 <- index2 + 1
            if(index2 > length(non_fixed_params)) index2 <- 1
            proposal <- proposalfunction(proposal, lower_bounds, upper_bounds, steps,j)
        } else proposal <- mvr_proposal(current_params, non_fixed_params, scaledCovMat)
        names(proposal) <- names(current_params)
        
        ## Propose new parameters and calculate posterior
        if(!any(proposal[non_fixed_params] < lower_bounds[non_fixed_params] | proposal[non_fixed_params] > upper_bounds[non_fixed_params])){
            newprobab <- posterior(t_pars, proposal, par_names, par_labels, startDays, endDays, buckets, microCeph, births, data_locals, incDat, peakTimes, allPriors, ...)
            ## Calculate log difference in posteriors and accept/reject
            difflike <- newprobab - probab
            if ((!is.nan(difflike) & !is.infinite(newprobab)) & (runif(1) < exp(difflike))){
                current_params <- proposal
                probab <- newprobab
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else tempaccepted <- tempaccepted + 1
            }
        }
        if(is.null(mvrPars)) tempiter[j] <- tempiter[j] + 1
        else tempiter <- tempiter + 1
        
        ## If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        if(sampno %% thin ==0){
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-1)] <- current_params
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }
        
        ## Update step sizes based on acceptance rate
        ## Note that if opt_freq is 0, then no tuning will take place
        if(i < burnin & i%%opt_freq == 0){
            pcur <- tempaccepted/tempiter
            print(paste("Pcur: ",pcur[non_fixed_params],sep=""))
        }
        
        if(opt_freq != 0 & i > burnin & i <= (adaptive_period+burnin)){
            chain[chain_index,1] <- sampno
            chain[chain_index,2:(ncol(chain)-1)] <- current_params
            chain[chain_index,ncol(chain)] <- probab
            chain_index <- chain_index + 1
            if(chain_index - opt_freq > 0) {
                pcur <- tempaccepted/tempiter
                if(is.null(mvrPars)){
                    print(pcur[non_fixed_params])
                    tempaccepted <- tempiter <- reset
                    tmp_transform <- steps
                    for(x in non_fixed_params){
                        if(pcur[x] < popt - (TUNING_ERROR*popt) | pcur[x] > popt + (TUNING_ERROR*popt)){
                            tmp_transform[x] <- scaletuning(tmp_transform[x],popt,pcur[x])
                        }
                    }
                    print("Step sizes:")
                    print(tmp_transform[non_fixed_params])
                    steps <- tmp_transform
                } else {
                    print(paste("Pcur: ",pcur,sep=""))
                    scale <- scaletuning(scale,popt,pcur)
                    print(paste("New scale: ",scale,sep=""))
                    oldCov <- covMat
                    covMat <- cov(chain[,2:(ncol(chain)-1)])
                    covMat <- (1-w)*oldCov + w*covMat
                    scaledCovMat <- covMat*scale
                    tmpiter <- tmpaccepted <- 0
                }
                chain <- empty_chain
                chain_index <- 1
            }
        }
        
        if(no_recorded > save_block){
            print(i)
            write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            save_chain <- empty_save_chain
            no_recorded <- 1
        }
        sampno <- sampno + 1
    }
    
    ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
    ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
    ## rather than an array
    if(no_recorded > 2){
        write.table(save_chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
    }

    if(is.null(mvrPars)){
        covMat <- NULL
        scale <- NULL
    } else {
        steps <- NULL
    }
    return(list("file"=mcmc_chain_file,"covMat"=covMat,"scale"=scale, "steps"=steps))
}








#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @export
#' @useDynLib zikaProj
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}





#' MCMC diagnostic tests
#'
#' Runs some basic MCMC diagnostics on the given chain and saves a few plots. The diagnostics are:
#' \itemize{
#' \item{Gelman Diagnostics: }{Saves the Gelman plots at the given file location}
#' \item{Auto-correlation: }{Saves autocorrelation plots at the given file location}
#' }
#' @param mcmc_chains the entire MCMC chain to be tested
#' @param filename the full file path at which to save the diagnostics. _gelman.pdf will be appended, for example
#' @param param_table the parameter table
#' @param VERBOSE boolean flag for additional output. Defaulst to FALSE
#' @return returns any error messages raised during the tests
#' @export
mcmc_diagnostics <- function(mcmc_chains, filename, param_table,VERBOSE=FALSE){
    errors <- NULL
    final <- NULL
    if(length(mcmc_chains) > 1){
        gelman.error <- tryCatch({
            if(VERBOSE) print("Saving Gelman diagnostics")
            gelman.filename <- paste(filename,"_gelman.pdf",sep="")
            pdf(gelman.filename)
            gelman.plot(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            if(VERBOSE) print("Gelman diagnostics:")
            gelman.diag(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war1){
            if(VERBOSE) print(paste("A warnings occured in gelman diagnostics: ", war1))
            final <- war1
        }, error=function(err1){
            if(VERBOSE) print(paste("An error occured in gelman diagnostics: ", err1))
            final <- err1
        }, finally = {
            dev.off()
            errors <- c(errors, final)
        })
    }

    for(i in 1:length(mcmc_chains)){
        autocorr.error <- tryCatch({
            if(VERBOSE) print("Saving auto-correlation plot")
            autocorr.filename <- paste(filename,"_autocor_",i,".pdf",sep="")
            pdf(autocorr.filename)
            autocorr.plot(mcmc_chains[[i]][,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war2){
            if(VERBOSE) print(paste("A warnings occured in autocorr plot: ", war2))
            final <- war2
        }, error=function(err2){
            if(VERBOSE) print(paste("An error occured in autocorr plot diagnostics: ", err2))
            final <- err2
        }, finally = {
            dev.off()
            errors <- c(errors, final)            
        })
    }
    
    return(errors)
}


#' @export
#' @useDynLib zikaProj
make_optim_r0 <- function(t_pars, pars, incDat){
    function(x){
        pars["density"] <- x[1]
        pars["constSeed"] <- x[2]
        pars["incPropn"] <- x[3]
        pars["baselineInc"] <- x[4]
        y0s <- generate_y0s(pars["N_H"],pars["density"])
        
        y <- solveModelSimple(t_pars,y0s,pars)
        tmpY <- y[y[,"times"] >= min(incDat[,"startDay"]) & y[,"times"] <= max(incDat[,"endDay"]),]
        N_H <- average_buckets(rowSums(tmpY[,c("I_H","S_H","E_H","R_H")]), incDat[,"buckets"])
        inc <- average_buckets(tmpY[,"I_H"], incDat[,"buckets"])
        perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
        return(-incidence_likelihood(perCapInc, incDat[,"inc"],incDat[,"N_H"]))
    }
}

#' @export
#' @useDynLib zikaProj
make_optim_micro <- function(t_pars, pars, parNames,microDat){
    function(x){
        pars[names(pars) %in% parNames] <- x

        y0s <- generate_y0s(pars["N_H"],pars["density"])
        
        y <- solveModelSimple(t_pars,y0s,pars)
        probs <- generate_micro_curve(pars)
        probM <- generate_probM(y[,"I_M"], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)*pars["propn"]
        probM <- probM[which(y[,"times"] >= min(microDat$startDay) & y[,"times"] <= max(microDat$endDay))]
        probM <- average_buckets(probM, microDat$buckets)
        
        return(-likelihood_probM(microDat$microCeph, microDat$births, probM))
    }
}

#' @export
#' @useDynLib zikaProj
optimise <- function(state, parTab, t_pars, incDat, testDat){
    incPars <- c("density","constSeed","incPropn","baselineInc")
    microcephPars <- c("baselineProb","mean","var","c","p1","p2","p3","p4","p5","p6","p7","p8")
    checkLn <- NaN
    tmpDat <- incDat[incDat$local == state,]
    tmpTab <- parTab[parTab$local %in% c("all",state),]
    tmpMicro <- testDat[testDat$local == state,]
    indices <- parTab$local %in% c("all",state)
    while(is.nan(checkLn) | is.infinite(checkLn)){
        startPars <- generate_start_pars(tmpTab)
        if(!is.null(incDat)){
            lik_r0 <- make_optim_r0(t_pars, startPars, tmpDat)
            opti1 <- optim(startPars[names(startPars) %in% incPars],lik_r0,control=list("maxit"=5000))
            startPars[names(startPars) %in% incPars] <- opti1$par
        }
        lik_micro <- make_optim_micro(t_pars, startPars, startPars[names(startPars) %in% microcephPars], tmpMicro)
        opti2 <- optim(startPars[names(startPars) %in% microcephPars], lik_micro, control=list("maxit"=5000))
        startPars[names(startPars) %in% microcephPars] <- opti2$par
        startPars <- bounds_check(startPars, parTab)
        checkLn <- posterior(t_pars, startPars, parTab$names, parTab$local,testDat$startDay,testDat$endDay,testDat$buckets,testDat$microCeph,testDat$births,testDat$local,incDat)
    }
    tmpPars <- parTab$values
    tmpPars[indices] <- startPars
    names(tmpPars) <- parTab$names
    return(tmpPars)
}

#' @export
#' @useDynLib zikaProj
generate_start_pars <- function(parTab){
    startPars <- parTab$values
    names(startPars) <- parTab$names
    for(i in which(parTab$fixed==0)){
        startPars[i] <- runif(1,parTab[i,"lower_bounds"],parTab[i,"upper_bounds"])
    }
    return(startPars)
}

#' @export
#' @useDynLib zikaProj
bounds_check <- function(x, parTab){
    pars <- x
    for(i in 1:length(x)){
            if(x[i] < parTab[i,"lower_bounds"]) pars[i] <- parTab[i,"lower_bounds"]
            if(x[i] > parTab[i,"upper_bounds"]) pars[i] <- parTab[i,"upper_bounds"]
    }
    return(pars)
}
