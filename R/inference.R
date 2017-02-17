#' Protect function
#'
#' Wrapper function to protect calls to the posterior function. If posterior does not compute correctly, returns -Inf.
#' @param f the function to be protected
#' @return the protected function
#' @export
protect <- function(f){
    function(...){
        tryCatch(f(...),error=function(e){
            message("caught error: ", e$message)
            -Inf
        })
    }
}

    
#' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{proposalfunction}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param data the microcephaly data against which the likelihood is calculated
#' @param ts vector of times to solve the ODE model over
#' @param param_table a table of parameter data used for information such as bounds and prior function pointers.
#' @param mcmcPars a named vector with parameters for the MCMC procedure. Iterations, popt, opt_freq, thin, burnin, adaptive_period and save_block.
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param mvrPars a list of parameters if using a multivariate proposal. Must contain an initial covariance matrix, weighting for adapting cov matrix, and an initial scaling parameter (0-1)
#' @param incDat optional data frame of incidence data if including incidence data in the likelihood function
#' @param peakTimes optional parameter - data frame of peak times for Zika incidence for each state
#' @param allPriors user function of prior for model parameters. Should take values, names and local from param_table
#' @return a list with: 1) full file path at which the MCMC chain is saved as a .csv file; 2) the last used covarianec matrix; 3) the last used scale size
#' @export
#' @seealso \code{\link{posterior_complex_buckets}}, \code{\link{proposalfunction}}
#' @useDynLib zikaProj
run_metropolis_MCMC <- function(data=NULL,
                                ts,
                                param_table,
                                mcmcPars,
                                filename,
                                mvrPars=NULL,
                                incDat=NULL,
                                peakTimes=NULL,
                                allPriors=NULL,
                                truePars=NULL
                                ){
                                        # MCMC par setup ---------------------------------------------------------- 
    ## Allowable error in scale tuning
    TUNING_ERROR <- 0.1
    OPT_TUNING  <- 0.2
    
    ## Extract MCMC parameters
    iterations <- mcmcPars["iterations"]
    popt <- mcmcPars["popt"]
    opt_freq<- mcmcPars["opt_freq"]
    thin <- mcmcPars["thin"]
    adaptive_period<- mcmcPars["adaptive_period"]
    save_block <- mcmcPars["save_block"]
    
    
                                        # Parameter par setup -------------------------------------
    ## Time vector for ODE model
    ## Extract parameter table into individual vectors for speed
    tmp_non_fixed_params <- non_fixed_params <- which(param_table$fixed==0)
    par_names <- param_table$names
    current_params <- start_values <- param_table$values
    par_labels <- param_table$local
    lower_bounds <- param_table$lower_bounds
    upper_bounds <- param_table$upper_bounds
    steps <- param_table$steps
    fixed <- param_table$fixed
    names(current_params) <- names(start_values) <- par_names
    unique_states <- unique(par_labels)
    unique_states <- unique_states[unique_states != "all"]
    all_states <- unique(param_table[param_table$fixed==0,"local"])
    ## Get some parameters useful for indexing
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- nrow(param_table)

  
    
    ## Arrays to store acceptance rates
    ## If univariate proposals
    if(is.null(mvrPars)){
        tempaccepted <- tempiter <- integer(all_param_length)
        reset <- integer(all_param_length)
        reset[] <- 0
    } else { # If multivariate proposals
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]][non_fixed_params,non_fixed_params]
        scale <- mvrPars[[3]]
    }
    
                                        # Data extraction ---------------------------------------------------------
    ## Extract data into vectors
    startDays <- NULL
    endDays <- NULL
    buckets <- NULL
    microCeph <- NULL
    births <- NULL
    data_locals <- NULL

    ## Microceph data
    if(!is.null(data)){
        startDays <- data[,"startDay"]
        endDays <- data[,"endDay"]
        buckets <- data[,"buckets"]
        microCeph <- data[,"microCeph"]
        births <- data[,"births"]
        data_locals <- data[,"local"]
    }

    ## Peak times
    peak_startDays <- NULL
    peak_endDays <- NULL
    peak_locals <- NULL
    if(!is.null(peakTimes)){
        peak_locals <- peakTimes$local
        peak_startDays <- peakTimes$start
        peak_endDays <- peakTimes$end
    }

    ## Incidence data
    inc_startDays <- NULL
    inc_endDays <- NULL
    inc_locals <- NULL
    inc_buckets <- NULL
    inc_ZIKV <- NULL
    inc_NH <- NULL
    if(!is.null(incDat)){
        inc_startDays <- incDat[,"startDay"]
        inc_endDays <- incDat[,"endDay"]
        inc_locals <- incDat[,"local"]
        inc_buckets <- incDat[,"buckets"]
        inc_ZIKV <- incDat[,"inc"]
        inc_NH <- incDat[,"N_H"]
    }
    
                                        # Posterior setup ---------------------------------------------------------
    ## Create posterior function with closures for neatness
 
    posterior_new <- create_posterior(ts, current_params, par_names, par_labels, 
                                       startDays, endDays, buckets, microCeph, births, 
                                       data_locals, inc_startDays,inc_endDays,inc_locals,
                                       inc_buckets,inc_ZIKV,inc_NH, peak_startDays, 
                                       peak_endDays,peak_locals, unique_states, allPriors, microDat, incDat)
    posterior_simp <- protect(posterior_new)
  
                                        # Chain setups ------------------------------------------------------------
    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",param_table$names,"lnlike")
    
    ## Create empty chain to store every iteration for the adaptive period
    opt_chain <- matrix(nrow=adaptive_period,ncol=non_fixed_params_length)
    chain_index <- 1
    
    ## Create empty chain to store "save_block" iterations at a time
    save_chain <- empty_save_chain <- matrix(nrow=save_block,ncol=all_param_length+2)
                                        # Initial conditions ------------------------------------------------------
    ## Initial likelihood
    probab <- posterior_simp(current_params)
    log_probab <- 0
    ## If available, find the true parameter posterior for comparison
    true_probab <- NULL
    if(!is.null(truePars)){
        true_probab <- posterior_simp(truePars)
        message(cat("True parameter posterior: ",true_probab,sep="\t"))
    }
    
    ## Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,as.numeric(current_params),probab)
    colnames(tmp_table) <- chain_colnames
    
    ## Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    ## Initial indexing parameters
    no_recorded <- 1
    sampno <- 2
    par_i <- 1
    
                                        # Main MCMC algorithm -----------------------------------------------------
    ## Go through chain
    for (i in 1:(iterations+adaptive_period)){
        ## If using univariate proposals
        if(is.null(mvrPars)) {
            ## For each parameter (Gibbs)
            j <- non_fixed_params[par_i]
            par_i <- par_i + 1
            if(par_i > length(non_fixed_params)) par_i <- 1
            proposal <- proposalfunction(current_params, lower_bounds, upper_bounds, steps,j)
            tempiter[j] <- tempiter[j] + 1
            ## If using multivariate proposals
        } else {
            proposal <- mvr_proposal(current_params, non_fixed_params, covMat)
            tempiter <- tempiter + 1
        }
        names(proposal) <- names(current_params)
        ## Propose new parameters and calculate posterior
        ## Check that all proposed parameters are in allowable range
        if(!any(
                proposal[non_fixed_params] < lower_bounds[non_fixed_params] |
                proposal[non_fixed_params] > upper_bounds[non_fixed_params]
            )
           ){
            ## Calculate new likelihood and find difference to old likelihood
            new_probab <- posterior_simp(proposal)
            log_prob <- min(new_probab-probab,0)
            if(!is.finite(log_prob)){
                #message("Not finite")
                #message(cat(proposal[c("density","constSeed")]," "))
            }
            ## Accept with probability 1 if better, or proportional to
            ## difference if not
            if(is.finite(log_prob) && log(runif(1)) < log_prob){
                current_params <- proposal
                probab <- new_probab
                
                ## Store acceptances
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else {
                    tempaccepted <- tempaccepted + 1
                }
            }
        }
        
        ## If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
        ## save this block of chain to file and reset chain
        if(i %% thin ==0){
            save_chain[no_recorded,1] <- sampno
            save_chain[no_recorded,2:(ncol(save_chain)-1)] <- current_params
            save_chain[no_recorded,ncol(save_chain)] <- probab
            no_recorded <- no_recorded + 1
        }
        
        ## If within adaptive period, need to do some adapting!

        if(i <= adaptive_period){
            ## Current acceptance rate
            pcur <- tempaccepted/tempiter
            ## Save each step
            opt_chain[chain_index,] <- current_params[non_fixed_params]
            
            ## If in an adaptive step
            if(chain_index %% opt_freq == 0){
                ## If using univariate proposals
                if(is.null(mvrPars)){
                    ## For each non fixed parameter, scale the step size
                    for(x in non_fixed_params) steps[x] <- scaletuning(steps[x],popt,pcur[x])
                    message(cat("Optimisation iteration: ", i,sep="\t"))
                    message(cat("Pcur: ", pcur[non_fixed_params],sep="\t"))
                    message(cat("Step sizes: ", steps[non_fixed_params],sep="\t"))
                    tempaccepted <- tempiter <- reset
                } else {       ## If using multivariate proposals
                    if(chain_index > OPT_TUNING*adaptive_period & chain_index < (0.9*adaptive_period)){
                        covMat <- scale*cov(opt_chain[1:chain_index,])
                        tempiter <- tempaccepted <- 0
                        message(cat("Optimisation iteration: ", i,sep="\t"))
                        ## Print acceptance rate
                        message(cat("Pcur: ", pcur,sep="\t"))
                        message(cat("Step size: ", scale,sep="\t"))
                    }
                    if(chain_index > (0.9)*adaptive_period){
                        scale <- scaletuning(scale, popt,pcur)
                        message(cat("Scale: ",scale,sep=""))
                    }
                }
            }
            chain_index <- chain_index + 1
        }
        if(i %% save_block == 0) message(cat("Current iteration: ", i, sep="\t"))
        if(no_recorded == save_block){
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
    return(list("file"=mcmc_chain_file,"covMat"=covMat,"scale"=scale, "steps"=steps,"truePosterior"=true_probab))
}






#' Scale step sizes
#'
#' Scales the given step size (between 0 and 1) based on the current acceptance rate to get closed to the desired acceptance rate
#' @param step the current step size
#' @param popt the desired acceptance rate
#' @param pcur the current acceptance rate
#' @return the scaled step size
#' @useDynLib zikaProj
scaletuning <- function(step, popt,pcur){
    if(pcur ==1) pcur <- 0.99
    if(pcur == 0) pcur <- 0.01
    step = (step*qnorm(popt/2))/qnorm(pcur/2)
    if(step > 1) step <- 1
    return(step)
}




#' Step size scaler
#'
#' Given current step scale, desired optimum step size etc, adapts the MCMC step size
#' @param step_scale the current step size
#' @param mc the number of MCMC iterations run
#' @param popt the desired acceptance rate
#' @param log_prob the current log probability
#' @param N_adapt the total number of adaptive steps to be taken
#' @export
rm_scale <- function(step_scale, mc, popt,log_prob, N_adapt)
{
    dd <- exp(log_prob)
    if( dd < -30 ){ dd <- 0 }
    dd <- min( dd, 1 )

    rm_temp <- ( dd - popt )/( (mc+1)/(0.01*N_adapt+1) )
    
    out <- step_scale*exp(rm_temp)
    
    out <- max( out, 0.02 )
    out <- min( out, 2)
    out
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
    proposed[fixed] <- MASS::mvrnorm(n=1,mu=proposed[fixed],Sigma=(5.6644/length(fixed))*covMat)
    return(proposed)
}

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


#' Generate allowable starting parameters
#'
#' Given a desirable SEIR peak time and time window, generates a set of R0 and t0 parameters that satisfy this peak time
#' @param peakTime the desired central peak time. Defaults to 927
#' @param peakTimeRange the desired allowable peak time range (ie. square prior)
#' @param stateNames the set of Brazilian states for which allowable parameters should be generated. Must match those in the microDatFile
#' @param microDatFile the location of the microcephaly data file which is needed to obtain N_H and L_H for each state. Note that this is only needed if no parTab is provided
#' @param parTab the parameter table as returned by \code{\link{setupParTable}}
#' @param allowableParsFile file location for the allowable parameters table, if it exists.
#' @return a table of allowable parameters with columns for t0, mosquito density (R0), corresponding state (as this will vary by N_H and life expectancy), and corresponding peak time
#' @export
generate_allowable_params <- function(peakTimings=927, peakTimeRange=60, stateNames,parTab=NULL,allowableParsFile="allowablePars.csv", R0max=6.5){
    if(!is.null(allowableParsFile) & file.exists(allowableParsFile)) allowablePars <- read.table(allowableParsFile)
    else {
        allowablePars <- NULL
        peakTimes <- matrix(nrow=100,ncol=100)
        for(local in stateNames){
            tmpTab <- parTab[parTab$local %in% c(local, "all"),]
            for(i in 1:nrow(peakTimes)){
                for(j in 1:ncol(peakTimes)){
                    pars <- tmpTab$values
                    names(pars) <- tmpTab$names
#                    print(pars)
                    pars["density"] <- j/10
                    pars["constSeed"] <- i*10
                    y0s <- generate_y0s(as.numeric(pars["N_H"]),as.numeric(pars["density"]))
                    t_pars <- seq(0,3003,by=1)
                    y <- solveModelSimple_rlsoda(t_pars, y0s,pars,TRUE)
                    peakTimes[i,j] <- y[which.max(diff(y[,"incidence"])),"time"]
                    R0 <- r0.calc(pars)
                    if(R0 > 1 & R0 < R0max & peakTimes[i,j] > (peakTimings[local] - peakTimeRange/2) & peakTimes[i,j] < (peakTimings[local] + peakTimeRange/2)){
                        allowablePars <- rbind(allowablePars,data.frame(i*10,j/10,local,peakTimes[i,j],R0))
                    }
                }
            }
        }
        allowabledPars <- allowablePars[complete.cases(allowablePars),]
        colnames(allowablePars) <- c("constSeed","density","local","peak","r0")
    }
    return(allowablePars)
}

