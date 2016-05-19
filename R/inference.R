#' MCMC parameter setup
#'
#' Sets up the parameter table for use in \code{\link{run_metropolis_MCMC}}
#' @param pars the vector of parameters that would be used for the ODE model
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, boolean for log scale, initial step sizes, log proposal and whether or not the parameter should be fixed.
#' @export
setupParTable <- function(version=1, realDat=NULL){
    names <- c("sampFreq","sampPropn","mu_I","sd_I","mu_N","sd_N","probMicro","baselineProb","burnin","epiStart","L_M","D_EM","L_H","D_C","D_F","D_EH","D_IH","b","p_HM","p_MH","constSeed","mean","var","scale","tstep")    
    paramTable <- matrix(0, ncol=9, nrow=length(names))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local","use_log", "lower_bounds","upper_bounds","steps","log_proposal","fixed")
    paramTable[,"names"] <- names
    paramTable$names <- as.character(paramTable$names)
    

    paramTable[paramTable[,"names"]=="sampFreq",2:ncol(paramTable)] <- c(7,"all",0,0,30,0.1,0,1)
    paramTable[paramTable[,"names"]=="sampPropn",2:ncol(paramTable)] <- c(1,"all",0,0,1,0.1,0,1)
    paramTable[paramTable[,"names"]=="mu_I",2:ncol(paramTable)] <- c(30,"all",0,0,60,0.1,0,1)
    paramTable[paramTable[,"names"]=="sd_I",2:ncol(paramTable)] <- c(2,"all",0,0,10,0.1,0,1)
    paramTable[paramTable[,"names"]=="mu_N",2:ncol(paramTable)] <- c(30,"all",0,0,60,0.1,0,1)
    paramTable[paramTable[,"names"]=="sd_N",2:ncol(paramTable)] <- c(2,"all",0,0,10,0.1,0,1)
    paramTable[paramTable[,"names"]=="probMicro",2:ncol(paramTable)] <- c(0.1,"all",9,0,1,0.1,0,1)
    paramTable[paramTable[,"names"]=="baselineProb",2:ncol(paramTable)] <- c(0.002,"all",0,0,1,0.1,0,0)
    paramTable[paramTable[,"names"]=="burnin",2:ncol(paramTable)] <- c(0,"all",0,0,1000,0.1,0,1)
    paramTable[paramTable[,"names"]=="epiStart",2:ncol(paramTable)] <- c(0,"all",0,0,1000,0.1,0,0)
    paramTable[paramTable[,"names"]=="L_M",2:ncol(paramTable)] <- c(14,"all",0,0,100,0.1,0,1)
    paramTable[paramTable[,"names"]=="D_EM",2:ncol(paramTable)] <- c(4,"all",0,0,100,0.1,0,1)
    paramTable[paramTable[,"names"]=="L_H",2:ncol(paramTable)] <- c(365*70,"all",0,0,200*365,0.1,0,1)
    paramTable[paramTable[,"names"]=="D_C",2:ncol(paramTable)] <- c(365*18,"all",0,0,25*365,0.1,0,1)
    paramTable[paramTable[,"names"]=="D_F",2:ncol(paramTable)] <- c(0.75*365,"all",0,0,365,0.1,0,1)
    paramTable[paramTable[,"names"]=="D_EH",2:ncol(paramTable)] <- c(4,"all",0,0,100,0.1,0,1)
    paramTable[paramTable[,"names"]=="D_IH",2:ncol(paramTable)] <- c(5,"all",0,0,100,0.1,0,1)
    paramTable[paramTable[,"names"]=="b",2:ncol(paramTable)] <- c(0.25,"all",0,0,100,0.1,0,0)
    paramTable[paramTable[,"names"]=="p_HM",2:ncol(paramTable)] <- c(0.5,"all",0,0,1,0.1,0,1)
    paramTable[paramTable[,"names"]=="p_MH",2:ncol(paramTable)] <- c(0.5,"all",0,0,1,0.1,0,1)
    paramTable[paramTable[,"names"]=="constSeed",2:ncol(paramTable)] <- c(0,"all",0,0,100,0.1,0,1)
    paramTable[paramTable[,"names"]=="mean",2:ncol(paramTable)] <- c(12,"all",0,0,100,0.1,0,0)
    paramTable[paramTable[,"names"]=="var",2:ncol(paramTable)] <- c(10,"all",0,0,100,0.1,0,0)
    paramTable[paramTable[,"names"]=="scale",2:ncol(paramTable)] <- c(1,"all",0,0,100,0.1,0,0)
    paramTable[paramTable[,"names"]=="tstep",2:ncol(paramTable)] <- c(7,"all",0,0,100,0.1,0,1)
    
    paramTable <- paramTable[paramTable[,"names"] %in% names(setupParsLong(version)),]
    if(!is.null(realDat)) paramTable <- rbind(paramTable, setupStateParTable(realDat))

    paramTable[,c("values","use_log","lower_bounds","upper_bounds","steps","log_proposal","fixed")] <- lapply(paramTable[,c("values","use_log","lower_bounds","upper_bounds","steps","log_proposal","fixed")], FUN=as.numeric)
    
    return(paramTable)
    

}

#' Setup state specific parameter table
#'
#' Takes the entire data set and returns a useable format parameter table for the place specific parameters
#' @param stateDat the data frame of all data
#' @return a parameter table
#' @export
setupStateParTable <- function(stateDat){
    places <- as.character(unique(stateDat$local))
    paramTable <- matrix(0, ncol=9, nrow=4*length(places))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local","use_log", "lower_bounds","upper_bounds","steps","log_proposal","fixed")
    index <- 1
    for(place in places){
        tmpDat <- stateDat[stateDat$local==place,]
        paramTable[index,] <- c("L_H",tmpDat[1,"L_H"],place,0,0,200*365,0.1,0,1)
        index <- index + 1
        paramTable[index,] <- c("N_H", tmpDat[1,"N_H"], place, 0,0,100000000,0.1,0,1)
        index <- index + 1
        paramTable[index,] <- c("density",3,place,0,0,100,0.1,0,0)
        index <- index + 1
        paramTable[index,] <- c("epiStart",50,place,0,0,600,0.1,0,0)
        index <- index + 1
    }
    return(paramTable)
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

#' MCMC proposal function
#'
#' Proposal function for MCMC random walk, taking random steps of a given size. Random walk may be on a linear or log scale
#' @param param a vector of the parameters to be explored
#' @param param_table a matrix of flags and bounds that specify how the proposal should be treated. Crucial elements are the 3rd column (lower bound), 4th column (upper bound), 5th column )step size) and 6th column (flag for log proposal)
#' @param index numeric value for the index of the parameter to be moved from the param table and vector
#' @return the parameter vector after step
#' @export
#' @useDynLib zikaProj
proposalfunction <- function(param_table,index){
                                        # 4th index is upper bound, 3rd is lower
                                        # 1st and 2nd index used in other functions
    mn <- param_table[index,"lower_bounds"]
    mx <- param_table[index,"upper_bounds"]

    rtn <- param_table[index,"values"]
    
    x <- rtn[index]
    
    x <- toUnitScale(param_table[index,"values"],mn,mx)

                                        # 5th index is step size
    stp <- param_table[index,"steps"]

    rv <- runif(1)
    rv <- (rv-0.5)*stp
    x <- x + rv

                                        # Bouncing boundary condition
    if (x < 0) x <- -x
    if (x > 1) x <- 2-x

                                        # Cyclical boundary conditions
                                        #if (x < 0) x <- 1 + x	
                                        #if (x > 1) x <- x - 1
    
    if(x < 0 | x > 1) print("Stepped outside of unit scale. Something went wrong...")

    rtn[index] <- fromUnitScale(x,mn,mx)
    rtn
}

#' Multivariate proposal function
#'
#' Given the current parameters, bounds and a covariance matrix, returns a vector for a proposed jump from a multivariate normal distribution
#' @param current the vector of current parameter values
#' @param param_table a matrix of flags and bounds that specify how the proposal should be treated. Crucial elements are the 3rd column (lower bound), 4th column (upper bound), 5th column )step size) and 6th column (flag for log proposal)
#' @param covMat the 2D covariance matrix for all of the parameters
#' @return a parameter vector of a proposed move. Note that these may fall outside the allowable ranges.
#' @export
#' @useDynLib zikaProj
mvr_proposal <- function(param_table, covMat){
    proposed <- param_table$values
    fixed <- param_table$fixed==0
    proposed[fixed] <- proposed[fixed] + mvrnorm(n=1,mu=rep(0,length(proposed[fixed])),covMat[fixed,fixed])
                                        #proposed[fixed] <- proposed[fixed] + mvrnorm(n=1,mu=rep(0,length(proposed[fixed])),covMat)
                                        #proposed[fixed] <- mvrnorm(n=1,mu=proposed[fixed],covMat)
    return(proposed)
}

#' Posterior function for complex ODE model
#'
#' Given the time vector, initial conditions and ODE parameters, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @param threshold optional paramter. If not null, should be the threshold for microcephaly diagnosis
#' @param times optional parameter. If not null, then this should be a matrix of times over which data is recorded
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_complex <- function(t_pars, y0s, pars, dat, threshold=NULL, times=NULL){
    y <- solveModel(list(t_pars,y0s,pars))
    if(length(y) <= 1 && y=="Error") return(-Inf)
    sampFreq <- pars["sampFreq"]
    sampPropn <- pars["sampPropn"]
    mu_I <- pars["mu_I"]
    sd_I <- pars["sd_I"]
    mu_N <- pars["mu_N"]
    sd_N <- pars["sd_N"]
    probMicro <- pars["probMicro"]
    baselineProb <- pars["baselineProb"]

    if(!is.null(threshold)){
        if(!is.null(times)){
            alphas <- calculate_alphas_buckets(
                as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),
                probMicro,
                times)
            lik <- likelihood_threshold(
                as.matrix(unname(dat[,c("microCeph","births")])),
                unname(cbind(alphas,1-alphas)),
                c(mu_I,mu_N),
                c(sd_I,sd_N),
                threshold)
        } else {
            alphas <- calculate_alphas(
                as.matrix(unname(y[,c("I_F","S_F","E_F","R_F")])),
                probMicro,
                sampFreq)
            lik <- likelihood(
                as.matrix(unname(dat)),
                unname(cbind(alphas,1-alphas)),
                c(mu_I,mu_N),
                c(sd_I,sd_N)
            )
        }
    }
    else{
        if(!is.null(times)){
            alphas <- calculate_alphas_prob_buckets(
                as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),
                probMicro,
                baselineProb,
                times
            )
        } else {
            alphas <- calculate_alphas_prob_sampfreq(
                as.matrix(unname(y[,c("I_F","S_F","E_F","R_F")])),
                probMicro,
                baselineProb,
                sampFreq)
        }
        lik <- likelihood_prob(
            as.matrix(unname(dat[,c("microCeph","births")])),
            alphas
        )
        
    }
    return(lik)
}

#' Posterior function for ODE model
#'
#' Given the time vector, initial conditions and ODE parameters, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @param version optional parameter to specify which model we are solving. If 1, this is the simple SEIR model
#' @param threshold optional paramter. If not null, should be the threshold for microcephaly diagnosis
#' @param times optional parameter. If not null, then this should be a matrix of times over which data is recorded
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior <- function(t_pars, y0s, paramTable, dat, version = 1, threshold=NULL, times=NULL, incDat=NULL,allPriors=NULL,actualPeakTimes=NULL){
    tmpPars <- paramTable$value
    names(tmpPars) <- paramTable$names
    if(version==1) return(posterior_simple_buckets(t_pars,y0s,tmpPars,dat, incDat))
    if(version==2) return(posterior_complex(t_pars, y0s, tmpPars, dat, threshold, times))
    if(version==3) return(posterior_complex_buckets(t_pars,paramTable,dat,incDat,allPriors,actualPeakTimes))
}

#' Posterior function for the simple SEIR model
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_simple <- function(t_pars, y0s, pars, dat){
    y <- solveModelSimple(list(t_pars,y0s,pars))
    if(length(y) <= 1 && y=="Error"){
        print("Wow")
        return(-Inf)
    }
    NH <- sum(y0s[c("S_H","E_H","I_H","R_H")])
    b <- pars["b"]
    pHM <- pars["p_HM"]
    tstep <- pars["tstep"]
    bp <- pars["baselineProb"]
    shape <- pars["shape"]
    rate <- pars["rate"]
    scale <- pars["scale"]
    
    probs <- dgamma(0:39,shape,rate)*scale
    probs[probs > 1] <- 1
    lik <- likelihood_probM_all(dat[,1],dat[,2],y[,"I_M"], probs, NH, b, pHM, bp, tstep)
    return(lik)
}

#' Posterior function for the simple SEIR model with bucketed data
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_simple_buckets <- function(t_pars, y0s, pars, dat, incDat = NULL, allPriors=NULL, actualPeakTime=NULL){
    y <- solveModelSimple(list(t_pars,y0s,pars))
    peakTime <- y[which.max(y[,"I_H"]),"times"]
    
    y <- y[y[,"times"] >= min(dat[,"startDay"]) & y[,"times"] <= max(dat[,"endDay"]),]
    
    if(length(y) <= 1 && y=="Error"){
        print("Wow")
        return(-Inf)
    }
    NH <- sum(y0s[c("S_H","E_H","I_H","R_H")])
    b <- pars["b"]
    pMH <- pars["p_MH"]
    tstep <- pars["tstep"]
    bp <- pars["baselineProb"]
    scale <- pars["scale"]
    
    
    gammaMean <- pars["mean"]
    gammaVar <- pars["var"]
    
    rate <- gammaMean/gammaVar
    shape <- gammaMean*rate
    
    probs <- dgamma(0:39,shape,rate)*scale
    probs[probs > 1] <- 1
    probs <- rep(probs,each=tstep)

    probM <- generate_probM(y[,"I_M"],probs, NH, b, pMH, bp, 1)
    
    lik <- 0
    if(!is.null(dat)){
        buckets <- dat[,"buckets"]
        microDat <- dat[,"microCeph"]
        births <- dat[,"births"]
        probM <- average_buckets(probM, buckets)
        lik <- likelihood_probM(microDat, births, probM)
    }
    if(!is.null(allPriors)) lik <- lik + priors(allPriors)
    if(!is.null(actualPeakTime)) lik <- lik + dunif(peakTime, actualPeakTime["start"],actualPeakTime["end"],1)
    if(!is.null(incDat)) lik <- lik + incidence_likelihood(y[,"I_H"]/rowSums(y[,c("I_H","S_H","E_H","R_H")]),incDat)

    return(lik)
}


#' Posterior function for the simple SEIR model with bucketed data for multiple states
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @param 
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_complex_buckets <- function(t_pars, paramTable, dat, incDat = NULL, allPriors = NULL, peakTimes=NULL){
    places <- as.character(unique(dat$local))
    lik <- 0
    for(place in places){
        tmpDat <- dat[dat$local==place,c("startDay","endDay","buckets","microCeph","births")]
        tmpPars <- paramTable[paramTable$local=="all" | paramTable$local==place,"values"]
        tmpNames <- paramTable[paramTable$local=="all" | paramTable$local==place,"names"]
        names(tmpPars) <- tmpNames
        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))
        tmpIncDat <- tmpPriors <- tmpPeaks <- NULL
        if(!is.null(incDat)) tmpIncDat <- incDat[incDat$local==place,]
        if(!is.null(priors)) tmpPriors <- allPriors[allPriors$local==place,]
        if(!is.null(peakTimes)) tmpPeaks <- peakTimes[peakTimes$local==place,]
        lik <- lik + posterior_simple_buckets(t_pars, tmpY0s, tmpPars, tmpDat, tmpIncDat, tmpPriors, tmpPeaks)
    }
    return(lik)
}

#' All model priors
#'
#' Takes the vector of model parameters and returns a single value for the log prior probability
#' @param pars the vector of parameters
#' @return a single log prior value
#' @export
#' @useDynLib zikaProj
priors <- function(pars){
    meanPrior <- dnorm(pars["mean"],12,5,1)
    return(meanPrior)
}

incidence_likelihood <- function(perCapInc, dat){
    return(sum(dbinom(dat[,1],dat[,2],perCapInc,log=1)))    
}



    #' Adaptive Metropolis-within-Gibbs Random Walk Algorithm.
#'
#' The Adaptive Metropolis-within-Gibbs algorithm. Given a starting point and the necessary MCMC parameters as set out below, performs a random-walk of the posterior space to produce an MCMC chain that can be used to generate MCMC density and iteration plots. The algorithm undergoes an adaptive period, where it changes the step size of the random walk for each parameter to approach the desired acceptance rate, popt. After this, a burn in period is established, and the algorithm then uses \code{\link{proposalfunction}} to explore the parameter space, recording the value and posterior value at each step. The MCMC chain is saved in blocks as a .csv file at the location given by filename.
#' @param startvalue a vector of parameter values used as the starting point for the MCMC chain. MUST be valid parameters for the model function
#' @param iterations number of iterations to run the MCMC chain for. Note that each parameter is moved once independently for each iteration. Defaults to 1000
#' @param data the data against which the likelihood is calculated
#' @param ts vector of times to solve the ODE model over
#' @param y0s starting conditions for the ODE model
#' @param param_table a table of parameter data used for information such as bounds and prior function pointers.
#' @param popt the desired acceptance rate. Defaults to 0.44
#' @param opt_freq how frequently the acceptance rate is adapted. Defaults to 50
#' @param thin thinning value for the MCMC chain. Default is 1
#' @param burnin the length of the burn in period. Defaults to 100
#' @param adaptive_period length of the adaptive period. Defaults to 1
#' @param filename the full filepath at which the MCMC chain should be saved. "_chain.csv" will be appended to the end of this, so filename should have no file extensions
#' @param save_block the number of iterations that R will keep in memory before writing to .csv. Defaults to 500
#' @param _threshold boolean value that is true if using threshold data
#' @param VERBOSE boolean flag for additional output. Defaults to FALSE
#' @return the full file path at which the MCMC chain is saved as a .csv file
#' @export
#' @seealso \code{\link{posterior}}, \code{\link{proposalfunction}}
#' @useDynLib zikaProj
run_metropolis_MCMC <- function(iterations=1000,
                                data,
                                t_pars,
                                y0s,
                                N_H,
                                N_M,
                                version = 1,
                                param_table,
                                popt=0.44,
                                opt_freq=50,
                                thin=1,
                                burnin=100,
                                adaptive_period=1,
                                filename,
                                save_block = 500,
                                VERBOSE=FALSE,
                                threshold=NULL,
                                buckets=NULL,
                                mvrPars=NULL,
                                incDat=NULL,
                                allPriors=NULL,
                                peakTimes=NULL
                                ){
    TUNING_ERROR<- 0.1

    if(opt_freq ==0 && VERBOSE){ print("Not running adaptive MCMC - opt_freq set to 0")}
    else if(VERBOSE){ print("Adaptive MCMC - will adapt step size during specified burnin period")}

    ## Set up quicker tables
    ## Turns the data frame table into a matrix that can allow faster indexing
    ##    param_transform_table <- as.matrix(param_table[,c("names","value","use_log","lower_bounds","upper_bounds","steps","log_proposal")])
    
    ## Get those parameters which should be optimised
    non_fixed_params <- which(param_table$fixed==0)
    non_fixed_params_length <- length(non_fixed_params)
    all_param_length <- nrow(param_table)

    ## Setup MCMC chain file with correct column names
    mcmc_chain_file <- paste(filename,"_chain.csv",sep="")
    chain_colnames <- c("sampno",param_table$names,"r0","lnlike")
    
    ## Arrays to store acceptance rates
    if(is.null(mvrPars)) tempaccepted <- tempiter <- reset <- integer(all_param_length)
    else {
        tempaccepted <- tempiter <- 0
        covMat <- mvrPars[[1]]
        #covMat <- covMat[non_fixed_params,non_fixed_params]
        w <- mvrPars[[2]]
        scale <- mvrPars[[3]]
        #epsilon <- mvrPars[[3]]
        #scale <- (2.38)^2/non_fixed_params_length
        scaledCovMat <- covMat*scale
        #scaledCovMat <- scaledCovMat + epsilon*diag(nrow(scaledCovMat))
    }
    reset <- integer(all_param_length)
    reset[] <- 0

    # Create empty chain to store "save_block" iterations at a time
#    empty_chain <- chain <- matrix(nrow=save_block,ncol=all_param_length+2)
    empty_chain <- chain <- matrix(nrow=save_block,ncol=all_param_length+3)
    
    # Set starting value and params
    current_params <- paramTable$values
    
    # Create array to store values
    empty_values <- values <- sample <- numeric(save_block)

    probab <- posterior(t_pars, y0s, param_table, data,version, threshold, buckets, incDat,allPriors, peakTimes)
                                           # Set up initial csv file
    tmp_table <- array(dim=c(1,length(chain_colnames)))
    tmp_table <- as.data.frame(tmp_table)
    tmp_table[1,] <- c(1,as.numeric(param_table$values),r0.calc(param_table$values,N_H,N_M),probab)
    colnames(tmp_table) <- chain_colnames
    
    # Write starting conditions to file
    write.table(tmp_table,file=mcmc_chain_file,row.names=FALSE,col.names=TRUE,sep=",",append=FALSE)

    no_recorded <- 1
    sampno <- 2
                                            # Go through chain
    for (i in 1:(iterations+adaptive_period+burnin)){
        current_params <- param_table$values
        if(is.null(mvrPars)) {
            ## For each parameter (Gibbs)
            j <- sample(non_fixed_params,1)
            param_table$values <- proposalfunction(param_table,j)
            
        } else {
            param_table$values <- mvr_proposal(param_table,scaledCovMat)
        }
        ## Propose new parameters and calculate posterior
        if(!any(param_table[non_fixed_params,"values"]< param_table$lower_bounds[non_fixed_params] | proposal[non_fixed_params] > param_table$upper_bounds[non_fixed_params])){
            newprobab <- posterior(t_pars, y0s, param_table, data,version, threshold,buckets, incDat, allPriors,peakTimes, placePars)
            ## Calculate log difference in posteriors and accept/reject
            difflike <- newprobab - probab
            if ((!is.nan(difflike) & !is.infinite(newprobab)) & (runif(1) < exp(difflike) |  difflike > 0)){
                probab <- newprobab
                if(is.null(mvrPars)){
                    tempaccepted[j] <- tempaccepted[j] + 1
                } else {
                    tempaccepted <- tempaccepted + 1
                }
            } else { param_table$values <- current_params }
            
        }

        ##proposal[j] <- proposal_function(current_params[j],param_transform_table[j,"upper_bounds"],param_transform_table[j,"lower_bounds"],param_transform_table[j,"steps"])
     
        
        if(is.null(mvrPars)) tempiter[j] <- tempiter[j] + 1
        else tempiter <- tempiter + 1


                                        # If current iteration matches with recording frequency, store in the chain. If we are at the limit of the save block,
                                        # save this block of chain to file and reset chain
        if(sampno %% thin ==0){
            r0 <- r0.calc(current_params,N_H,N_M)
            chain[no_recorded,1] <- sampno
            chain[no_recorded,2:(ncol(chain)-2)] <- current_params
            chain[no_recorded,ncol(chain)-1] <- r0
            chain[no_recorded,ncol(chain)] <- probab
            no_recorded <- no_recorded + 1
          
        }
        sampno <- sampno + 1
                                        #        }
        
                                        # Update step sizes based on acceptance rate
                                        # Note that if opt_freq is 0, then no tuning will take place
        if(opt_freq != 0 & i <= adaptive_period & i%%opt_freq== 0) {
            pcur <- tempaccepted/tempiter
            if(is.null(mvrPars)){
                print(pcur[non_fixed_params])
                tempaccepted <- tempiter <- reset
                tmp_transform <- param_table[,"steps"]
                for(x in non_fixed_params){
                    if(pcur[x] < popt - (TUNING_ERROR*popt) | pcur[x] > popt + (TUNING_ERROR*popt)){
                        tmp_transform[x] <- scaletuning(tmp_transform[x],popt,pcur[x])
                                        #tmp_transform[x] <- scaletuning2(tmp_transform[x],popt,pcur[x])
                    }
                }
                print("Step sizes:")
                print(tmp_transform[non_fixed_params])
                param_table[,"steps"] <- tmp_transform
            } else {
                print(paste("Pcur: ",pcur,sep=""))
                scale <- scaletuning(scale,popt,pcur)
                print(paste("New scale: ",scale,sep=""))
                oldCov <- covMat
                covMat <- cov(chain[,2:(ncol(chain)-2)])
                #covMat <- cov(chain[1:i,non_fixed_params+1])

                covMat <- (1-w)*oldCov + w*covMat
                scaledCovMat <- covMat*scale
                #scaledCovMat <- covMat*(2.38)^2/non_fixed_params_length + epsilon*diag(nrow(covMat))
                tmpiter <- tmpaccepted <- 0
            }
        }
        if(no_recorded > save_block){
            print(i)
            write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,col.names=FALSE,row.names=FALSE,sep=",",append=TRUE)
            chain <- empty_chain
            no_recorded <- 1
          }
    }
        
        ## If there are some recorded values left that haven't been saved, then append these to the MCMC chain file. Note
        ## that due to the use of cbind, we have to check to make sure that (no_recorded-1) would not result in a single value
        ## rather than an array
        if(no_recorded > 2){
            write.table(chain[1:(no_recorded-1),],file=mcmc_chain_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
        }
        
        if(VERBOSE){
            print("Final step sizes:")
            print(param_table$step)
        }
        return(mcmc_chain_file)
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
