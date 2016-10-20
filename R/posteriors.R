
#' Posterior function for the simple SEIR model with bucketed data for multiple states
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param inc_startDays start day vector for incidence data
#' @param inc_endDays end day vector for incidence data
#' @param inc_locals vector of state names for incidence data
#' @param inc_buckets sampling window vector for incidence data
#' @param inc_ZIKV vector of ZIKV incidence for each time point
#' @param inc_NH vector of population size over time
#' @param peak_startDays vector of initial allowable day of incidence peak time
#' @param peak_endDays vector of end allowable day of incidence peak time
#' @param peak_locals vector of corresponding state name
#' @param unique_states vector of all unique states
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_complex_buckets <- function(ts, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals,inc_startDays=NULL,inc_endDays=NULL,inc_locals=NULL,inc_buckets=NULL,inc_ZIKV=NULL,inc_NH=NULL, peak_startDays=NULL, peak_endDays=NULL,peak_locals=NULL, unique_states, allPriors=NULL){
    lik <- 0

    ## For each state considered here
    for(place in unique_states){
        ## Get the indices for this state in the data and parameter table
        indices <- data_locals == place | data_locals == "all"
        indices_pars <- local == place | local == "all"

        ## Isolate these vectors
        tmpMicro <- microCeph[indices]
        tmpBirths <- births[indices]
        tmpStart <- startDays[indices]
        tmpEnd <- endDays[indices]
        tmpBuckets <- buckets[indices]
        tmpPars <- values[indices_pars]
        names(tmpPars) <- names[indices_pars]

        ## Generate starting y0 for the current parameter values
        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))

        ## Optional extra data
        tmpInc_ZIKV <- NULL
        tmpInc_NH<- NULL
        tmpInc_buckets <- NULL
        tmpInc_start <- NULL
        tmpInc_end <- NULL

        peak_start <- NULL
        peak_end <- NULL

        if(!is.null(inc_ZIKV)){
            indices_inc <- which(inc_locals == place)
            if(length(indices_inc) > 0){
                tmpInc_ZIKV<- inc_ZIKV[indices_inc]
                tmpInc_NH<- inc_NH[indices_inc]
                tmpInc_buckets <- inc_buckets[indices_inc]
                tmpInc_start <- inc_startDays[indices_inc]
                tmpInc_end <- inc_endDays[indices_inc]
            }
        } 
        if(!is.null(peak_startDays)){
            indices_peak <- which(peak_locals == place)
            if(length(indices_peak) > 0){
                peak_start <- peak_startDays[indices_peak]
                peak_end <- peak_endDays[indices_peak]
            }
        }
        
        lik <- lik + posterior_simple_buckets(ts,
                                              tmpY0s,
                                              tmpPars,
                                              tmpStart,
                                              tmpEnd,
                                              tmpBuckets,
                                              tmpMicro,
                                              tmpBirths,
                                              tmpInc_ZIKV,
                                              tmpInc_NH,
                                              tmpInc_buckets,
                                              tmpInc_start,
                                              tmpInc_end,
                                              peak_start,
                                              peak_end,
                                              allPriors)
    }
    return(lik)
}

#' Posterior function for the simple SEIR model with bucketed data for a single state
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param zikv vector of ZIKV incidence
#' @param nh vector of population sizes
#' @param inc_buckets vector of bucket size for incidence data
#' @param inc_start vector of inc start days
#' @param inc_end vector of inc end days
#' @param peak_start vector of peak time start days
#' @param peak_end vector of peak time end days
#' @param allPriors defaults to FALSE. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
posterior_simple_buckets <- function(ts, y0s, pars, startDays, endDays, buckets, microCeph, births, zikv=NULL,nh=NULL,inc_buckets=NULL,inc_start=NULL,inc_end=NULL,peak_start=NULL,peak_end=NULL, allPriors=NULL){
    lik <- 0

    ## Solve the ODE model with current parameter values
    y <- solveModelSimple_rlsoda(ts, y0s, pars,FALSE)
    y["I_M",][y["I_M",] < 0] <- 0
    
    ## Extract peak time. Need to add 1 as rlsoda does not return the first time point.
    peakTime <- y["time",which.max(y["I_H",])]
    
    if(!is.null(zikv)){
        tmpY <- y[,which(y["time",] >= min(inc_start) & y["time",] <= max(inc_end))]
        N_H <- average_buckets(colSums(tmpY[5:8,]), inc_buckets)
        
        inc <- diff(tmpY["incidence",])
        inc <- sum_buckets(inc,inc_buckets)
        inc[inc < 0] <- 0
        
        perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
        lik <- lik + incidence_likelihood(perCapInc, zikv,nh)
    }

    probs <- generate_micro_curve(pars)
    probM <- generate_probM(y["I_M",], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)*pars["propn"]
    probM <- probM[which(y["time",] >= min(startDays) & y["time",] <= max(endDays))]
    probM <- average_buckets(probM, buckets)
   
    lik <- lik + likelihood_probM(microCeph, births, probM)

    if(!is.null(allPriors)) lik <- lik + allPriors(pars)
    if(!is.null(peak_start)) lik <- lik + log(dunif(peakTime, peak_start,peak_end))
    
    return(lik)
}



#' Posterior function for the simple SEIR model with incidence only
#'
#' Given the time vector, initial conditions ODE parameters and deconstructed matrix of incidence data, calculates the posterior value for a given data set. 
#' @param ts time vector over which to solve the ODE model
#' @param values model parameters
#' @param names names of the model parameters
#' @param local the vector of state names corresponding to the parameter table
#' @param inc_startDays vector of all starting days for reporting buckets
#' @param inc_endDays vector of all end days for reporting buckets
#' @param inc_local vector of all state names corresponding to the incidence table
#' @param inc_buckets vector of all reporting bucket sizes
#' @param inc_ZIKV vector of reported ZIKV caes
#' @param inc_NH vector of all population sizes
#' @param unique_states vector of all unique states to be tested
#' @param allPriors defaults to FALSE. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
posterior_inc <- function(ts, values, names, local,inc_startDays,inc_endDays,inc_locals,
                          inc_buckets,inc_ZIKV,inc_NH, unique_states, allPriors=NULL){
    lik <- 0

    ## For each state considered here
    for(place in unique_states){
        indices_pars <- local == place | local == "all"
        indices_inc <- which(inc_locals == place)
        
        tmpPars <- values[indices_pars]
        names(tmpPars) <- names[indices_pars]

        ## Generate starting y0 for the current parameter values
        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))

        tmpInc_ZIKV<- inc_ZIKV[indices_inc]
        tmpInc_NH<- inc_NH[indices_inc]
        tmpInc_buckets <- inc_buckets[indices_inc]
        tmpInc_start <- inc_startDays[indices_inc]
        tmpInc_end <- inc_endDays[indices_inc]
        
        ## Solve the ODE model with current parameter values
        y <- solveModelSimple_rlsoda(ts, tmpY0s, tmpPars,FALSE)
        
        tmpY <- y[,which(y["time",] >= min(tmpInc_start) & y["time",] <= max(tmpInc_end))]
        N_H <- average_buckets(colSums(tmpY[5:8,]), tmpInc_buckets)
        
        inc <- diff(tmpY["incidence",])
        inc <- sum_buckets(inc,tmpInc_buckets)
        
        perCapInc <- (1-(1-(inc/N_H))*(1-tmpPars["baselineInc"]))*tmpPars["incPropn"]
        lik <- lik + incidence_likelihood(perCapInc, tmpInc_ZIKV,tmpInc_NH)
    }

    if(!is.null(allPriors)) lik <- lik + allPriors(values)
    return(lik)
}


#' Incidence likelihood
#'
#' Given time-varying probabilities of observing Zika incidence and actual data, returns a single likelihood
#' @param perCapInc expected per capita incidence of Zika
#' @param inc vector of actual incidence
#' @param N_H vector of population sizes (ie. per capita inc)
#' @return a single log likelihood
#' @export
#' @useDynLib zikaProj
incidence_likelihood <- function(perCapInc, inc, N_H){
    return(sum(dbinom(inc,N_H,perCapInc,log=1)))
}


#' Simplify posterior
#'
#' Simplifies the call to the posterior function using closures
#' @param ts time vector over which to solve the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param inc_startDays start day vector for incidence data
#' @param inc_endDays end day vector for incidence data
#' @param inc_locals vector of state names for incidence data
#' @param inc_buckets sampling window vector for incidence data
#' @param inc_ZIKV vector of ZIKV incidence for each time point
#' @param inc_NH vector of population size over time
#' @param peak_startDays vector of initial allowable day of incidence peak time
#' @param peak_endDays vector of end allowable day of incidence peak time
#' @param peak_locals vector of corresponding state name
#' @param unique_states vector of all unique state names
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
create_posterior <- function(ts, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, inc_startDays=NULL, inc_endDays=NULL,inc_locals=NULL,inc_buckets=NULL,inc_ZIKV=NULL,inc_NH=NULL, peak_startDays=NULL, peak_endDays=NULL,peak_locals=NULL,unique_states, allPriors=NULL){
    if(!is.null(microCeph)){
        f <- function(values){
            return(posterior_complex_buckets(ts, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, inc_startDays,inc_endDays,inc_locals,inc_buckets,inc_ZIKV,inc_NH, peak_startDays, peak_endDays,peak_locals, unique_states,allPriors))
        }
    } else {
        f <- function(values){
            return(posterior_inc(ts, values, names, local, inc_startDays,inc_endDays,inc_locals,inc_buckets,inc_ZIKV,inc_NH, unique_states,allPriors))
            
        }
    }
    return(f)
}

