
#' Posterior function for ODE model
#'
#' Given the time vector, initial conditions and ODE parameters, calculates the posterior value for a given data set. Note that no priors are used here. It seems like a lot of parameters, but these are basically the vectorised versions of the data frames provided by \code{\link{setupParTable}} and \code{\link{generate_microceph_dat}}
#' @param t_pars time vector over which to solve the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior <- function(t_pars, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat=NULL,actualPeakTimes=NULL,allPriors=NULL, ...){
    return(posterior_complex_buckets(t_pars,values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat, actualPeakTimes,allPriors))
   ## return(posterior_simple_buckets(t_pars,y0s,values, startDays, endDays, buckets, microCeph, births, incDat))
    ##if(version==2) return(posterior_complex(t_pars, y0s, values, startDays, endDays, buckets, microCeph, births, threshold, times))

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
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param allPriors defaults to FALSE. Arguments for parameter priors, if desired.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_simple_buckets <- function(t_pars, y0s, pars, startDays, endDays, buckets, microCeph, births, incDat = NULL, actualPeakTime=NULL, allPriors=NULL){
    lik <- 0
    
    y <- solveModelSimple(t_pars,y0s,pars)
    peakTime <- y[which.max(y[,"I_H"]),"times"]
    if(!is.null(incDat)){
        tmpY <- y[y[,"times"] >= min(incDat[,"startDay"]) & y[,"times"] <= max(incDat[,"endDay"]),]
        N_H <- average_buckets(rowSums(tmpY[,c("I_H","S_H","E_H","R_H")]), incDat[,"buckets"])
        inc <- average_buckets(tmpY[,"I_H"], incDat[,"buckets"])
        perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
        lik <- lik + incidence_likelihood(perCapInc, incDat[,"inc"],incDat[,"N_H"])
    }

    if(length(y) <= 1 && y=="Error") return(-Inf)
    
    probs <- generate_micro_curve(pars)
    probM <- generate_probM(y[,"I_M"], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)*pars["propn"]
    probM <- probM[which(y[,"times"] >= min(startDays) & y[,"times"] <= max(endDays))]
    probM <- average_buckets(probM, buckets)
   
    lik <- lik + likelihood_probM(microCeph, births, probM)

    if(!is.null(allPriors)) lik <- lik + allPriors(pars)
    if(!is.null(actualPeakTime)) lik <- lik + log(dunif(peakTime, actualPeakTime["start"],actualPeakTime["end"]))
    
    return(lik)
}


#' Posterior function for the simple SEIR model with bucketed data for multiple states
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param pars ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_complex_buckets <- function(t_pars, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat = NULL, peakTimes=NULL, allPriors=NULL){
    lik <- 0
    places <- unique(local)
    places <- places[places != "all"]
    
    for(place in places){
        indices <- data_locals == place | data_locals == "all"
        indices_pars <- local == place | local == "all"
        tmpMicro <- microCeph[indices]
        tmpBirths <- births[indices]
        tmpStart <- startDays[indices]
        tmpEnd <- endDays[indices]
        tmpBuckets <- buckets[indices]
        tmpPars <- values[indices_pars]
        names(tmpPars) <- names[indices_pars]
        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))

        tmpIncDat <- tmpPeaks <- NULL
        if(!is.null(incDat)) tmpIncDat <- incDat[incDat[,"local"] == place,]
        if(!is.null(peakTimes)){
            tmpPeaks <- as.numeric(peakTimes[peakTimes[,"local"] == place,c("start","end")])
            names(tmpPeaks) <- c("start","end")
        }
        lik <- lik + posterior_simple_buckets(t_pars, tmpY0s, tmpPars, tmpStart, tmpEnd, tmpBuckets, tmpMicro, tmpBirths, tmpIncDat, tmpPeaks,allPriors)
    }
    return(lik)
}


#' Incidence likelihood
#'
#' Given time-varying probabilities of observing Zika incidence and actual data, returns a single likelihood
#' @param perCapInc expected per capita incidence of Zika
#' @param dat two column matrix. First column is observed Zika cases, second is total possible zika cases
#' @return a single log likelihood
#' @export
#' @useDynLib zikaProj
incidence_likelihood <- function(perCapInc, inc, N_H){
   #return(likelihood_probM(inc,N_H,perCapInc))
    return(sum(dbinom(inc,N_H,perCapInc,log=1)))
}
