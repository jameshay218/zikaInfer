
#' Posterior function for ODE model
#'
#' Given the time vector, initial conditions and ODE parameters, calculates the posterior value for a given data set. Note that no priors are used here. It seems like a lot of parameters, but these are basically the vectorised versions of the data frames provided by \code{\link{setupParTable}} and \code{\link{generate_microceph_dat}}
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param version which version of the model is being used. 1 is for the simple single state SEIR, 3 is for the multi-state version. 2 is for an outdated version of the model.
#' @param threshold Defaults to NULL. Used for an outdated version of the model.
#' @param times defaults to NULL. Used for an outdated version of the model.
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior <- function(t_pars, y0s, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, version = 1, threshold=NULL, times=NULL, incDat=NULL,allPriors=FALSE,actualPeakTimes=NULL){
    if(version==1) return(posterior_simple_buckets(t_pars,y0s,values, startDays, endDays, buckets, microCeph, births, incDat))
    if(version==2) return(posterior_complex(t_pars, y0s, values, startDays, endDays, buckets, microCeph, births, threshold, times))
    if(version==3) return(posterior_complex_buckets(t_pars,values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat, allPriors, actualPeakTimes))
}

#' Posterior function for the simple SEIR model
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set. Note that no priors are used here.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_simple <- function(t_pars, y0s, pars, microCeph, births){
    y <- solveModelSimple(t_pars,y0s,pars)
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
    lik <- likelihood_probM_all(microCeph,births,y[,"I_M"], probs, NH, b, pHM, bp, tstep)
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
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param allPriors defaults to FALSE. Arguments for parameter priors, if desired.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_simple_buckets <- function(t_pars, y0s, pars, startDays, endDays, buckets, microCeph, births, incDat = NULL, allPriors=FALSE, actualPeakTime=NULL){
    y <- solveModelSimple(t_pars,y0s,pars)
    
    peakTime <- y[which.max(y[,"I_H"]),"times"]
    
    y <- y[y[,"times"] >= min(startDays) & y[,"times"] <= max(endDays),]

    if(length(y) <= 1 && y=="Error") return(-Inf)
    
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

    probM <- average_buckets(probM, buckets)

    lik <- likelihood_probM(microCeph, births, probM)
    
    if(allPriors){
        lik <- lik + dnorm(pars["mean"], 12, 10,1)
    }
    if(!is.null(actualPeakTime)){
        lik <- lik + log(dunif(peakTime, actualPeakTime["start"],actualPeakTime["end"]))
        #lik <- lik <-  dnorm(peakTime, mean(c(actualPeakTime["start"],actualPeakTime["end"])), 10,1)
    }
    if(!is.null(incDat)) lik <- lik + incidence_likelihood(y[,"I_H"]/rowSums(y[,c("I_H","S_H","E_H","R_H")]),incDat)

    return(lik)
}


#' Posterior function for the simple SEIR model with bucketed data for multiple states
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
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
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_complex_buckets <- function(t_pars, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat = NULL, allPriors = NULL, peakTimes=NULL){
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
            
        
        tmpIncDat <- tmpPriors <- tmpPeaks <- NULL
        if(!is.null(incDat)) tmpIncDat <- incDat[incDat[,"local"] == place,]
        if(!is.null(priors)) tmpPriors <- allPriors[allPriors[,"local"] == place,]
        if(!is.null(peakTimes)){
            tmpPeaks <- as.numeric(peakTimes[peakTimes[,"local"] == place,c("start","end")])
            names(tmpPeaks) <- c("start","end")
        }
        lik <- lik + posterior_simple_buckets(t_pars, tmpY0s, tmpPars, tmpStart, tmpEnd, tmpBuckets, tmpMicro, tmpBirths, tmpIncDat, tmpPriors, tmpPeaks)
    }
    return(lik)
}

#' @export
posterior_SA <- function(values, fixed_values, t_pars, names_fixed, names_unfixed, local_fixed, local_unfixed, startDays, endDays, buckets, microCeph, births, data_locals, peakTimes=NULL){
    lik <- 0
    places <- unique(data_locals)
    places <- places[places != "all"]
    
    for(place in places){
        indices <- data_locals == place | data_locals == "all"
        indices_pars <- local_fixed== place | local_fixed== "all"
        tmpParsFixed <- fixed_values[indices_pars]
        names(tmpParsFixed) <- names_fixed[indices_pars]
        

        indices_pars_unfixed <- local_unfixed == place | local_unfixed == "all"

        tmpParsUnfixed <- values[indices_pars_unfixed]

        names(tmpParsUnfixed) <- names_unfixed[indices_pars_unfixed]

        tmpPars <- c(tmpParsFixed,tmpParsUnfixed)
        
        tmpMicro <- microCeph[indices]
        tmpBirths <- births[indices]
        tmpStart <- startDays[indices]
        tmpEnd <- endDays[indices]
        tmpBuckets <- buckets[indices]

        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))
            
        
        tmpIncDat <- tmpPriors <- tmpPeaks <- NULL
             if(!is.null(peakTimes)){
            tmpPeaks <- as.numeric(peakTimes[peakTimes[,"local"] == place,c("start","end")])
            names(tmpPeaks) <- c("start","end")
        }
        lik <- lik + posterior_simple_buckets(t_pars, tmpY0s, tmpPars, tmpStart, tmpEnd, tmpBuckets, tmpMicro, tmpBirths, tmpIncDat, tmpPriors, tmpPeaks)
    }
    return(-lik)
}


#' Oh wow
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param y0s initial conditions for the ODE model
#' @param pars ODE parameters
#' @param dat unnamed matrix over which to calculate likelihoods
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
posterior_apply <- function(place, t_pars, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals,incDat = NULL, allPriors = NULL, peakTimes=NULL){
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

    tmpIncDat <- tmpPriors <- tmpPeaks <- NULL
    if(!is.null(incDat)) tmpIncDat <- incDat[incDat[,"local"] == place,]
    if(!is.null(priors)) tmpPriors <- allPriors[allPriors[,"local"] == place,]
    if(!is.null(peakTimes)){
        tmpPeaks <- as.numeric(peakTimes[peakTimes[,"local"] == place,])
        print(tmpPeaks)
    }
    
    lik <- posterior_simple_buckets(t_pars, tmpY0s, tmpPars, tmpStart, tmpEnd, tmpBuckets, tmpMicro, tmpBirths, tmpIncDat, tmpPriors, tmpPeaks)
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
    y <- solveModel(t_pars,y0s,pars)
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
