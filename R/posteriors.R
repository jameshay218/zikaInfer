#' Posterior function for the simple SEIR model with bucketed data for multiple locations
#'
#' Given the time vector, initial conditions ODE parameters and a matrix of microcephaly data, calculates the posterior value for a given data set.
#' @param ts time vector over which to solve the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the locations that are being considered in microcephaly data corresponding to parTab (ie. gets indices for each location from the parameter table)
#' @param startDays vector of the start day for each bucket of microcephaly data
#' @param endDays vector of the end day for each bucket of microcephaly data
#' @param buckets bucket sizes for microcephaly sampling periods
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of location names for the microcephaly data
#' @param inc_startDays start day vector for incidence data
#' @param inc_endDays end day vector for incidence data
#' @param inc_locals vector of location names for incidence data
#' @param inc_buckets sampling window vector for incidence data
#' @param inc_ZIKV vector of ZIKV incidence for each time point
#' @param inc_NH vector of population size over time
#' @param peak_startDays vector of initial allowable day of incidence peak time
#' @param peak_endDays vector of end allowable day of incidence peak time
#' @param peak_locals vector of corresponding location names
#' @param unique_locations vector of all unique locations
#' @param microceph_valid_days vector of all days included in the microcephaly data
#' @param microceph_valid_local vector of locations corresponding to the days in microceph_valid_days
#' @param inc_valid_days vector of all days included in the incidence data
#' @param inc_valid_local vector of locations corresponding to the days in inc_valid_days
#' @param solver specify which ODE solver to use, between "rlsoda" and "lsoda"
#' @return a single value for the posterior
#' @export
posterior_complex_buckets <- function(ts, values, names, local, ## Inputs related to model parameters
                                      startDays, endDays, buckets, microCeph, births, data_locals, ## Inputs related to microcephaly data 
                                      inc_startDays=NULL,inc_endDays=NULL,inc_locals=NULL,inc_buckets=NULL,inc_ZIKV=NULL,inc_NH=NULL, ## Inputs related to incidence data
                                      peak_startDays=NULL, peak_endDays=NULL,peak_locals=NULL, ## Inputs related to peak time priors
                                      unique_locations,
                                      microceph_valid_days=NULL, microceph_valid_local=NULL, inc_valid_days=NULL, inc_valid_local=NULL,## Inputs related to indexing data by specific days
                                      solver="rlsoda"
                                      ){
    lik <- 0
   
    ## For each location considered here
    for(place in unique_locations){
        ## Get the indices for this location in the data and parameter tables
        indices <- data_locals == place | data_locals == "all"
        indices_pars <- local == place | local == "all"

        ## Isolate these vectors
        tmpMicro <- microCeph[indices]
        tmpBirths <- births[indices]
        tmpStart <- startDays[indices]
        tmpEnd <- endDays[indices]
        tmpBuckets <- buckets[indices]
        tmpPars <- values[indices_pars]
        tmpValidDaysMicro <- NULL
        if(!is.null(microceph_valid_days)) tmpValidDaysMicro <- microceph_valid_days[microceph_valid_local == place]

        ## get model parameters related to this location
        names(tmpPars) <- names[indices_pars]

        ## Generate starting y0 for the current parameter values
        tmpY0s <- generate_y0s(as.numeric(tmpPars["N_H"]),as.numeric(tmpPars["density"]))

        ## Optional extra ZIKV incidence data
        tmpInc_ZIKV <- NULL
        tmpInc_NH<- NULL
        tmpInc_buckets <- NULL
        tmpInc_start <- NULL
        tmpInc_end <- NULL
        tmpValidDaysInc <- NULL

        ## Optional ZIKV epidemic peak times for prior
        peak_start <- NULL
        peak_end <- NULL

        ## If using incidence data, include
        if(!is.null(inc_ZIKV)){
            indices_inc <- which(inc_locals == place)
            if(length(indices_inc) > 0){
                tmpInc_ZIKV<- inc_ZIKV[indices_inc]
                tmpInc_NH<- inc_NH[indices_inc]
                tmpInc_buckets <- inc_buckets[indices_inc]
                tmpInc_start <- inc_startDays[indices_inc]
                tmpInc_end <- inc_endDays[indices_inc]
                if(!is.null(inc_valid_days)) tmpValidDaysInc <- inc_valid_days[inc_valid_local == place]
            }
        }

        ## If using peak time prior, include
        if(!is.null(peak_startDays)){
            indices_peak <- which(peak_locals == place)
            if(length(indices_peak) > 0){
                peak_start <- peak_startDays[indices_peak]
                peak_end <- peak_endDays[indices_peak]
            }
        }
        ## Calculate posterior probability for this location
        tmpLik <- posterior_simple_buckets(ts,tmpY0s,tmpPars,tmpStart,tmpEnd,tmpBuckets,tmpMicro,
                                           tmpBirths,tmpInc_ZIKV,tmpInc_NH,tmpInc_buckets,tmpInc_start,
                                           tmpInc_end,peak_start,peak_end, tmpValidDaysMicro,tmpValidDaysInc,
                                           solver)

        ## Add to total likelihood, but multiply by weighting given to this location
        lik <- lik + tmpLik*tmpPars["location_weight"]
    }
    return(lik)
}



#' Posterior function for forecast
#'
#' Using the microcephaly risk model, produces a forecast of microcephaly affected births for all time periods and calculates a likelihood of observing given microcephaly cases
#' @param pars model parameters
#' @param startDays vector of all starting days for reporting buckets
#' @param endDays vector of all end days for reporting buckets
#' @param buckets vector of all reporting bucket sizes
#' @param microCeph vector of reported microcephaly cases
#' @param births vector of birth numbers
#' @param zikv vector of ZIKV incidence
#' @param nh vector of population sizes over time
#' @param inc_buckets vector of reporting windows for ZIKV incidence
#' @param inc_start vector of all start days for ZIKV reporting buckets
#' @param inc_end vector of all end days for reporting ZIKV buckets
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @return a single value for the posterior
posterior_known_inc <- function(pars, startDays, endDays,
                                buckets, microCeph, births,
                                zikv, nh, inc_buckets,
                                inc_start, inc_end
                                ){
    ## Times at which reporting rates chance
    switch_time_i <- pars["switch_time_i"]
    switch_time_m <- pars["switch_time_m"]
    
    ## Generate microcephaly curve for these parameters
    probs <- generate_micro_curve(pars)
    
    ## Use incidence data to get daily incidence
    ## Get incidence before and after the switch time - different reporting rates
    inc_1 <- zikv[which(inc_start < switch_time_i)]/pars["incPropn"]
    inc_2 <- zikv[inc_start >= switch_time_i]/pars["incPropn2"]
    inc <- c(inc_1,inc_2)

    inc <- rep(inc/nh/inc_buckets, inc_buckets)

    ## Generate probability of observing a microcephaly case on a given day
    probM <- generate_probM_aux(inc, probs, pars["baselineProb"])
    probM[startDays < switch_time_m] <- probM[startDays < switch_time_m]*pars["propn"]
    probM[startDays >= switch_time_m] <- probM[startDays >= switch_time_m]*pars["propn2"]
    
    births[startDays >= switch_time_m] <- as.integer(births[startDays >= switch_time_m]*(1-pars["birth_reduction"]))

    probM <- average_buckets(probM,buckets)
    lik <- likelihood_probM(microCeph,births,probM)

    return(lik)
}


#' Posterior function for the simple SEIR model with bucketed data for a single location
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
#' @param valid_days_micro a vector of days for which data were collected, allowing non-contiguous comparisons. Microcephaly
#' @param valid_days_inc a vector of days for which data were collected, allowing non-contiguous comparisons. Incidence
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @return a single value for the posterior
#' @export
posterior_simple_buckets <- function(ts, y0s, pars,
                                     startDays, endDays, buckets, microCeph, births,
                                     zikv=NULL,nh=NULL,inc_buckets=NULL,inc_start=NULL,inc_end=NULL,
                                     peak_start=NULL,peak_end=NULL,
                                     valid_days_micro=NULL,valid_days_inc=NULL, solver="rlsoda"){
    lik <- 0

    ## Solve the ODE model with current parameter values
    y <- solveSEIRModel(ts, y0s, pars,solver)

    ## Make sure we haven't created neglibly small negative numbers
    y["I_M",][y["I_M",] < 0] <- 0
    
    ## Extract ZIKV epidemic peak time
    peakTime <- y["time",which.max(diff(y["incidence",]))]

    ## If including ZIV incidence in posterior
    if(!is.null(zikv)){
        inc <- diff(y["incidence",])
        inc[inc < 0] <- 0
        N_H <- colSums(y[5:8,])
        
        if(!is.null(valid_days_inc)){
            inc <- inc[y["time",] %in% valid_days_inc]
            N_H <- N_H[y["time",] %in% valid_days_inc]
        } else {
            inc <- inc[which(y["time",]  >= min(inc_start) & y["time",] <= max(inc_end))]
            N_H <- N_H[which(y["time",]  >= min(inc_start) & y["time",] <= max(inc_end))]
        }

        ## Bucket data by sampling windows
        inc <- sum_buckets(inc,inc_buckets)
        N_H <- average_buckets(N_H, inc_buckets)

        ## Modify with baseline reporting rate and reporting accuracy
        perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]

        ## Using normal or binomial error?
        if(is.na(pars["inc_sd"])){
            lik <- lik + pars["inc_weight"]*incidence_likelihood(perCapInc, zikv,nh)
        } else {
            lik <- lik + pars["inc_weight"]*incidence_likelihood_norm(perCapInc, zikv,nh,pars["inc_sd"])
        }
        
    }

    ## Generate the microcephaly risk curve
    probs <- generate_micro_curve(pars)

    ## Expected proportion of microcephaly affected births
    probM <- generate_probM(y["I_M",], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)*pars["propn"]

    if(is.null(valid_days_micro)){
        probM <- probM[which(y["time",] >= min(startDays) & y["time",] <= max(endDays))]
    } else {
        probM <- probM[y["time",] %in% valid_days_micro]
    }
    probM <- average_buckets(probM, buckets)

    if(is.na(pars["micro_sd"])) {
        lik <- lik + (1-pars["inc_weight"])*likelihood_probM(microCeph, births, probM)
    } else {
        lik <- lik + (1-pars["inc_weight"])*likelihood_probM_norm(microCeph, births, probM, unname(pars["micro_sd"]))
    }

    ## If using a prior on epidemic peak time, include here

    if(!is.null(peak_start)){
        lik <- lik + dunif(peakTime, peak_start,peak_end,1)
    }

    return(unname(lik))
}



#' Posterior function for the simple SEIR model with incidence only
#'
#' Given the time vector, ODE parameters and deconstructed matrix of incidence data, calculates the posterior probability for a given data set. 
#' @param ts time vector over which to solve the ODE model
#' @param values model parameters
#' @param names names of the model parameters
#' @param local the vector of location names corresponding to the parameter table
#' @param inc_startDays vector of all starting days for reporting buckets
#' @param inc_endDays vector of all end days for reporting buckets
#' @param inc_local vector of all location names corresponding to the incidence table
#' @param inc_buckets vector of all reporting bucket sizes
#' @param inc_ZIKV vector of reported ZIKV case
#' @param inc_NH vector of all population sizes
#' @param unique_locations vector of all unique locations to be tested
#' @param inc_valid_days vector of all days included in the incidence data
#' @param inc_valid_local vector of locals corresponding to the days in inc_valid_days
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @return a single value for the posterior
#' @export
posterior_inc <- function(ts, values, names, local,inc_startDays,inc_endDays,inc_locals,
                          inc_buckets,inc_ZIKV,inc_NH, unique_locations,
                          inc_valid_days=NULL,inc_valid_local=NULL,solver="rlsoda"){
    lik <- 0

    ## For each location considered here
    for(place in unique_locations){
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
        tmpValidDaysInc <- NULL
        if(!is.null(inc_valid_days)) tmpValidDaysInc <- inc_valid_days[inc_valid_local == local]
        
        ## Solve the ODE model with current parameter values
        y <- solveSEIRModel(ts, tmpY0s, tmpPars,solver)

        inc <- diff(y["incidence",])
        inc[inc < 0] <- 0
        N_H <- colSums(y[5:8,])
        
        if(!is.null(tmpValidDaysInc)){
            inc <- inc[y["time",] %in% valid_days_inc]
            N_H <- N_H[y["time",] %in% valid_days_inc]
        } else {
            inc <- inc[which(y["time",]  >= min(inc_start) & y["time",] <= max(inc_end))]
            N_H <- N_H[which(y["time",]  >= min(inc_start) & y["time",] <= max(inc_end))]
        }
        
        inc <- sum_buckets(inc,inc_buckets)
        N_H <- average_buckets(N_H, inc_buckets)
        
        
        perCapInc <- (1-(1-(inc/N_H))*(1-tmpPars["baselineInc"]))*tmpPars["incPropn"]
        if(is.na(tmpPars["inc_sd"])){
            lik <- lik + incidence_likelihood(perCapInc, tmpInc_ZIKV,tmpInc_NH)
        } else {
            lik <- lik + incidence_likelihood_norm(perCapInc, tmpInc_ZIKV,tmpInc_NH,tmpPars["inc_sd"])
        }
    }

    return(lik)
}

#' Posterior function for forecasting analysis
#'
#' Using the microcephaly risk model, produces a forecast of microcephaly affected births for all time periods and calculates a likelihood of observing given microcephaly cases
#' NOTE that baselineProb is on a log scale in this analysis
#' @param pars model parameters
#' @param startDays vector of all starting days for reporting buckets
#' @param endDays vector of all end days for reporting buckets
#' @param buckets vector of all reporting bucket sizes
#' @param microCeph vector of reported microcephaly cases
#' @param births vector of birth numbers
#' @param zikv vector of ZIKV incidence
#' @param nh vector of population sizes over time
#' @param inc_buckets vector of reporting windows for ZIKV incidence
#' @param inc_start vector of all start days for ZIKV reporting buckets
#' @param inc_end vector of all end days for reporting ZIKV buckets
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @return a single value for the posterior
#' @export
posterior_known_inc_seir <- function(pars, startDays, endDays,
                                     buckets, microCeph, births,
                                     zikv, nh, inc_buckets,
                                     inc_start, inc_end,
                                     solver="rlsoda"){
    lik <- 0

    ## Solve model from start to end of incidence reporting window
    ts <- seq(min(inc_start), max(inc_end)-1,by=1)
    
    ## Times after which reporting/birth behaviour changes
    switch_time_i <- pars["switch_time_i"]
    switch_time_m <- pars["switch_time_m"]
    switch_time_behaviour <- pars["switch_time_behaviour"]
    
    ## If this parameter is set to 1, then we're assuming that ZIKV reporting did not change
    check_par <- pars["zikv_reporting_change"]
    
    ## Generate microcephaly curve for these parameters
    probs <- generate_micro_curve(pars)

    ## Solve SEIR model
    y0s <- generate_y0s(as.numeric(pars["N_H"]),as.numeric(pars["density"]))
    y <- solveSEIRModel(seq(0,3003,by=1), y0s, pars,solver)

    ## Calculate SEIR generated incidence for before switch time
    calc_inc <- diff(y["incidence",])
    calc_inc[calc_inc < 0] <- 0
    N_H <- floor(colSums(y[5:8,]))

    calc_inc <- calc_inc[which(y["time",] >= min(inc_start) & y["time",] < switch_time_i)]
    
    ## Population size before switch time
    N_H <- N_H[which(y["time",] >= min(inc_start) & y["time",] < switch_time_i)]

    ## Get average buckets for this
    calc_inc <- sum_buckets(calc_inc,inc_buckets[which(inc_start < switch_time_i)])
    N_H <- average_buckets(N_H,inc_buckets[which(inc_start < switch_time_i)])
    
    ## Calculate per capita incidence from SEIR model
    perCapInc <- (1-(1-(calc_inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]

    ## Get subset of incidence for these times
    inc_1 <- zikv[which(inc_start < switch_time_i)]

    ## Calculate likelihood of SEIR model parameters given incidence up to this time
    lik <- lik + pars["inc_weight"]*sum(dbinom(x=inc_1,size=N_H,prob=perCapInc,log=TRUE))

    ## Convert to true underlying incidence
    inc_1 <- inc_1/pars["incPropn"]

    ## Get incidence after switch time with new reporting rate
    ## Otherwise use old reported rate
    if(check_par == 1){
        inc_2 <- zikv[which(inc_start >= switch_time_i)]/pars["incPropn2"]
    } else {
        inc_2 <- zikv[which(inc_start >= switch_time_i)]/pars["incPropn"]
    }
    inc <- c(inc_1,inc_2)
    
    ## Convert to daily incidence
    inc <- rep(inc/nh/inc_buckets, inc_buckets)

    ## Generate probability of observing a microcephaly case on a given day
    ## Note that this is on a log scale here
    bp <- exp(pars["baselineProb"])
    
    tmp_buckets <- buckets[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    tmp_births <- births[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    tmp_microCeph <- microCeph[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    ## Getting non-aborted births, (1-a)p_m(t)
    probM_a <- generate_probM_forecast(inc, probs, bp,
                                       pars["abortion_rate"], pars["birth_reduction"],
                                       which(ts == pars["switch_time_behaviour"]), pars["abortion_cutoff"], FALSE)
    ## Getting aborted births, ap_m(t)
    probM_b <- generate_probM_forecast(inc, probs, bp,
                                       pars["abortion_rate"], pars["birth_reduction"],
                                       which(ts == pars["switch_time_behaviour"]), pars["abortion_cutoff"], TRUE)

    ## So probability of becoming a microcephaly case and being observed is
    ## the proportion of births that become microcephaly cases of those that
    ## were not aborted microcephaly cases
    probM <- probM_a/(1-probM_b)
    
    probM[which(ts < switch_time_m)] <- probM[which(ts < switch_time_m)]*pars["propn"]
    probM[which(ts >= switch_time_m)] <- probM[which(ts >= switch_time_m)]*pars["propn2"]

    probM <- average_buckets(probM,tmp_buckets)
    
    ## Births after prediction time are not realised births - it may be that some of these
    ## were avoided at time t-40 from switch_time_behaviour
    ## Comment this out if we don't want to use avoided births parameter
    prediction_time <- pars["predicted_births"]
    
    ## Births after a certain time are predicted and not actual births
    predicted_births <- tmp_births[which(startDays >= prediction_time)]
    
    ## Can calculate the number of aborted births from microcephaly measurements
    ## and estimated microcephaly risk
    probM_b[which(ts < switch_time_behaviour)] <- 0

    probM_abortions <- probM_b/probM_a
    probM_abortions[is.nan(probM_abortions)] <- 0
    probM_abortions <- average_buckets(probM_abortions,tmp_buckets)
    aborted_births <- tmp_microCeph*probM_abortions
    
    ## Aborted births for the predicted birth time  
    aborted_births <- aborted_births[which(startDays >= prediction_time)]
    
    ## From this, we can infer the true number of births that happened assuming that some births
    ## were aborted and a proportion of the forecasted births did no occur
    inferred_births <- (1-pars["avoided_births"])*predicted_births - aborted_births
    
    ## The actual number of births after the prediction time are inferred
    births2 <- tmp_births
    births2[which(startDays >= prediction_time)] <- as.integer(inferred_births)
    lik <- lik + (1-pars["inc_weight"])*likelihood_probM(tmp_microCeph,births2,probM)

    return(lik)
}


#' Incidence likelihood (binomial)
#'
#' Given time-varying probabilities of observing Zika incidence and actual data, returns a single likelihood assuming binomially distributed observations
#' @param perCapInc expected per capita incidence of Zika
#' @param inc vector of actual incidence
#' @param N_H vector of population sizes (ie. per capita inc)
#' @return a single log likelihood
#' @export
incidence_likelihood <- function(perCapInc, inc, N_H){
    return(sum(dbinom(x=inc,size=N_H,prob=perCapInc,log=1)))
}

#' Incidence likelihood (normal)
#'
#' Given time-varying probabilities of observing Zika incidence and actual data, returns a single likelihood assuming normally distributed observations
#' @param perCapInc expected per capita incidence of Zika
#' @param inc vector of actual incidence
#' @param N_H vector of population sizes (ie. per capita inc)
#' @param lik_sd standard deviation of the observation distribution
#' @return a single log likelihood
#' @export
incidence_likelihood_norm <- function(perCapInc, inc, N_H, lik_sd){
    return(sum(dnorm(inc, perCapInc*N_H, lik_sd, 1)))
}

#' Posterior function for CRS prediction
#'
#' Using the microcephaly risk model, produces a forecast of CRS affected births for all time periods and calculates a likelihood of observing given CRS cases. Uses reported incidence as per capita risk (after dividing my reporting rate).
#' @export
posterior_CRS <- function(pars, startDays, endDays,
                                buckets, microCeph, births,
                                inc, nh, inc_buckets,
                                inc_start, inc_end){
    switch_time <- pars["switch_time"]
     ## Generate microcephaly curve for these parameters
    probs <- generate_micro_curve(pars)
    ## Use incidence data to get daily incidence
    ## Get incidence before and after the switch time - different reporting rates
    inc <- inc/pars["incPropn"]
    inc <- rep(inc/nh/inc_buckets, inc_buckets)
    
    ## Generate probability of observing a microcephaly case on a given day
    probM <- generate_probM_aux(inc, probs, pars["baselineProb"])
    probM[startDays < switch_time] <- probM[startDays < switch_time]*pars["propn"]
    probM[startDays >= switch_time] <- probM[startDays >= switch_time]*pars["propn2"]
   
    probM <- average_buckets(probM,buckets)
    lik <- likelihood_probM(microCeph,births,probM)

    return(lik)
}


#' Simplify posterior call
#'
#' Simplifies the call to the posterior function using closures. This function is set up to work with the lazymcmc package, found at github.com/jameshay218/lazymcmc
#' @param parTab the parameter table. See \code{\link{exampleParTab}}
#' @param data the microcephaly data. See \code{\link{exampleMicroDat}}
#' @param PRIOR_FUNC pointer to a function to calculate prior values. Can take "pars" and "..." as arguments, where "pars" is the model parameters.
#' @param incDat should also include incidence data, though this can be left as NULL. See \code{\link{exampleIncDat}}.
#' @param peakTimes include data on epidemic peak time ranges. Can be left as NULL. See \code{\link{examplePeakTimes}}
#' @param ts time vector over which to solve the ODE model
#' @param version usually just leave this as "binomial", indicating assumed binomially distributed errors. I've added a "forecast" version which expects parameters to do with the lack of second wave.
#' @return a single value for the posterior
#' @export
create_posterior <- function(parTab, data, PRIOR_FUNC, ...){
    args <- list(...)
    exists <- sapply(c("incDat","peakTimes","ts"), function(x) x %in% names(args))
    if(!any(exists)) stop("Error - insufficient arguments passed to posterior function. Should have: incDat, peakTimes, and ts. incDat and peakTimes can be NULL")

    peakTimes <- args[["peakTimes"]]
    incDat <- args[["incDat"]]
    ts <- args[["ts"]]
    version <- "binomial"
    if("version" %in% names(args)) version <- args[["version"]]
    solver <- "rlsoda"
    if("solver" %in% names(args)) version <- args[["solver"]]
    
    ## Extract main data into vectors
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


    ## Extract peak times
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

    ## Parameter data
    names <- parTab$names
    local <- parTab$local
    unique_locations <- unique(local)
    unique_locations <- unique_locations[unique_locations != "all"]
    
    ## I have included some code to explicitly extract all days included within the data sampling periods.
    ## this allows for non-contiguous data to be included, as we have to index the model predicted
    ## incidences by the days in the data (which might not be contiguous).
    ## This bit of code extracts these days and creates a vector with each day explicitly included, along
    ## with the corresponding location that it refers to. This means that we can pass these as vectors
    ## to the likelihood function rather than data frames which will aid computational speed.
    ## Must be done for incidence and microcephaly data seperately.

    ## Microcephaly data
    times_microceph<- NULL
    for(place in unique(data_locals)){
        tmpStartDays <- startDays[which(data_locals == place)]
        tmpEndDays <- endDays[which(data_locals == place)]
        days <- NULL
        for(i in 1:length(tmpStartDays)) days <- c(days, tmpStartDays[i]:tmpEndDays[i])
        days <- unique(days)
        times_microceph<- rbind(times_microceph, data.frame("local"=place,"valid_days"=days))
    }
    microceph_valid_days <- times_microceph$valid_days
    microceph_valid_local <- as.character(times_microceph$local)

    ## Incidence data
    inc_valid_days <- NULL
    inc_valid_local <- NULL
    if(!is.null(inc_startDays)){
        times_inc <- NULL
        for(place in unique(inc_locals)){
            tmpStartDays <- inc_startDays[which(inc_locals == place)]
            tmpEndDays <- inc_endDays[which(inc_locals == place)]
            days <- NULL
            for(i in 1:length(tmpStartDays)) days <- c(days, tmpStartDays[i]:tmpEndDays[i])
            days <- unique(days)
            times_inc<- rbind(times_inc, data.frame("local"=place,"valid_days"=days))

        }
        inc_valid_days <- times_inc$valid_days
        inc_valid_local <- as.character(times_inc$local)
    }
   
    ## If assuming normally distributed error
    if(version == "binomial"){
        ## May be using microcephaly data, or just ZIKV incidence data
        if(!is.null(microCeph)){
            f <- function(pars){
                lik <- posterior_complex_buckets(ts, pars, names, local,
                                                 startDays, endDays, buckets, microCeph, births, data_locals,
                                                 inc_startDays,inc_endDays,inc_locals,inc_buckets,inc_ZIKV,inc_NH,
                                                 peak_startDays, peak_endDays,peak_locals,
                                                 unique_locations,
                                                 microceph_valid_days, microceph_valid_local,
                                                 inc_valid_days, inc_valid_local, solver)
                ## If using priors, include here
                if(!is.null(PRIOR_FUNC)){
                    lik <- lik + PRIOR_FUNC(pars,...)
                }
                return(lik)
            }
        } else {
            ## If only using incidence data (ie. we might just want to fit R0)
            f <- function(values){
                lik <- posterior_inc(ts, pars, names, local,
                                     inc_startDays,inc_endDays,inc_locals,inc_buckets,inc_ZIKV,inc_NH,
                                     unique_locations,
                                     inc_valid_days, inc_valid_local, solver)
                ## If using priors, include here
                if(!is.null(PRIOR_FUNC)){
                    lik <- lik + PRIOR_FUNC(pars,...)
                }
                return(lik)
            }
        }
        ## This function is for using the version of the model that assumes changing behaviour etc, related to the forecasting analysis.
    } else if(version == "forecast"){
        f <- function(pars){
            lik <- posterior_known_inc_seir(pars, startDays, endDays,
                                            buckets, microCeph, births,
                                            inc_ZIKV, inc_NH, inc_buckets,
                                            inc_startDays, inc_endDays, solver)
            ## If using priors, include here
            if(!is.null(PRIOR_FUNC)){
                lik <- lik + PRIOR_FUNC(pars,...)
            }
            return(lik)        }
    } else {
        return(NULL)
    }
        
    return(f)
}

