#' R0 calculation
#'
#' Calculates the R0 of the SEIR model given a vector of parameters. R0 defined as number of expected human cases given introduction of 1 infected human into a totally naive population of humans and mosquitoes.
#' @param params Vector of parameters matching those returned by \code{\link{setupListPars}}
#' @return A single value for R0
#' @export
#' @seealso \code{\link{b.calc}}
r0.calc <- function(params){
    NH <- params["N_H"]
    NM <- params["N_H"]*params["density"]
    muM <- 1/params["L_M"]
    sigmaM <- 1/params["D_EM"]

    muH <- 1/params["L_H"]
    gammaH <- 1/params["D_IH"]

    b <- params["b"]
    pHM <- params["p_HM"]
    pMH <- params["p_MH"]
    R0 <- (b^2*pHM*pMH*NM*sigmaM)/((sigmaM+muM)*muM*(gammaH+muH)*NH)
    return(unname(R0))
}


#' R0 vector calculation
#'
#' Calculates the R0 of the SEIR model given a matrix of parameters. R0 defined as number of expected human cases given introduction of 1 infected human into a totally naive population of humans and mosquitoes.
#' @param params Matrix of parameters matching those returned by \code{\link{setupListPars}}
#' @return A vector of values for R0
#' @export
#' @seealso \code{\link{b.calc}}
r0.vector <- function(params){
    NH <- params[,"N_H"]
    NM <- params[,"N_H"]*params[,"density"]
    muM <- 1/params[,"L_M"]
    sigmaM <- 1/params[,"D_EM"]

    muH <- 1/params[,"L_H"]
    gammaH <- 1/params[,"D_IH"]

    b <- params[,"b"]
    pHM <- params[,"p_HM"]
    pMH <- params[,"p_MH"]

    R0 <- (b^2*pHM*pMH*NM*sigmaM)/((sigmaM+muM)*muM*(gammaH+muH)*NH)
    return(R0)
}

#' Bite rate calculation
#'
#' Calculates the bite rate needed to generate a given R0 value, assuming that all other parameters are fixed.
#' @param params Vector of parameters matching those returned by \code{\link{setupListPars}}
#' @param R0 desired R0 value
#' @return A single value for bite rate
#' @export
#' @seealso \code{\link{r0.calc}}
b.calc <- function(params,R0){
    NH <- params["N_H"]
    NM <- params["N_H"]*params["density"]
    muM <- 1/params["L_M"]
    sigmaM <- 1/params["D_EM"]

    muH <- 1/params["L_H"]
    gammaH <- 1/params["D_IH"]

    b <- params["b"]
    pHM <- params["p_HM"]
    pMH <- params["p_MH"]
    
    b <- sqrt(((sigmaM+muM)*muM*(gammaH+muH)*NH)/((pHM * pMH * NM * sigmaM))*R0)
    return(unname(b))
    
}


#' Density calculation
#'
#' Calculates the density needed to generate a given R0 value, assuming that all other parameters are fixed.
#' @param params Vector of parameters matching those returned by setupListPars
#' @param R0 desired R0 value
#' @return A single value for mosquito density
#' @export density.calc
density.calc <- function(params,R0){
    NH <- params["N_H"]
    NM <- params["N_H"]*params["density"]
    muM <- 1/params["L_M"]
    sigmaM <- 1/params["D_EM"]

    muH <- 1/params["L_H"]
    gammaH <- 1/params["D_IH"]

    b <- params["b"]
    pHM <- params["p_HM"]
    pMH <- params["p_MH"]

    bot <- muM*(sigmaM+muM)*(gammaH + muH)*NH
    top <- (b^2)*pMH*pHM*sigmaM

    NM <- R0*bot/top
    density <- NM/NH

    return(unname(density))
}

#' Generates y0s
#'
#' Generates initial values for the simple SEIR model given population size and mosquito density
#' @param N_H human population size
#' @param density number of mosquitoes per person
#' @return a vector of initial population sizes
#' @export
generate_y0s <- function(N_H, density, iniI=10){
    N_M <- unname(N_H)*unname(density)
    S_M = 1*(N_M)
    E_M = 0
    I_M = 0

    S_H = unname(N_H) - iniI
    E_H = 0
    I_H = iniI
    R_H = 0

    incidence = 0
    
    return(c("S_M" = S_M, "E_M"=E_M,"I_M"=I_M, "S_H"=S_H, "E_H"=E_H,"I_H"=I_H,"R_H"=R_H, "incidence"=incidence))
}

#' Generate microceph curve
#'
#' Generates the microcephaly risk curve from model parameters
#' @param pars the model parameters
#' @return a vector of risks
#' @seealso \code{\link{microceph_v1}}, \code{\link{microceph_v2}}, \code{\link{microceph_v3}}, \code{\link{microceph_v4}}
#' @export
generate_micro_curve <- function(pars){
    if(!is.na(pars["mean"])) return(microceph_v1(pars))
    else if(!is.na(pars["p8"])) return(microceph_v4(pars))
    else if(!is.na(pars["p6"])) return(microceph_v3(pars))
    else return(microceph_v2(pars))
}

#' Microcephaly risk V1
#' 
#' Microcephaly risk curve under gamma distribution. 0 - 279 (all days)
#' @param pars the model parameters, specifying "mean" (gamma mean),"var" (gamma variance), and "c" (additional scaling constant)
#' @return the vector of risks
#' @export
microceph_v1 <- function(pars){
    mean <- pars["mean"]
    var <- pars["var"]

    scale <- var/mean
    shape <- mean/scale
 
    probs <- dgamma(0:279,shape=shape,scale=scale)*pars["c"]
    probs[probs > 1] <- 1

    return(probs)
}

#' Microcephaly risk V2
#'
#' Microcephaly risk curve with 3 distinct periods (14, 14 and 12 weeks)
#' @param pars the model parameters, "p1","p2",and "p3"
#' @return the vector of risks
#' @export
microceph_v2 <- function(pars){
    probs <- c(rep(pars["p1"], 14*7),rep(pars["p2"],14*7),rep(pars["p3"],12*7))
    return(unname(probs))
}

#' Microcephaly risk V3
#'
#' Microcephaly risk curve with 6 distinct periods (7, 7, 7, 7, 7 and 5 weeks)
#' @param pars the model parameters, "p1"-"p6"
#' @return the vector of risks
#' @export
microceph_v3 <- function(pars){
    probs <- c(rep(pars["p1"], 7*7),rep(pars["p2"],7*7),rep(pars["p3"],7*7),rep(pars["p4"],7*7),rep(pars["p5"],7*7),rep(pars["p6"],5*7))
    return(unname(probs))
}

#' Microcephaly risk V4
#'
#' Microcephaly risk curve with 8 distinct periods
#' @param pars the model parameters, "p1"-"p8"
#' @return the vector of risks
#' @export
microceph_v4 <- function(pars){
    probs <- c(rep(pars["p1"], 5*7),rep(pars["p2"],5*7),rep(pars["p3"],5*7),rep(pars["p4"],5*7),rep(pars["p5"],5*7),rep(pars["p6"],5*7),rep(pars["p7"],5*7),rep(pars["p8"],5*7))
    return(unname(probs))
}

#' Solve simple SEIR model
#'
#' Solves the SEIR model using the specified ode solver
#' @param ts vector of times to solve the ODE model over
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @param compatible if TRUE, makes the rlsoda output the same as lsoda output
#' @return a data frame of the solved ODE model
#' @export
solveSEIRModel <- function(ts, y0s, pars, solver="rlsoda", compatible=FALSE){
    if(solver=="rlsoda") return(solveSEIRModel_rlsoda(ts, y0s, pars, compatible))
    else return(solveSEIRModel_lsoda(ts, y0s, pars, TRUE))
}

#' Solve simple SEIR model (lsoda)
#'
#' Given a list of parameters as generated by \code{\link{setupListPars}}, solves the simple SEIR model and returns a named data frame. Uses lsoda from ODE
#' @param ts vector of times to solve the ODE model over
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @param makenames indicates if solved data frame should be given column names
#' @return a data frame of the solved ODE model
#' @export
solveSEIRModel_lsoda <- function(ts, y0s, pars,makenames=FALSE){
    ## Package ODE pars
    pars <- pars[c("L_M","L_H","D_EM","D_EH","D_IH","b","p_HM","p_MH","t0")]
    y <- deSolve::ode(y0s, ts, func="SEIR_model_lsoda",parms=pars,dllname="zikaProj",initfunc="initmodSEIR",nout=0, rtol=1e-5,atol=1e-5)
    if(makenames) colnames(y) <- c("time","S_M","E_M","I_M","S_H","E_H","I_H","R_H","incidence")
    return(y)
}

#' Solve simple SEIR model (rlsoda)
#'
#' Given a list of parameters as generated by \code{\link{setupListPars}}, solves the simple SEIR model and returns a named data frame. Uses rlsoda from rlsoda
#' @param ts vector of times to solve the ODE model over
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @param compatible indicates if solved data frame should be given in deSolve return format
#' @return a data frame of the solved ODE model
#' @export
solveSEIRModel_rlsoda <- function(ts, y0s, pars,compatible=FALSE){
    pars <- pars[c("L_M","L_H","D_EM","D_EH","D_IH","b","p_HM","p_MH","t0")]
    rlsoda::rlsoda(y0s, ts, C_SEIR_model_rlsoda, parms=pars, dllname="zikaProj", deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-5,rtol=1e-5)
}


#' Get days per month of the year
#'
#' Returns a vector of days per month
#' @param breaks total number of months to consider (defaults to 12 to give raw days per each month)
#' @return vector of days
#' @export
getDaysPerMonth <- function(breaks=12){
    days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
    days <- colSums(matrix(days,ncol=breaks))
    return(days)
}


#' Best pars
#'
#' Given an MCMC chain, returns the set of best fitting parameters
#' @param chain the MCMC chain
#' @return a name vector of the best parameters
#' @export
get_best_pars <- function(chain){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    bestPars <- as.numeric(chain[which.max(chain[,"lnlike"]),2:(ncol(chain)-1)])
    names(bestPars) <- tmpNames
    return(bestPars)
}

#' Index pars
#'
#' Given an MCMC chain, returns the parameters at the specified index
#' @param chain the MCMC chain
#' @param index the index
#' @return a named vector of the best parameters
#' @export
get_index_pars <- function(chain, index){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    pars <- as.numeric(chain[index,2:(ncol(chain)-1)])
    names(pars) <- tmpNames
    return(pars)
}
    
#' Risk window names
#'
#' In case you forget which version corresponds to which parameters!
#' @seealso \code{generate_micro_curve}, \code{\link{microceph_v1}}, \code{\link{microceph_v2}}, \code{\link{microceph_v3}}, \code{\link{microceph_v4}}
#' @export
print_version_names <- function(){
    print("Version 1 = Gamma curve")
    print("Version 2 = 3 pregnancy risk windows")
    print("Version 3 = 6 pregnancy risk windows")
    print("Version 4 = 8 pregnancy risk windows")
}

#' Sums nth elements
#'
#' Sums every n values of a numeric vector
#' @param v the vector to me summed
#' @param n the number of elements to sum
#' @export
sum_n_vector <- function(v,n){
    nv <- length(v)
    if (nv %% n)
      v[ceiling(nv / n) * n] <- NA
    colSums(matrix(v, n), na.rm=TRUE)
}

#' Find peak times
#'
#' Given a table of parameters as returned by \code{\link{setupParTable}}, finds the time of peak ZIKV incidence for each state
#' @param parTab table of parameters with names, values and corresponding state
#' @param ts vector of time in days in solve model over
#' @param solver which ODE solver to use
#' @return a vector of peak times for each state in the parameter table
#' @export
find_peak_times <- function(parTab, ts=seq(0,3003,by=1),solver="rlsoda"){
    unique_states <- unique(parTab$local)
    unique_states <- unique_states[unique_states != "all"]
    peakTimes <- NULL
    for(place in unique_states){
        ## Get the indices for this state in the data and parameter table
        indices_pars <- parTab$local == place | parTab$local == "all"

        ## Isolate these vectors
        pars <- parTab$values[indices_pars]
        names(pars) <- parTab$names[indices_pars]
        
        ## Generate starting y0 for the current parameter values
        y0s <- generate_y0s(as.numeric(pars["N_H"]),as.numeric(pars["density"]))
        
        y <- solveSEIRModel(ts,y0s,pars,solver)

        ## Extract peak time
        peakTime <- y[1,which.max(diff(y["incidence",]))]
        message(peakTime)
        peakTimes[place] <- peakTime
    }
    peakTimes
    
}

#' Microcephaly incidence forecast
#'
#' Given an MCMC chain, generates the incidence of microcephaly affected births
#' @param chain the MCMC chain to simulate from
#' @param microParChain the MCMC chain for microcephaly parameters
#' @param parTab the parameter table matching the generated MCMC chain
#' @param microDat the microcephaly data used for model fitting
#' @param incDat the ZIKV incidence data used for model fitting
#' @param ts the vector of times to solve the model over
#' @param runs number of draws to make to get uncertainty
#' @param origin the date from which to forecast from
#' @return a table of incidence and dates
#' @export
forecast_microceph <- function(chain,microParChain=NULL,parTab,
                               microDat, incDat, ts,runs,
                               origin="2013-01-01"){
    unique_states <- unique(parTab$local)
    unique_states <- unique_states[unique_states != "all"]
    
    allMicroBounds <- NULL
    allIncBounds <- NULL
    for(local in unique_states){
        allMicro <- NULL
        allInc <- NULL
        ## Need random samples from recent chain and microcephaly saved chain
        samples <- sample(nrow(chain),runs)
        if(!is.null(microParChain)) sample_micro <- sample(nrow(microParChain),runs)
        index <- 1

        ## For each sample
        for(i in samples){
            tmpMicro <- microDat[microDat$local == local,]
            tmpInc <- incDat[incDat$local == local,]
            
            ## Get these pars
            pars <- get_index_pars(chain, i)
            ## Also get the sampled microcephaly pars
            if(!is.null(microParChain)){
                micro_pars <- as.numeric(microParChain[sample_micro[index],])
                pars[colnames(microParChain)] <- micro_pars
            }
            ## This is a bit too hard coded - FIX
            number <- which(unique(parTab$local)==local)-2
            
            ## Format parameter vector correctly
            state_pars <- parTab[parTab$local==local,"names"]
            if(number >= 1) state_pars <- paste(state_pars,".",number,sep="")
            state_pars <- pars[state_pars]
            names(state_pars) <- parTab[parTab$local==local,"names"]
            all_pars <- pars[parTab[parTab$local=="all","names"]]
            
            new_pars<- c(all_pars,state_pars)
            
            tmp <- forecast_model_normal(new_pars, tmpMicro, tmpInc, ts, TRUE)
            microPred <- tmp$microCeph
            incPred <- tmp$ZIKV
            allMicro <- rbind(allMicro,microPred)
            allInc <- rbind(allInc, incPred)
            index <- index + 1
        }
        ## Save and get bounds for microcephaly data
        colnames(allMicro) <- c("microCeph","time")
        microBounds <- as.data.frame(reshape2::melt(sapply(unique(allMicro$time),function(x) c(quantile(allMicro[allMicro$time==x,"microCeph"],c(0.025,0.5,0.975)),"mean"=mean(allMicro[allMicro$time==x,"microCeph"])))))
        colnames(microBounds) <- c("quantile","time","microCeph")
        microBounds[,"time"] <- times[microBounds[,"time"]]
        microBounds$state <- local
        allMicroBounds <- rbind(microBounds, allMicroBounds)

        ## Save and get bounds for incidence data
        colnames(allInc) <- c("time","inc")
        incBounds <- as.data.frame(reshape2::melt(sapply(unique(allInc$time),function(x) c(quantile(allInc[allInc$time==x,"inc"],c(0.025,0.5,0.975)),"mean"=mean(allInc[allInc$time==x,"inc"])))))
        colnames(incBounds) <- c("quantile","time","inc")
        incBounds[,"time"] <- times[incBounds[,"time"]]
        incBounds$state <- local
        allIncBounds <- rbind(incBounds, allIncBounds)
      
    }
    ## Consolidate results for microcephaly incidence
    means <- allMicroBounds[allMicroBounds$quantile=="mean",c("microCeph","time","state")]
    medians <- allMicroBounds[allMicroBounds$quantile=="50%",c("microCeph","time","state")]
    lower <- allMicroBounds[allMicroBounds$quantile=="2.5%",c("microCeph","time","state")]
    upper <- allMicroBounds[allMicroBounds$quantile=="97.5%",c("microCeph","time","state")]
    
    res <- plyr::join(means,medians,by=c("state","time"))
    res <- plyr::join(res, lower, by=c("state","time"))
    res <- plyr::join(res, upper, by=c("state","time"))

    colnames(res) <- c("mean","time","state","median","lower","upper")
    res$time <- as.Date(res$time,origin=origin)
    
    microCephRes <- res[,c("state","time","mean","median","lower","upper")]


    ## Consolidate results for ZIKV incidence
    means <- allIncBounds[allIncBounds$quantile=="mean",c("inc","time","state")]
    medians <- allIncBounds[allIncBounds$quantile=="50%",c("inc","time","state")]
    lower <- allIncBounds[allIncBounds$quantile=="2.5%",c("inc","time","state")]
    upper <- allIncBounds[allIncBounds$quantile=="97.5%",c("inc","time","state")]
    
    res <- plyr::join(means,medians,by=c("state","time"))
    res <- plyr::join(res, lower, by=c("state","time"))
    res <- plyr::join(res, upper, by=c("state","time"))

    colnames(res) <- c("mean","time","state","median","lower","upper")
    res$time <- as.Date(res$time,origin=origin)
    
    incRes <- res[,c("state","time","mean","median","lower","upper")]
    
    return(list("microCeph"=microCephRes,"ZIKV"=incRes))
}


#' @export
forecast_known_inc_seir <- function(pars, startDays, endDays,
                                    buckets, births,microCeph,
                                    zikv, nh, inc_buckets,
                                    inc_start, inc_end,
                                    perCap=FALSE, solver="rlsoda"){
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
    ## THIS IS WHAT WE PLOT AS THE MODEL FIT
    perCapInc <- (1-(1-(calc_inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
    if(!perCap) perCapInc <- perCapInc*N_H
    tmpIncStart <- inc_start[which(inc_start < switch_time_i)]
    tmpIncEnd <- inc_end[which(inc_start < switch_time_i)]
    tmpIncMean <- (tmpIncStart + tmpIncEnd)/2
    modelDat <- data.frame("time"=tmpIncMean,"inc"=perCapInc)

    
    ## Get subset of incidence for these times
    inc_1 <- zikv[which(inc_start < switch_time_i)]

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
    ## THIS IS WHAT WE PLOT AS FORECASTED INCIDENCE
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

    ## THIS IS THE FORECASTED MICROCEPHALY CASES
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
    aborted_births_predicted <- aborted_births[which(startDays >= prediction_time)]
    
    ## From this, we can infer the true number of births that happened assuming that some births
    ## were aborted and a proportion of the forecasted births did no occur
    inferred_births <- (1-pars["avoided_births"])*predicted_births - aborted_births_predicted
    
    ## The actual number of births after the prediction time are inferred
    births2 <- tmp_births
    births2[which(startDays >= prediction_time)] <- as.integer(inferred_births)

    if(!perCap) probM <- births2*probM
    
    tmpStart <- startDays[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    tmpEnd <- endDays[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    ## Reported on mean of start and end report day
    meanDays <- (tmpStart + tmpEnd)/2
    
    microCephData <- data.frame("time"=meanDays,"microCeph"=probM)

    aborted_data <- data.frame("time"=meanDays,"aborted"=aborted_births)
    
    return(list("microCeph"=microCephData,"ZIKV"=modelDat,"aborted"=aborted_data))
}


#' @export
forecast_model_normal <- function(pars, microDat,incDat, ts, perCap=FALSE, solver="rlsoda"){
    ## Generate actual incidence data for these parameters
    y0s <- generate_y0s(pars["N_H"],pars["density"])
    y <- solveSEIRModel(ts, y0s, pars,solver)
    
    ## Generate predicted microcephaly incidence for these parameters - need to restrict to predictions within the data range
    probs <- generate_micro_curve(pars)
    probM <- generate_probM(y["I_M",], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)
    probM <- probM[which(y["time",] >= min(microDat[,"startDay"]) & y["time",] <= max(microDat[,"endDay"]))]
    probM <- average_buckets(probM, microDat[,"buckets"])
    
    ## Generate predicted microcephaly cases or per birth incidence depending on what was provided
    if(perCap){
        predicted <- probM*pars["propn"]
    } else {
        predicted <- probM*microDat[,"births"]*pars["propn"]
    }
    meanDay <- (microDat$startDay + microDat$endDay)/2
    predicted <- data.frame(time=meanDay,microCeph=predicted)
    
    ## Generate predicted incidence cases
    N_H <- average_buckets(colSums(y[5:8,]), incDat$buckets)
    tmpY <- y[,which(y["time",] >= min(incDat[,"startDay"]) & y["time",] <= max(incDat[,"endDay"]))]
    inc <- diff(tmpY["incidence",])
    ## At the moment this really needs to be in weeks
    inc <- sum_buckets(inc, incDat$buckets)
    perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
    if(!perCap) perCapInc <- perCapInc*N_H
    
    meanDayInc <- (incDat$startDay + incDat$endDay)/2
    y <- data.frame(time=meanDayInc,inc=perCapInc)
    return(list("microCeph"=predicted,"ZIKV"=y))
}



#' @export
create_forecast_function <- function(parTab,microDat,incDat, ts=seq(0,3003,by=1), singleWave=TRUE){

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
    if(singleWave){
        f <- function(values, perCap=FALSE){
            names(values) <- names_pars
            y <- forecast_model_normal(values, microDat, incDat, ts, perCap)
            return(y)
        }

    } else {        
        f <- function(values, perCap=FALSE){
            names(values) <- names_pars
            y <- forecast_known_inc_seir(values, startDays, endDays,
                                         buckets, births,microCeph,
                                         zikv, nh, inc_buckets,
                                         inc_start, inc_end, perCap)
            return(y)
        }
    }
    return(f)
}

