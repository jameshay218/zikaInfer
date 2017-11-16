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
#' @param params Vector of parameters matching those returned by \code{\link{setupListPars}}
#' @param R0 desired R0 value
#' @return A single value for mosquito density
#' @export
#' @seealso \code{\link{r0.calc}}
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
#' @return a vector of peak times for each state in the parameter table
#' @export
find_peak_times <- function(parTab, ts=seq(0,3003,by=1)){
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
        
        y <- solveModelSimple_rlsoda(ts,y0s,pars,TRUE)
        ## Extract peak time
        peakTime <- y[which.max(diff(y[,"incidence"])),1]
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
#' @param ts the vector of times to solve the model over
#' @param runs number of draws to make to get uncertainty
#' @param origin the date from which to forecast from
#' @return a table of incidence and dates
#' @export
forecast_microceph <- function(chain,microParChain=NULL,parTab, ts,runs,origin="2013-01-01"){
    unique_states <- unique(parTab$local)
    unique_states <- unique_states[unique_states != "all"]
    
    allMicroBounds <- NULL
    for(local in unique_states){
        allMicro <- NULL

        ## Need random samples from recent chain and microcephaly saved chain
        samples <- sample(nrow(chain),runs)
        if(!is.null(microParChain)) sample_micro <- sample(nrow(microParChain),runs)
        index <- 1

        ## For each sample
        for(i in samples){
            ## Get these pars
            pars <- get_index_pars(chain, i)
            ## Also get the sampled microcephaly pars
            if(!is.null(microParChain)){
                micro_pars <- as.numeric(microParChain[sample_micro[index],])
                pars[colnames(microParChain)] <- micro_pars
            }
            tmp <- forecast_microceph_indiv(pars, local,parTab,ts)
            times <- tmp$time
            allMicro <- rbind(allMicro, tmp)
            index <- index + 1
        }
        colnames(allMicro) <- c("day","number")
        microBounds <- as.data.frame(reshape2::melt(sapply(unique(allMicro$day),function(x) c(quantile(allMicro[allMicro$day==x,"number"],c(0.025,0.5,0.975)),"mean"=mean(allMicro[allMicro$day==x,"number"])))))
        colnames(microBounds) <- c("quantile","time","micro")
        microBounds[,"time"] <- times[microBounds[,"time"]]
        microBounds$state <- local
        allMicroBounds <- rbind(microBounds, allMicroBounds)
    }
    means <- allMicroBounds[allMicroBounds$quantile=="mean",c("micro","time","state")]
    medians <- allMicroBounds[allMicroBounds$quantile=="50%",c("micro","time","state")]
    lower <- allMicroBounds[allMicroBounds$quantile=="2.5%",c("micro","time","state")]
    upper <- allMicroBounds[allMicroBounds$quantile=="97.5%",c("micro","time","state")]
    
    res <- plyr::join(means,medians,by=c("state","time"))
    res <- plyr::join(res, lower, by=c("state","time"))
    res <- plyr::join(res, upper, by=c("state","time"))

    colnames(res) <- c("mean","time","state","median","lower","upper")
    res$time <- as.Date(res$time,origin=origin)
    
    res <- res[,c("state","time","mean","median","lower","upper")]
    
    return(res)
    
}

#' Microcephaly forecast single pars
#'
#' Forecasts monthly microcephaly birth proportion for a single set of model parameters
#' @param pars the model parameters to solve
#' @param local the name of the Brazilian state
#' @param parTab the parameter table to draw names etc. from
#' @param ts the vector of times to solve the model over. Should be in days, but will be converted to months.
#' @param weeks if preferred, will return weekly proportions rather than months
#' @return a data frame of months and proportion microcephaly affected births
#' @export
forecast_microceph_indiv <- function(pars,local, parTab, ts, weeks=FALSE){
    noYears <- floor(max(ts)/365)

    ## Enumerate out either months or weeks
    if(weeks){
        weeks <- floor(max(ts)/7)
        ts <- seq(0,weeks*7,by=1)
        times <- rep(7, weeks)
    } else {
        times <- rep(getDaysPerMonth(),noYears)
        ts <- seq(0,noYears*365,by=1)
    }
    buckets <- times
    ## Get start and end days for each month - using mid point
    startDay <- cumsum(times) - times
    endDay <- cumsum(times)
    meanDay <- rowMeans(cbind(startDay,endDay))


    ## This is a bit too hard coded - FIX
    number <- which(unique(parTab$local)==local)-2
    
    ## Format parameter vector correctly
    state_pars <- parTab[parTab$local==local,"names"]
    if(number >= 1) state_pars <- paste(state_pars,".",number,sep="")
    state_pars <- pars[state_pars]
    names(state_pars) <- parTab[parTab$local==local,"names"]
    all_pars <- pars[parTab[parTab$local=="all","names"]]
    
    new_pars<- c(all_pars,state_pars)

     ## Generate actual incidence data for these parameters
    y0s <- generate_y0s(new_pars["N_H"],new_pars["density"])
    
    y <- solveSEIRModel_rlsoda(ts, y0s, new_pars,TRUE)
    inc <- diff(y[,"incidence"])
    inc[inc < 0] <- 0
    N_H <- rowSums(y[,5:8])
    inc_weeks <- floor(max(ts)/7)
    inc_buckets <- inc_times <- rep(7, inc_weeks)
    startDayInc <- cumsum(inc_times) - inc_times
    endDayInc <- cumsum(inc_times)
    meanDayInc <- rowMeans(cbind(startDayInc,endDayInc))

    inc <- sum_buckets(inc,inc_buckets)
    N_H <- average_buckets(N_H, inc_buckets)
    perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]

    ## Generate predicted microcephaly incidence for these parameters -
    ## need to restrict to predictions within the data range
    probs <- generate_micro_curve(new_pars)
    probM <- generate_probM(y[,"I_M"], new_pars["N_H"], probs, new_pars["b"], new_pars["p_MH"], new_pars["baselineProb"], 1)
    probM <- average_buckets(probM, buckets)
    ## Generate predicted microcephaly cases or per birth incidence depending on what was provided
    
    predicted <- probM
    return(list("microCeph"=data.frame("time"=meanDay,"microceph"=predicted),
                "ZIKV"=data.frame("time"=meanDayInc,"inc"=perCapInc)))
}

#' @export
create_forecast_normal <- function(local, parTab, ts, weeks=FALSE){
    f <- function(values){
        forecast_microceph_indiv(values,local, parTab, ts, weeks=FALSE)
    }
    return(f)
}

#' @export
forecast_known_inc_seir <- function(pars, startDays, endDays,
                           buckets, microCeph, births,
                           zikv, nh, inc_buckets,
                           inc_start, inc_end,
                           valid_days_micro, valid_days_inc,
                           perCap=FALSE){
    ## Generate microcephaly curve for these parameters
    probs <- generate_micro_curve(pars)
    ts <- seq(min(inc_start), max(inc_end)-1,by=1)
    switch_time_i <- pars["switch_time_i"]
    switch_time_m <- pars["switch_time_m"]
    switch_time_behaviour <- pars["switch_time_behaviour"]
    
    ## If this is one, then we're assuming that ZIKV reporting did not change
    check_par <- pars["zikv_reporting_change"]
    
    ## Get daily actual incidence based on data
    inc_1 <- zikv[which(inc_start < switch_time_i)]
    inc_1 <- inc_1/pars["incPropn"]
    if(check_par == 1){
        inc_2 <- zikv[which(inc_start >= switch_time_i)]/pars["incPropn2"]
    } else {
        inc_2 <- zikv[which(inc_start >= switch_time_i)]/pars["incPropn"]
    }
    inc <- c(inc_1,inc_2)
    
    inc <- rep(inc/nh/inc_buckets, inc_buckets)
    ## Generate probability of observing a microcephaly case on a given day
    probM_a <- generate_probM_forecast(inc, probs, pars["baselineProb"],
                                       pars["abortion_rate"], pars["birth_reduction"],
                                       which(ts == pars["switch_time_behaviour"]), pars["abortion_cutoff"], FALSE)

    ## Getting aborted births, ap_m(t)
    probM_b <- generate_probM_forecast(inc, probs, pars["baselineProb"],
                                       pars["abortion_rate"], pars["birth_reduction"],
                                       which(ts == pars["switch_time_behaviour"]), pars["abortion_cutoff"], TRUE)
    ## Get subset of times for the data we have incidence data
    probM <- probM_a/(1-probM_b)
    
    probM[which(ts < switch_time_m)] <- probM[which(ts < switch_time_m)]*pars["propn"]
    probM[which(ts >= switch_time_m)] <- probM[which(ts >= switch_time_m)]*pars["propn2"]
    probM <- average_buckets(probM,buckets)
    probM <- probM*births
    tmpStart <- startDays[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    tmpEnd <- endDays[which(startDays >= min(inc_start) & endDays <= max(inc_end))]
    ## Reported on mean of start and end report day
    meanDays <- (tmpStart + tmpEnd)/2
    return(data.frame("x"=meanDays,"y"=probM))
}

#' @export
create_f_forecast <- function(parTab,microDat,incDat){

    startDays <- microDat$startDay
    endDays <- microDat$endDays
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
        lik <- forecast_known_inc_seir(values, startDays, endDays,
                                       buckets, microCeph, births,
                                       zikv, nh, inc_buckets,
                                       inc_start, inc_end)
        return(lik)
    }
}
