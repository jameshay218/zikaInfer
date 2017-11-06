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
#' @useDynLib
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
#' @param pars the model parameters
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
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v2 <- function(pars){
    probs <- c(rep(pars["p1"], 14*7),rep(pars["p2"],14*7),rep(pars["p3"],12*7))
    return(unname(probs))
}

#' Microcephaly risk V3
#'
#' Microcephaly risk curve with 6 distinct periods (7, 7, 7, 7, 7 and 5 weeks)
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v3 <- function(pars){
    probs <- c(rep(pars["p1"], 7*7),rep(pars["p2"],7*7),rep(pars["p3"],7*7),rep(pars["p4"],7*7),rep(pars["p5"],7*7),rep(pars["p6"],5*7))
    return(unname(probs))
}

#' Microcephaly risk V4
#'
#' Microcephaly risk curve with 8 distinct periods
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v4 <- function(pars){
    probs <- c(rep(pars["p1"], 5*7),rep(pars["p2"],5*7),rep(pars["p3"],5*7),rep(pars["p4"],5*7),rep(pars["p5"],5*7),rep(pars["p6"],5*7),rep(pars["p7"],5*7),rep(pars["p8"],5*7))
    return(unname(probs))
}

#' Solve simple SEIR model 
#'
#' Given a list of parameters as generated by \code{\link{setupListPars}}, solves the simple SEIR model and returns a named data frame. Uses lsoda from ODE
#' @param ts vector of times to solve the ODE model over
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @param makenames indicates if solved data frame should be given column names
#' @return a data frame of the solved ODE model
#' @export
solveModelSimple_lsoda<- function(ts, y0s, pars,makenames=FALSE){
    ## Package ODE pars
    pars <- pars[c("L_M","L_H","D_EM","D_EH","D_IH","b","p_HM","p_MH","constSeed")]
    y <- deSolve::ode(y0s, ts, func="simpleSEIR",parms=pars,dllname="zikaProj",initfunc="initmodSEIR",nout=0, rtol=1e-5,atol=1e-5)
    if(makenames) colnames(y) <- c("time","S_M","E_M","I_M","S_H","E_H","I_H","R_H","incidence")
    return(y)
}

#' Solve simple SEIR model 
#'
#' Given a list of parameters as generated by \code{\link{setupListPars}}, solves the simple SEIR model and returns a named data frame. Uses rlsoda from rlsoda
#' @param ts vector of times to solve the ODE model over
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @param compatible indicates if solved data frame should be given in deSolve return format
#' @return a data frame of the solved ODE model
#' @export
solveModelSimple_rlsoda<- function(ts, y0s, pars,compatible=FALSE){
    pars <- pars[c("L_M","L_H","D_EM","D_EH","D_IH","b","p_HM","p_MH","constSeed")]
    rlsoda::rlsoda(y0s, ts, CsimpleSEIR_rich, pars, dllname="zikaProj", deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-6,rtol=1e-6)
}


#' Get days per month of the year
#'
#' Returns a vector of days per month
#' @param breaks total number of months to consider (defaults to 12 to give raw days per each month)
#' @return vector of days
#' @export
#' @useDynLib zikaProj
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
#' @useDynLib zikaProj
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
#' @useDynLib zikaProj
get_index_pars <- function(chain, index){
    tmpNames <- colnames(chain)[2:(ncol(chain)-1)]
    pars <- as.numeric(chain[index,2:(ncol(chain)-1)])
    names(pars) <- tmpNames
    return(pars)
}
    
#' Risk window names
#'
#' In case you forget which version corresponds to which parameters!
#' @export
#' @useDynLib zikaProj
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
#' @return a vector of peak times for each state in the parameter table
#' @export
find_peak_times <- function(parTab){
    ts <- seq(0,3003,by=1)
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
forecast_microceph_indiv <- function(pars,local, parTab, ts, weeks=FALSE){
    noYears <- floor(max(ts)/365)

    if(weeks){
        weeks <- floor(max(ts)/7)
        times <- seq(0,weeks*7,by=7)
        ts <- seq(0,weeks*7,by=1)
    } else {
        times <- rep(getDaysPerMonth(),noYears)
        ts <- seq(0,noYears*365,by=1)
    }

    startDay <- cumsum(times) - times
    endDay <- cumsum(times)
    meanDay <- rowMeans(cbind(startDay,endDay))
    buckets <- times
    
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
  
    y <- solveModelSimple_rlsoda(ts, y0s, new_pars,TRUE)
    ## Generate predicted microcephaly incidence for these parameters - need to restrict to predictions within the data range
    probs <- generate_micro_curve(new_pars)
    probM <- generate_probM(y[,"I_M"], new_pars["N_H"], probs, new_pars["b"], new_pars["p_MH"], new_pars["baselineProb"], 1)
    probM <- average_buckets(probM, buckets)
    ## Generate predicted microcephaly cases or per birth incidence depending on what was provided
    
    predicted <- probM
    return(data.frame("time"=meanDay,"microceph"=predicted))
}



#' @export
forecast_known_inc_seir <- function(pars, startDays, endDays,
                                    buckets, microCeph, births,
                                    zikv, nh, inc_buckets,
                                    inc_start, inc_end,
                                    valid_days_micro, valid_days_inc){
    ## Times after which reporting/birth behaviour changes
    switch_time_i <- pars["switch_time_i"]
    switch_time_m <- pars["switch_time_m"]
    switch_time_behaviour <- pars["switch_time_behaviour"]
    
    ## Generate microcephaly curve for these parameters
    probs <- generate_micro_curve(pars)

    ## Solve SEIR model
    y0s <- generate_y0s(as.numeric(pars["N_H"]),as.numeric(pars["density"]))
    y <- solveModelSimple_rlsoda(seq(0,3003,by=1), y0s, pars,FALSE)

    ## Calculate SEIR generated incidence for before switch time
    calc_inc <- diff(y["incidence",])
    calc_inc[calc_inc < 0] <- 0
    N_H <- colSums(y[5:8,])

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

    ## Convert to true underlying incidence
    inc_1 <- inc_1/pars["incPropn"]
    #inc_1 <- perCapInc*N_H/pars["incPropn"]

    ## Get incidence after switch time with new reporting rate
    inc_2 <- zikv[which(inc_start >= switch_time_i)]/pars["incPropn2"]
    inc <- c(inc_1,inc_2)
    
    ## Convert to daily incidence
    inc <- rep(inc/nh/inc_buckets, inc_buckets)

    ## Generate probability of observing a microcephaly caes on a given day
    #probM <- generate_probM_aux(inc, probs, pars["baselineProb"])
    ts <- seq(min(inc_start), max(inc_start),by=1)
    probM <- generate_probM_forecast(inc, probs, pars["baselineProb"],
                                     pars["abortion_rate"], pars["birth_reduction"],
                                     which(ts == pars["switch_time_behaviour"]), 12*7)
    
    
    probM[which(ts < switch_time_m)] <- probM[which(ts < switch_time_m)]*pars["propn"]
    probM[which(ts >= switch_time_m)] <- probM[which(ts >= switch_time_m)]*pars["propn2"]
    ## Births after switch time have reduction
                                        #births[which(startDays >= switch_time_m)] <- as.integer(births[which(startDays >= switch_time_m)]*(1-pars["birth_reduction"]))

    probM <- average_buckets(probM,buckets)
    
    return(probM*births)
}

#' @export
create_f_wow <- function(parTab,microDat,incDat){

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
