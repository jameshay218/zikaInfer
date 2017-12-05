#' Simulation multi-state microcephaly and incidence data
#'
#' Given all of the parameters needed to solve the SEIR model, produces a data frame of simulated microcephaly affected births for multiple states. Produces data for each unique state in the parTab argument.
#' @param t_pars vector of two components - the total run time in days and the time step for the ODE solver
#' @param parTab Table of all parameters needed to solve the SEIR model and microcephaly simulation for each unique state. See \code{\link{setupParTable}}
#' @param weeks Defaults to FALSE. If true. returns weekly microcephaly births rather than by month.
#' @param dataRangeMicro a vector with lower and upper bounds on days for simulated microcephaly data
#' @param dataRangeInc a vector with lower and upper bounds on days for simulated ZIKV incidence data
#' @param noise bool indicating if noise should be added to the data
#' @param peakTimeRange the allowable prior range for ZIKV incidence peak time
#' @return a data frame of microcephaly births, total births and corresponding times; a data frame of ZIKV incidence data; a data frame describing the ZIKV epidemic peak times for each location included in the simulation
#' @param solver which ODE solver to use, either "rlsoda" or "lsoda"
#' @export
generate_multiple_data <- function(ts=seq(0,3003,by=1), parTab,weeks=FALSE,
                                   dataRangeMicro, dataRangeInc,noise=FALSE,
                                   peakTimeRange=NULL, solver="rlsoda", tstep=1){
    if(is.null(dataRangeMicro)) dataRangeMicro <- c(0,max(ts))
    if(is.null(dataRangeInc)) dataRangeInc <- c(0,max(ts))

    noYears <- floor(max(ts)/365)
    ## Create vector of how big the buckets should be - weeks or months
    if(!weeks) buckets <- rep(getDaysPerMonth(),noYears)
    else buckets <- rep(7, max(ts)/7)
    inc_buckets <- rep(7, max(ts)/7)
    
    ## Get all unique places that data should be created for
    places <- unique(parTab$local)
    places <- places[places!="all"]
    
    overallMicroDat <- NULL
    incDat <- NULL
    peakTimes <- NULL

    
    ## For each place
    for(place in places){
        ## Get state specific and universal parameters out of param table
        pars <- parTab[parTab$local %in% c("all",place),"values"]
        names(pars) <- parTab[parTab$local %in% c("all",place),"names"]

        ## Generate microcephaly curve
        probs <- generate_micro_curve(pars)

        ## Generate starting population sizes
        y0s <- generate_y0s(pars["N_H"], pars["density"])
        
        ## Solve the ODE model
        y <- solveSEIRModel(ts, y0s, pars, solver, compatible=TRUE)

        ######################################
        ## MICROCEPHALY DATA
        ######################################
        ## Generate simulated microcephaly data using solved ODE model and pregnancy risk curve
        y[y[,"I_M"] < 0,"I_M"] <- 0
        print(tstep)
        probM <- generate_probM(y[,"I_M"], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], tstep)*pars["propn"]
        probM <- average_buckets(probM, buckets/tstep)

        ## Total number of births each observation point is from populatoin size
        yearBirths <- pars["N_H"]/(pars["L_H"]/365)
        
        if(weeks) births <- yearBirths/52
        else births <- yearBirths/12
        
        births <- round(births)

        ## Certain proportion of all births are microcephaly, but only a given proportion of those are observed
        microDat <- probM*births

        if(noise){
            if(!any(parTab$names == "micro_sd")){
                microDat <- rbinom(length(probM), births, probM)
            } else {
                microDat <- rnorm(length(probM), births*probM, pars["micro_sd"])
                #microDat <- probM*births + err
            }
        }

        ## Needs to be to the nearest integer
        microDat <- sapply(microDat, function(x) max(0,x))
        microDat <- round(microDat)

######################################
        ## INCIDENCE DATA
######################################
        N_H <- as.integer(average_buckets(rowSums(y[,c("I_H","S_H","E_H","R_H")]), inc_buckets))
         ## Generate simulated incidence data
        inc <- diff(y[,"incidence"])
        inc[inc < 0] <- 0
        inc <- sum_buckets(inc, inc_buckets)
        perCapInc  <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]

        if(noise){
            if(!any(parTab$names=="inc_sd")){
                inc <- rbinom(length(perCapInc),pars["N_H"], perCapInc)
            } else {
                tmp <- rnorm(length(perCapInc),perCapInc*pars["N_H"],pars["inc_sd"])
                inc <- sapply(tmp, function(x) max(0,x))
            }
        } else {
            inc <- perCapInc*pars["N_H"]
        }
        
        ## Needs to be to the nearest integer
        inc <- round(inc)
######################################
        ## PEAK TIMES
######################################
        if(!is.null(peakTimeRange)){
            peakTime <- y[which.max(diff(y[,"incidence"])),"time"]
            lowerPeakTime <- peakTime - peakTimeRange/2
            upperPeakTime <- peakTime + peakTimeRange/2
            tmpPeakTime <- data.frame("actual"=peakTime,"start"=lowerPeakTime,"end"=upperPeakTime,"local"=as.character(place), stringsAsFactors=FALSE,row.names=NULL)
            peakTimes <- rbind(peakTimes, tmpPeakTime)
        }

######################################
        ## HOUSEKEEPING
######################################
        
        ## Add in times etc
        tmpMicroDat <- data.frame("startDay" = cumsum(buckets)-buckets,"endDay"=cumsum(buckets),"buckets"=buckets,"microCeph"=microDat,"births"=rep(births,length(buckets)),"local"=place,stringsAsFactors=FALSE)
        tmpIncDat <- data.frame("startDay" = cumsum(inc_buckets)-inc_buckets,"endDay" =cumsum(inc_buckets),"buckets"=inc_buckets,"inc"=inc,"N_H"=N_H,"local"=place,stringsAsFactors=FALSE)
        
        ## Only return data in the desired range
        tmpMicroDat <- tmpMicroDat[tmpMicroDat[,"startDay"] >= dataRangeMicro[1] & tmpMicroDat[,"endDay"] <= dataRangeMicro[2],]
        
        tmpIncDat <- tmpIncDat[tmpIncDat[,"startDay"] >= dataRangeInc[1] & tmpIncDat[,"endDay"] <= dataRangeInc[2],]
        
        incDat <- rbind(incDat, tmpIncDat)
        overallMicroDat <- rbind(overallMicroDat,tmpMicroDat)
    }
    return(list("micro"=overallMicroDat,"inc"=incDat,"peaks"=peakTimes))    
}
