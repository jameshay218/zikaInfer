#' Parameter setup
#'
#' Sets up the parameter table based on the version, locations and allowable parameters
#' @param locationNames the vector of locations to include
#' @param version Version 1-3 indicating which form of the microcephaly risk curve is to be used
#' @param locationInfo a table with location name, L_H (human lifespan in years) and N_H (human population size). See data(locationInfo)
#' @param useInc bool indicating whether or not ZIKV incidence data is considered here
#' @param allowablePars data frame of allowable parameters for density and t0
#' @param sharedProb bool if the version with shared microcephaly parameters should be used
#' @param normLik Using normal likelihood function for data. If false, uses binomial likelihood. If true, uses normal likelihood which incorporates two standard deviation parameters for ZIKV and microcephaly incidence data
#' @param locationWeights bool deciding if data from different locations should be weighted
#' @param forecast if TRUE, sets up the parameter table for the forecasting analysis
#' @return the modified parameter table
#' @export
parTabSetup <- function(locationNames, version,locationInfo,
                         useInc,allowablePars=NULL,sharedProb=TRUE,
                         normLik=FALSE,locationWeights=FALSE, forecast=FALSE){
    tmpLocations <- locationInfo[locationInfo$rawName %in% locationNames,]
    norm_location <- tmpLocations[which.min(tmpLocations$incidenceOrder),"rawName"]
    ## Get first location from the correct order found in passed locations
    ## Fix reporting proportion of microcephaly relative to this
    message(cat("Reference location: ",norm_location,"\n",sep=""))
    
    parTab <- parTabSetup_aux(locationInfo,version, sharedProb=sharedProb,
                              locationNormLik=normLik,locationWeights=locationWeights,
                              forecast=forecast)
    if(locationWeights) parTab[parTab$names=="location_weight","values"] <- 1/length(locationNames)
    else parTab[parTab$names=="location_weight","values"] <- 1
    ## If not using incidence data, need to fix incidence proportion and baseline inc parameters
    if(!useInc){
        incDat <- NULL
        parTab[parTab$names=="incPropn","fixed"] <- 1
        parTab[parTab$names=="incPropn","values"] <- 1
        parTab[parTab$names=="baselineInc","values"] <- 0
        parTab[parTab$names=="baselineInc","fixed"] <- 1
        ## If not using incidence data, give inc data zero weighting
        parTab[parTab$names=="inc_weight","values"] <- 0
    } else {
        parTab[parTab$names=="incPropn","fixed"] <- 0
        ## Incidence reporting rate likely very small
        parTab[parTab$names == "incPropn","values"] <- runif(length(parTab[parTab$names == "incPropn","values"]),0.0001,0.001)
        ## Give incidence data 50% weight
        parTab[parTab$names=="inc_weight","values"] <- 0.5
        if(normLik){
            parTab[parTab$names=="inc_sd","fixed"] <- 0
        }
    }
    ## Arbritrarily fix one location's reporting proportion
    parTab[parTab$names=="propn" & parTab$local== norm_location,c("fixed","values")] <- c(1,1)
    ## If there's only one location, fix propn, as this is the same as 'c'
    if(length(locationNames) == 1) parTab[parTab$names=="propn" & parTab$local== locationNames[1],c("fixed","values")] <- c(1,1)

    ## If using normal likelihood function, need to allow SD to vary.
    if(normLik) parTab[parTab$names=="micro_sd","fixed"] <- 0
        
      ## If there's a data frame of allowable starting parameters, use these to seed density and t0
    if(!is.null(allowablePars)){
        for(location in locationNames){
            tmpTab <- parTab[parTab$local %in% c("all",location),c("values","names")]
            pars <- tmpTab$values
            names(pars) <- tmpTab$names

            ## Limit density to give an R0 between 1 and 10 - anything above around 6 just gives the same likelihood
            min_density <- density.calc(pars,1)
            max_density <- density.calc(pars,10)
            parTab[parTab$local == location & parTab$names == "density",c("lower_bounds","start_lower")] <- min_density
            parTab[parTab$local == location & parTab$names == "density",c("upper_bounds","start_upper")] <- max_density
            
            tmpAllowable <- allowablePars[allowablePars$local == location,]
            index <- as.integer(runif(1,1,nrow(tmpAllowable)))
            parTab[parTab$names == "t0" & parTab$local == location,"values"] <- as.numeric(tmpAllowable[index,"t0"])
            parTab[parTab$names == "density" & parTab$local == location,"values"] <- as.numeric(tmpAllowable[index,"density"])
        }
    }
    parTab <- parTab[parTab$local %in% c("all",locationNames),]
    return(parTab)
}

#' MCMC start parTab
#'
#' Generates random starting points from lower_bounds and upper_bounds to seed the MCMC chain for free parameters
#' @param parTab the parameter table to generate starting positions from
#' @param peakTimes the ODE solver and likelihood function is a bit tempermental to the choice of starting parameters for R0 and the epidemic seed time, t0.
#' @param restrictedR0 if TRUE, sets starting points for density and t0 based on generated allowable values
#' @param allowableParsFile set to "" if not available. These allowable start points take a while to generate (\code{\link{generateAllowableParams}}), so it can sometimes be helpful to pre-compute these and then pass the file location
#' @param allowableStarts the pre-computed matrix of allowable parameters
#' @return a parameter table in the same form as that which was passed in, but with randomised values
#' @export
generateStartingParTab <- function(parTab, peakTimes=NULL,restrictedR0=TRUE,allowableParsFile="", allowableStarts=NULL){
    startTab <- parTab
    stateNames <- unique(startTab$local)
    stateNames <- stateNames[stateNames != "all"]
    
    ## Generate random starting values
    for(i in which(startTab$fixed==0)) startTab[i,"values"] <- runif(1,startTab$start_lower[i],startTab$start_upper[i])

    ## If restricting state points for R0, set these seperately
    if(restrictedR0){
        if(is.null(allowableStarts)){
            allowableStarts <- generateAllowableParams(peakTime=NULL, peakTimeRange=120, stateNames=stateNames,
                                                       parTab=parTab,allowableParsFile=allowableParsFile,R0max=7,peakTimings=peakTimes)
        }
        for(local in stateNames){
            tmp_starts <- allowableStarts[allowableStarts$local == local,]
            startTab[startTab$local == local & startTab$names == "density","values"] <- as.numeric(tmp_starts[runif(1,1,nrow(tmp_starts)),"density"])
            startTab[startTab$local == local & startTab$names == "t0","values"] <- as.numeric(tmp_starts[runif(1,1,nrow(tmp_starts)),"t0"])
        }
    }
    return(startTab)
}


#' Defaults ODE parameters
#'
#' Generates the default vector of parameters for the ODE model. Can specify which version of the model to use. Passing 1 gives the simple model with no age classes, otherwise gives the parameters for the model with age classes.
#' @param version the version of the model to be used (see \code{\link{print_version_names}})
#' @return vector of named parameter values
#' @export
setupParsODE <- function(version=1){
    D_EM=10.5
    D_EH=5.9
    D_IH=5
    L_M=14
    L=73.6*365
    D_C=18*365
    D_F=3/12*365
    b=100/365
    p_MH=0.5
    p_HM=0.5

    if(version <= 4) pars <- c("L_M"=L_M,"D_EM"=D_EM,"D_EH"=D_EH,"D_IH"=D_IH,"b"=b,"p_HM"=p_HM,"p_MH"=p_MH)
    else pars <- c("L_M"=L_M,"L_H"=L,"D_EM"=D_EM,"D_EH"=D_EH,"D_IH"=D_IH,"b"=b,"p_HM"=p_HM,"p_MH"=p_MH)
    return(pars)
}

#' Default simulation parameters
#'
#' Generates the default vector of parameters for the simulation. Can specify which version of the model to use. Passing 1 gives the simple model with no age classes, 2 gives the parameters for the model with age classes. 3 gives the parameters for the age class model with parameters for head circumference distributions.
#' @param version the version of the model being run
#' @return vector of names parameter values
#' @export
#' @seealso \code{\link{setupParsODE}}
setupParsLong <- function(version = 1){
    pars <- setupParsODE(version)
    baseline <- 0.0002
    mean <- 17
    var <- 3
    c <- 0.15
    p1 <- 0.1
    p2 <- 0.1
    p3 <- 0.1
    p4 <- 0.1
    p5 <- 0.1
    p6 <- 0.1
    p7 <- 0.1
    p8 <- 0.1
    inc_weight <- 0
    
    if(version ==1){
        pars <- c("location_weight"=1,"inc_weight"=inc_weight,"baselineProb"=baseline,pars,"mean"=mean,"var"=var,"c"=c)
    } else if(version ==2){
        pars <- c("location_weight"=1,"inc_weight"=inc_weight,"baselineProb"=baseline,pars, "p1"=p1, "p2"=p2, "p3"=p3)
    } else if(version==3){
        pars <- c("location_weight"=1,"inc_weight"=inc_weight,"baselineProb"=baseline,pars, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6)
    } else if(version==4){
        pars <- c("location_weight"=1,"inc_weight"=inc_weight,"baselineProb"=baseline,pars, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6, "p7"=p7, "p8"=p8)  
    } else {
        pars <- c("location_weight"=1,"inc_weight"=inc_weight,"baselineProb"=baseline,pars,"mean"=mean,"var"=var,"c"=c)
    }
    return(pars)   
}

#' Default parameter list setup
#'
#' Sets up the default list of parameters and initial conditions to be used by \code{\link{solveSEIRModel}}. Can specify which version of the model to use. Passing 1 gives the simple model with no age classes, otherwise gives the parameters for the model with age classes.
#' @param duration Duration in years of the simulation. Defaults to 10 years
#' @param version the version of the model that is going to be used (see \code{\link{print_version_names}})
#' @return list of times, initial conditions, and ODE parameters
#' @export
#' @seealso \code{\link{setupParsLong}}
setupListPars <- function(duration=2*365, version=1){
    pars <- setupParsLong(version)
    y0 <- generate_y0s(pars["N_H"], pars["density"])
    ts <- seq(0,duration,by=1)
    return(list(ts,y0,pars))
}


#' MCMC parameter setup
#'
#' Sets up the parameter table for use in MCMC
#' @param locationInfo data frame of location info, see data(locationInfo). Needs columns for L_H, N_H and rawName
#' @param version which model version is being run (see \code{\link{print_version_names}})
#' @param sharedProb boolean indicating if all locations share the same risk curve
#' @param parFile if already specified, where to find the parameter table file without location parameters
#' @param locationParFile as parFile, but with location specific parameters added
#' @param locationNormLik boolean deciding whether or not a normal likelihood function is used (otherwise binomial). If NULL, uses the binomial. If not NULL, true or false indicates whether the likelihood standard deviations are location specific or not.
#' @param locationWeights bool deciding if location weightings should be used (otherwise all locations have equal weighting)
#' @param forecast if TRUE, sets up the parameter table for the forecasting analysis
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, boolean for log scale, initial step sizes, log proposal and whether or not the parameter should be fixed.
#' @export
parTabSetup_aux <- function(locationInfo, version=1,sharedProb=FALSE,
                            parFile = "", locationParFile = "",
                            locationNormLik=FALSE,locationWeights=FALSE,
                            forecast=FALSE){
    ## Checks if given filename exists. If not, creates a fresh parameter table
    if(file.exists(parFile)){
        print(paste("Reading in: ",parFile,sep=""))
        paramTable <- read.table(parFile, header=TRUE,sep=",",stringsAsFactors=FALSE)
    }
    else paramTable <- createParTable(NULL)
    
    useNames <- names(setupParsLong(version))
    if(locationWeights) useNames <- useNames[useNames != "location_weight"]
    ## Gets the parameter names used for each version of the model
    allMicroPars <- c("mean","var","c","p1","p2","p3","p4","p5","p6","p7","p8")
    if(version==2){
        microPars <- c("p1","p2","p3")
    } else if(version==3) {
        microPars <- c("p1","p2","p3","p4","p5","p6")
    } else if(version==4) {
        microPars <- c("p1","p2","p3","p4","p5","p6","p7","p8")
    } else {
        microPars <- c("mean","var","c")
    }
    notMicroPars <- allMicroPars[!(allMicroPars %in% microPars)]
    if(locationNormLik){
        useNames <- c(useNames, "micro_sd","inc_sd")
    }

    if(forecast){
        useNames <- c(useNames, "switch_time_i","switch_time_m","switch_time_behaviour","predicted_births","birth_reduction", "abortion_rate","avoided_births","abortion_cutoff", "zikv_reporting_change")
    }
    
    ## Remove all unneccessary parameters
    paramTable <- paramTable[paramTable[,"names"] %in% useNames,]
    if(!sharedProb) paramTable <- paramTable[!(paramTable[,"names"] %in% microPars),]

    ## Adds location specific parameters
    locationParTable <- NULL
    if(!is.null(locationInfo) | file.exists(locationParFile)){
        locationParTable <- setupLocationParTable(locationInfo, locationParFile)
        locationParTable <- locationParTable[!(locationParTable[,"names"] %in% paramTable$names),]
        locationParTable <- locationParTable[!(locationParTable[,"names"] %in% notMicroPars),]
        if(!locationNormLik) locationParTable <- locationParTable[!(locationParTable$names %in% c("micro_sd","inc_sd")),]
        if(!forecast)  locationParTable <- locationParTable[!(locationParTable$names %in% c("incPropn2","propn2")),]
        if(sharedProb) locationParTable <- locationParTable[!(locationParTable[,"names"] %in% microPars),]
    }

    paramTable <- rbind(paramTable, locationParTable)
    
    paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")] <- lapply(paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")], FUN=as.numeric)
    
    return(paramTable)
}

#' MCMC parameter table creation
#'
#' Sets up the parameter table for use in MCMC
#' @param saveFile if provided, writes the parameter table to the given filename
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, initial step sizes and whether or not the parameter should be fixed.
#' @export
createParTable <- function(saveFile=NULL){
    names <- c("inc_weight","baselineProb","L_M","D_EM","L_H","D_C","D_F","D_EH","D_IH",
               "b","p_HM","p_MH","t0","mean","var","c","tstep","p1","p2","p3","p4","p5","p6","p7","p8",
               "inc_sd","micro_sd","location_weight",
               "switch_time_i","switch_time_m","switch_time_behaviour","predicted_births","birth_reduction",
               "abortion_rate","avoided_births","abortion_cutoff", "zikv_reporting_change")    
    paramTable <- matrix(0, ncol=9, nrow=length(names))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")
    paramTable[,"names"] <- names
    paramTable$names <- as.character(paramTable$names)
    
    paramTable[paramTable[,"names"]=="inc_weight",2:ncol(paramTable)] <- c(0,"all",0,1,0,1,0,1)
    paramTable[paramTable[,"names"]=="baselineProb",2:ncol(paramTable)] <- c(0.0002,"all",0,1,0.1,0,0,0.001)
    paramTable[paramTable[,"names"]=="L_M",2:ncol(paramTable)] <- c(5,"all",0,100,0.1,1,2,14)
    paramTable[paramTable[,"names"]=="D_EM",2:ncol(paramTable)] <- c(8.4,"all",0,100,0.1,1,7.3,9.3)
    paramTable[paramTable[,"names"]=="L_H",2:ncol(paramTable)] <- c(365*70,"all",0,200*365,0.1,1,0,200*365)
    paramTable[paramTable[,"names"]=="D_C",2:ncol(paramTable)] <- c(365*18,"all",0,25*365,0.1,1,0,25*365)
    paramTable[paramTable[,"names"]=="D_F",2:ncol(paramTable)] <- c(0.75*365,"all",0,365,0.1,1,0,365)
    paramTable[paramTable[,"names"]=="D_EH",2:ncol(paramTable)] <- c(4,"all",0,100,0.1,1,2.3,6)
    paramTable[paramTable[,"names"]=="D_IH",2:ncol(paramTable)] <- c(6,"all",0,100,0.1,1,4,8)
    paramTable[paramTable[,"names"]=="b",2:ncol(paramTable)] <- c(0.5,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="p_HM",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="p_MH",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="t0",2:ncol(paramTable)] <- c(0,"all",0,1000,0.1,1,400,600)
    paramTable[paramTable[,"names"]=="mean",2:ncol(paramTable)] <- c(100,"all",0,10000,0.1,0,60,200)
    paramTable[paramTable[,"names"]=="var",2:ncol(paramTable)] <- c(500,"all",0,100000,0.1,0,100,1000)
    paramTable[paramTable[,"names"]=="c",2:ncol(paramTable)] <- c(3,"all",0,100,0.1,0,0.1,5)
    paramTable[paramTable[,"names"]=="tstep",2:ncol(paramTable)] <- c(7,"all",0,100,0.1,1,0,7)
    paramTable[paramTable[,"names"]=="p1",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p2",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p3",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p4",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p5",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p6",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p7",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p8",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="inc_sd",2:ncol(paramTable)] <- c(1,"all",0,100,0.1,1,0.5,10)
    paramTable[paramTable[,"names"]=="micro_sd",2:ncol(paramTable)] <- c(1,"all",0,100,0.1,1,0.5,10)
    paramTable[paramTable[,"names"]=="location_weight",2:ncol(paramTable)] <- c(1,"all",0,1,0.1,1,0.1,1)

    ## Related to forecasting method
    paramTable[paramTable[,"names"]=="switch_time_i",2:ncol(paramTable)] <- c(1044,"all",860,1064,0.1,1,860,1064)
    paramTable[paramTable[,"names"]=="switch_time_m",2:ncol(paramTable)] <- c(1167,"all",1000,1300,0.1,1,1000,1300)
    paramTable[paramTable[,"names"]=="switch_time_behaviour",2:ncol(paramTable)] <- c(1044,"all",1000,2000,0.1,1,1000,1300)
    paramTable[paramTable[,"names"]=="predicted_births",2:ncol(paramTable)] <- c(1037,"all",1000,2000,0.1,1,1000,2000)
    paramTable[paramTable[,"names"]=="birth_reduction",2:ncol(paramTable)] <- c(0,"all",0,1,0.1,0,0,1)
    paramTable[paramTable[,"names"]=="abortion_rate",2:ncol(paramTable)] <- c(0,"all",0,1,0.1,0,0,1)
    paramTable[paramTable[,"names"]=="avoided_births",2:ncol(paramTable)] <- c(0,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="abortion_cutoff",2:ncol(paramTable)] <- c(168,"all",0,300,0.1,1,0,300)
    paramTable[paramTable[,"names"]=="zikv_reporting_change",2:ncol(paramTable)] <- c(0,"all",0,1,0.1,1,0,1)
    
    if(!is.null(saveFile)) write.table(paramTable,saveFile,row.names=FALSE,sep=",")
    
    return(paramTable)
}

#' MCMC location-specific parameter table creation
#'
#' Sets up the parameter table for all locations for use in \code{\link{run_MCMC}}
#' @param locationInfo data frame of locations - with L_H, N_H and rawName
#' @param saveFile if provided, writes the parameter table to the given filename
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, initial step sizes and whether or not the parameter should be fixed.
#' @export
createLocationParTable <- function(locationInfo, saveFile = NULL){
    places <- as.character(unique(locationInfo$rawName))
    numberPars <- 23
    paramTable <- matrix(0, ncol=9, nrow=numberPars*length(places))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local", "lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")
    index <- 1
    
    for(place in places){
        tmpDat <- locationInfo[locationInfo$rawName==place,]
        paramTable[index,] <- c("L_H",tmpDat[1,"L_H"]*365,place,0,200*365,0.1,1,0,200*365)
        index <- index + 1
        paramTable[index,] <- c("N_H", tmpDat[1,"N_H"], place, 0,100000000,0.1,1,0,1000000000)
        index <- index + 1
        paramTable[index,] <- c("density",3,place,0,100,0.1,0,3,6)
        index <- index + 1
        paramTable[index,] <- c("propn",1,place,0,1,0.1,0,0.01,0.2)
        index <- index + 1
        paramTable[index,] <- c("t0",400,place,0,2000,0.1,0,400,600)
        index <- index + 1
        paramTable[index,] <- c("mean",100,place,0,2000,0.1,0,20,200)
        index <- index + 1
        paramTable[index,] <- c("var", 500, place, 0,100000,0.1,0,50,10000)
        index <- index + 1
        paramTable[index,] <- c("c",3,place,0,100,0.1,0,0.1,5)
        index <- index + 1
        paramTable[index,] <- c("incPropn",0.0025,place,0,1,0.1,0,0.001,0.005)
        index <- index + 1
        paramTable[index,] <- c("baselineInc",0.0001,place,0,1,0.1,0,0.00005,0.001)
        index <- index + 1
        paramTable[index,] <- c("p1",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p2",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p3",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p4",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p5",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p6",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p7",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p8",0.1,place,0,1,0.1,0,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("inc_sd",1,place,0,100,0.1,1,0.5,10)
        index <- index + 1
        paramTable[index,] <- c("micro_sd",1,place,0,100,0.1,1,0.5,10)
        index <- index + 1
        paramTable[index,] <- c("location_weight",1,place,0,1,0.1,1,0.1,1)
        index <- index + 1

        ## Related to forecasting
        paramTable[index,] <- c("incPropn2",1,place,0,1,0.1,1,0,10)
        index <- index + 1
        paramTable[index,] <- c("propn2",1,place,0,1,0.1,1,0,1)
        index <- index + 1
        
    }
    if(!is.null(saveFile)) write.table(paramTable,saveFile,row.names=FALSE,sep=",")
    return(paramTable)
}




#' Setup location specific parameter table
#'
#' Takes the entire data set and returns a useable format parameter table for the place specific parameters
#' @param locationDat the data frame of all data
#' @param locationParFile file name for the parameter table if exists
#' @return a parameter table
#' @export
setupLocationParTable <- function(locationDat, locationParFile="locationParams.csv"){
    if(file.exists(locationParFile)){
        print(paste("Reading in: ",locationParFile,sep=""))
        paramTable <- read.table(locationParFile, header=TRUE,sep=",",stringsAsFactors=FALSE)
    }
    else paramTable <- createLocationParTable(locationDat)
    return(paramTable)
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
#' @param peakTimings Optional data frame providing a specific upper and lower peak time bound for each state
#' @return a table of allowable parameters with columns for t0, mosquito density (R0), corresponding state (as this will vary by N_H and life expectancy), and corresponding peak time
#' @export
generateAllowableParams <- function(peakTime=927, peakTimeRange=60, stateNames,parTab=NULL,allowableParsFile="allowablePars.csv", R0max=6.5,peakTimings=NULL){
    if(!is.null(allowableParsFile) & file.exists(allowableParsFile)) allowablePars <- read.table(allowableParsFile)
    else {
        allowablePars <- NULL
        peakTimes <- matrix(nrow=100,ncol=100)
        for(local in stateNames){
            tmpTab <- parTab[parTab$local %in% c(local, "all"),]
            for(i in 1:nrow(peakTimes)){
                for(j in 1:ncol(peakTimes)){
                    if(!is.null(peakTimings)){
                        tmpPeakLower <- peakTimings[peakTimings$local==local,"start"]
                        tmpPeakUpper <- peakTimings[peakTimings$local==local,"end"]
                    } else {
                        tmpPeakLower <- peakTime - peakTimeRange/2
                        tmpPeakUpper <- peakTime + peakTimeRange/2
                    }
                    pars <- tmpTab$values
                    names(pars) <- tmpTab$names

                    pars["density"] <- j/10
                    pars["t0"] <- i*10
                    
                    y0s <- generate_y0s(as.numeric(pars["N_H"]),as.numeric(pars["density"]))
                    t_pars <- seq(0,3003,by=1)
                    y <- solveSEIRModel_rlsoda(t_pars, y0s,pars,TRUE)
                    peakTimes[i,j] <- y[which.max(diff(y[,"incidence"])),"time"]
                    R0 <- r0.calc(pars)
                    if(R0 > 1 & R0 < R0max & peakTimes[i,j] > tmpPeakLower & peakTimes[i,j] < tmpPeakUpper){
                        allowablePars <- rbind(allowablePars,data.frame(i*10,j/10,local,peakTimes[i,j],R0))
                    }
                }
            }
        }
        allowablePars <- allowablePars[complete.cases(allowablePars),]
        colnames(allowablePars) <- c("t0","density","local","peak","r0")
    }
    return(allowablePars)
}

