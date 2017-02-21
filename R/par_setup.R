#' Parameter setup
#'
#' Sets up the parameter table based on the version, states and allowable parameters
#' @param stateNames the vector of states to include
#' @param version the version of the model to run (1-3)
#' @param realDat the data frame of actual microcephaly data. Need this to generate N_H and L_H
#' @param useInc bool indicating whether or not ZIKV incidence data is considered here
#' @param allowablePars data frame of allowable parameters for density and constSeed
#' @param sharedProb bool if the version with shared microcephaly parameters should be used
#' @param normLik Using normal likelihood function for data? If false, uses binomial likelihood. If true, uses normal likelihood which incorporates two standard deviation parameters for ZIKV and microcephaly incidence data
#' @param stateWeights bool deciding if data from different states should be weighted
#' @return the modified parameter table
#' @export
partab_setup <- function(stateNames, version, realDat, useInc,allowablePars=NULL,sharedProb=TRUE, normLik=NULL,stateWeights=FALSE){
    correct_order <- c("pernambuco", "sergipe", "paraiba", "bahia", "riograndedonorte", 
                       "acre", "piaui", "alagoas", "matogrosso", "maranhao", "rondonia", 
                       "ceara", "tocantins", "espiritosanto", "roraima", "riodejaneiro", 
                       "saopaulo", "para", "minasgerais", "matogrossdosul", "distritofederal", 
                       "goias", "riograndedosul", "amazonas", "parana", "santacatarina", 
                       "amapa")
    
    ## Get first state from the correct order found in passed states
    norm_state <- correct_order[correct_order %in% stateNames][1]
    message(cat("Reference state: ",norm_state,"\n",sep=""))
    
    parTab <- setupParTable(version,realDat,sharedProb=sharedProb,stateNormLik=normLik,stateWeights=stateWeights)
    if(stateWeights) parTab[parTab$names=="state_weight","values"] <- 1/length(stateNames)
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
        parTab[parTab$names == "incPropn","values"] <- runif(length(parTab[parTab$names == "incPropn","values"]),0.0001,0.001)
        parTab[parTab$names=="inc_weight","values"] <- 0.5
        if(!is.null(normLik)){
            parTab[parTab$names=="inc_sd","fixed"] <- 0
        }
    }
    ## Arbritrarily fix one state's reporting proportion
    parTab[parTab$names=="propn" & parTab$local== norm_state,c("fixed","values")] <- c(1,0.8)
    ## If there's only one state, fix propn, as this is the same as 'c'
    if(length(stateNames) == 1) parTab[parTab$names=="propn" & parTab$local== stateNames[1],c("fixed","values")] <- c(1,1)

    ## If using normal likelihood function, need to allow SD to vary.
    if(!is.null(normLik)) parTab[parTab$names=="micro_sd","fixed"] <- 0
        
      ## If there's a data frame of allowable starting parameters, use these to seed density and constSeed
    if(!is.null(allowablePars)){
        for(state in stateNames){
            tmpTab <- parTab[parTab$local %in% c("all",state),c("values","names")]
            pars <- tmpTab$values
            names(pars) <- tmpTab$names

            ## Limit density to give an R0 between 1 and 10 - anything above around 6 just gives the same likelihood
            min_density <- density.calc(pars,1)
            max_density <- density.calc(pars,10)
            parTab[parTab$local == state & parTab$names == "density",c("lower_bounds","start_lower")] <- min_density
            parTab[parTab$local == state & parTab$names == "density",c("upper_bounds","start_upper")] <- max_density
                
            tmpAllowable <- allowablePars[allowablePars$local == state,]
            index <- as.integer(runif(1,1,nrow(tmpAllowable)))
            parTab[parTab$names == "constSeed" & parTab$local == state,"values"] <- as.numeric(tmpAllowable[index,"constSeed"])
            parTab[parTab$names == "density" & parTab$local == state,"values"] <- as.numeric(tmpAllowable[index,"density"])
        }
    }
    parTab <- parTab[parTab$local %in% c("all",stateNames),]
    return(parTab)
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
    
    sampFreq <- 7
    sampPropn <- 0.9
    mu_I <- 28
    sd_I <- 2
    mu_N <- 35
    sd_N <- 2
    probMicro <- 1
    baseline <- 0.00002
    burnin <- 3
    epiStart <- 0
    mean <- 17
    var <- 3
    c <- 0.15
    tstep <- 7
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
        pars <- c("state_weight"=1,"inc_weight"=inc_weight,"burnin"=burnin,"baselineProb"=baseline,pars,"mean"=mean,"var"=var,"c"=c, "tstep"=tstep)
    } else if(version ==2){
        pars <- c("state_weight"=1,"inc_weight"=inc_weight,"burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3)
    } else if(version==3){
        pars <- c("state_weight"=1,"inc_weight"=inc_weight,"burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6)
    } else if(version==4){
        pars <- c("state_weight"=1,"inc_weight"=inc_weight,"burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6, "p7"=p7, "p8"=p8)  
    } else {
        pars <- c("state_weight"=1,"inc_weight"=inc_weight,"sampFreq"=sampFreq,"sampPropn"=sampPropn,"mu_I"=mu_I,"sd_I"=sd_I,"mu_N"=mu_N,"sd_N"=sd_N,"probMicro"=probMicro,"baselineProb"=baseline,"burnin"=burnin,"epiStart"=epiStart,pars,"mean"=mean,"var"=var,"c"=c,"tstep"=tstep)
    }
    return(pars)   
}

#' Default parameter list setup
#'
#' Sets up the default list of parameters and initial conditions to be used by \code{\link{generate_multiple_data}}. Can specify which version of the model to use. Passing 1 gives the simple model with no age classes, otherwise gives the parameters for the model with age classes.
#' @param duration Duration in years of the simulation. Defaults to 10 years
#' @param version the version of the model that is going to be used (see \code{\link{print_version_names}})
#' @return list of times, initial conditions, and ODE parameters
#' @export
#' @seealso \code{\link{setupParsLong}}
#' @useDynLib zikaProj
setupListPars <- function(duration=2*365, version=1){
    pars <- setupParsLong(version)
    y0 <- generate_y0s(pars["N_H"], pars["density"])
    ts <- seq(0,duration,by=1)
    return(list(ts,y0,pars))
}


#' MCMC parameter setup
#'
#' Sets up the parameter table for use in \code{\link{run_metropolis_MCMC}}
#' @param version which model version is being run (see \code{\link{print_version_names}})
#' @param realDat data frame of the real microcephaly data
#' @param sharedProb boolean indicating if all states share the same risk curve
#' @param parFile if already specified, where to find the parameter table file without state parameters
#' @param stateParFile as parFile, but with state specific parameters added
#' @param stateNormLik boolean deciding whether or not a normal likelihood function is used (otherwise binomial). If NULL, uses the binomial. If not NULL, true or false indicates whether the likelihood standard deviations are state specific or not.
#' @param stateWeights bool deciding if state weightings should be used (otherwise all states have equal weighting)
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, boolean for log scale, initial step sizes, log proposal and whether or not the parameter should be fixed.
#' @export
setupParTable <- function(version=1, realDat=NULL, sharedProb=FALSE, parFile = "", stateParFile = "",stateNormLik=NULL,stateWeights=FALSE){
    ## Checks if given filename exists. If not, creates a fresh parameter table
    if(file.exists(parFile)){
        print(paste("Reading in: ",parFile,sep=""))
        paramTable <- read.table(parFile, header=TRUE,sep=",",stringsAsFactors=FALSE)
    }
    else paramTable <- createParTable(NULL)
  
    useNames <- names(setupParsLong(version))
    if(stateWeights) useNames <- useNames[useNames != "state_weight"]
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
    if(!is.null(stateNormLik)){
        if(!stateNormLik) {
            useNames <- c(useNames, "micro_sd","inc_sd")
        }
    }
    
    ## Remove all unneccessary parameters
    paramTable <- paramTable[paramTable[,"names"] %in% useNames,]
    if(!sharedProb) paramTable <- paramTable[!(paramTable[,"names"] %in% microPars),]

    ## Adds state specific parameters
    stateParTable <- NULL
    if(!is.null(realDat) | file.exists(stateParFile)){
        stateParTable <- setupStateParTable(realDat, stateParFile)
        stateParTable <- stateParTable[!(stateParTable[,"names"] %in% paramTable$names),]
        stateParTable <- stateParTable[!(stateParTable[,"names"] %in% notMicroPars),]
        if(is.null(stateNormLik)) stateParTable <- stateParTable[!(stateParTable$names %in% c("micro_sd","inc_sd")),]
        if(sharedProb) stateParTable <- stateParTable[!(stateParTable[,"names"] %in% microPars),]
    }

    paramTable <- rbind(paramTable, stateParTable)
    
    paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")] <- lapply(paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")], FUN=as.numeric)
    
    return(paramTable)
}

#' MCMC parameter table creation
#'
#' Sets up the parameter table for use in \code{\link{run_metropolis_MCMC}}
#' @param saveFile if provided, writes the parameter table to the given filename
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, initial step sizes and whether or not the parameter should be fixed.
#' @export
#' @useDynLib zikaProj
createParTable <- function(saveFile=NULL){
    names <- c("sampFreq","sampPropn","mu_I","sd_I","mu_N","sd_N","probMicro","inc_weight","baselineProb","burnin","epiStart","L_M","D_EM","L_H","D_C","D_F","D_EH","D_IH","b","p_HM","p_MH","constSeed","mean","var","c","tstep","p1","p2","p3","p4","p5","p6","p7","p8","inc_sd","micro_sd","state_weight")    
    paramTable <- matrix(0, ncol=9, nrow=length(names))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")
    paramTable[,"names"] <- names
    paramTable$names <- as.character(paramTable$names)
    

    paramTable[paramTable[,"names"]=="sampFreq",2:ncol(paramTable)] <- c(7,"all",0,30,0.1,1,0,30)
    paramTable[paramTable[,"names"]=="sampPropn",2:ncol(paramTable)] <- c(1,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="mu_I",2:ncol(paramTable)] <- c(30,"all",0,60,0.1,1,0,60)
    paramTable[paramTable[,"names"]=="sd_I",2:ncol(paramTable)] <- c(2,"all",0,10,0.1,1,0,10)
    paramTable[paramTable[,"names"]=="mu_N",2:ncol(paramTable)] <- c(30,"all",0,60,0.1,1,0,60)
    paramTable[paramTable[,"names"]=="sd_N",2:ncol(paramTable)] <- c(2,"all",0,10,0.1,1,0,10)
    paramTable[paramTable[,"names"]=="probMicro",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="inc_weight",2:ncol(paramTable)] <- c(0,"all",0,1,0,1,0,1)
    paramTable[paramTable[,"names"]=="baselineProb",2:ncol(paramTable)] <- c(0.0002,"all",0,1,0.1,0,0,0.001)
    paramTable[paramTable[,"names"]=="burnin",2:ncol(paramTable)] <- c(0,"all",0,1000,0.1,1,0,1000)
    paramTable[paramTable[,"names"]=="epiStart",2:ncol(paramTable)] <- c(0,"all",0,700,0.1,1,0,700)
    paramTable[paramTable[,"names"]=="L_M",2:ncol(paramTable)] <- c(5,"all",2,20,0.1,1,2,14)
    paramTable[paramTable[,"names"]=="D_EM",2:ncol(paramTable)] <- c(8.4,"all",7.3,9.3,0.1,1,7.3,9.3)
    paramTable[paramTable[,"names"]=="L_H",2:ncol(paramTable)] <- c(365*70,"all",0,200*365,0.1,1,0,200*365)
    paramTable[paramTable[,"names"]=="D_C",2:ncol(paramTable)] <- c(365*18,"all",0,25*365,0.1,1,0,25*365)
    paramTable[paramTable[,"names"]=="D_F",2:ncol(paramTable)] <- c(0.75*365,"all",0,365,0.1,1,0,365)
    paramTable[paramTable[,"names"]=="D_EH",2:ncol(paramTable)] <- c(4,"all",2.3,6,0.1,1,2.3,6)
    paramTable[paramTable[,"names"]=="D_IH",2:ncol(paramTable)] <- c(6,"all",4,8,0.1,1,4,8)
    paramTable[paramTable[,"names"]=="b",2:ncol(paramTable)] <- c(0.5,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="p_HM",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="p_MH",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="constSeed",2:ncol(paramTable)] <- c(0,"all",0,1000,0.1,1,400,600)
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
    paramTable[paramTable[,"names"]=="state_weight",2:ncol(paramTable)] <- c(1,"all",0,1,0.1,1,0.1,1)
    if(!is.null(saveFile)) write.table(paramTable,saveFile,row.names=FALSE,sep=",")
    
    return(paramTable)
}

#' MCMC state-specific parameter table creation
#'
#' Sets up the parameter table for all states for use in \code{\link{run_metropolis_MCMC}}
#' @param stateDat data frame of states - had L_H, N_H and local name
#' @param saveFile if provided, writes the parameter table to the given filename
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, initial step sizes and whether or not the parameter should be fixed.
#' @export
#' @useDynLib zikaProj
createStateParTable <- function(stateDat, saveFile = NULL){
    places <- as.character(unique(stateDat$local))
    numberPars <- 21
    paramTable <- matrix(0, ncol=9, nrow=numberPars*length(places))
    paramTable <- as.data.frame(paramTable)
    colnames(paramTable) <- c("names", "values","local", "lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")
    index <- 1
    
    for(place in places){
        tmpDat <- stateDat[stateDat$local==place,]
        paramTable[index,] <- c("L_H",tmpDat[1,"L_H"]*365,place,0,200*365,0.1,1,0,200*365)
        index <- index + 1
        paramTable[index,] <- c("N_H", tmpDat[1,"N_H"], place, 0,100000000,0.1,1,0,1000000000)
        index <- index + 1
        paramTable[index,] <- c("density",3,place,0,100,0.1,0,3,6)
        index <- index + 1
        paramTable[index,] <- c("propn",1,place,0,1,0.1,0,0.01,0.2)
        index <- index + 1
        paramTable[index,] <- c("constSeed",400,place,0,2000,0.1,0,400,600)
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
        paramTable[index,] <- c("state_weight",1,place,0,1,0.1,1,0.1,1)
        index <- index + 1
        
    }
    if(!is.null(saveFile)) write.table(paramTable,saveFile,row.names=FALSE,sep=",")
    return(paramTable)
}




#' Setup state specific parameter table
#'
#' Takes the entire data set and returns a useable format parameter table for the place specific parameters
#' @param stateDat the data frame of all data
#' @param stateParFile file name for the parameter table if exists
#' @return a parameter table
#' @export
setupStateParTable <- function(stateDat, stateParFile="stateParams.csv"){
    if(file.exists(stateParFile)){
        print(paste("Reading in: ",stateParFile,sep=""))
        paramTable <- read.table(stateParFile, header=TRUE,sep=",",stringsAsFactors=FALSE)
    }
    else paramTable <- createStateParTable(stateDat)
    return(paramTable)
}
