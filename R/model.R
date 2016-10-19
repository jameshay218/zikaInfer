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

#' Generates y0s
#'
#' Generates initial values for the simple SEIR model given population size and mosquito density
#' @param N_H human population size
#' @param density number of mosquitoes per person
#' @return a vector of initial population sizes
#' @export
#' @useDynLib
generate_y0s <- function(N_H, density){
    N_M <- unname(N_H)*unname(density)
    S_M = 1*(N_M)
    E_M = 0
    I_M = 0

    S_H = unname(N_H) - 10
    E_H = 0
    I_H = 10
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
#' Microcephaly risk curve under gamma distribution
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v1 <- function(pars){
    mean <- pars["mean"]
    var <- pars["var"]

    scale <- var/mean
    shape <- mean/scale

    #shape <- pars["mean"]
    #scale <- pars["var"]
    
    #probs <- dgamma(0:39,shape=shape,scale=scale)*pars["c"]
    probs <- dgamma(0:279,shape=shape,scale=scale)*pars["c"]
    probs[probs > 1] <- 1
    #probs <- rep(probs, each=pars["tstep"])
    return(probs)
}

#' Microcephaly risk V1
#'
#' Microcephaly risk curve with 3 distinct periods
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v2 <- function(pars){
    probs <- c(rep(pars["p1"], 14),rep(pars["p2"],14),rep(pars["p3"],12))
    probs <- rep(probs, each=pars["tstep"])
    return(unname(probs))
}

#' Microcephaly risk V1
#'
#' Microcephaly risk curve with 6 distinct periods
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v3 <- function(pars){
    probs <- c(rep(pars["p1"], 7),rep(pars["p2"],7),rep(pars["p3"],7),rep(pars["p4"],7),rep(pars["p5"],7),rep(pars["p6"],5))
    probs <- rep(probs, each=pars["tstep"])
    return(unname(probs))
}

#' Microcephaly risk V1
#'
#' Microcephaly risk curve with 8 distinct periods
#' @param pars the model parameters
#' @return the vector of risks
#' @export
microceph_v4 <- function(pars){
    probs <- c(rep(pars["p1"], 5),rep(pars["p2"],5),rep(pars["p3"],5),rep(pars["p4"],5),rep(pars["p5"],5),rep(pars["p6"],5),rep(pars["p7"],5),rep(pars["p8"],5))
    probs <- rep(probs, each=pars["tstep"])
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
    rlsoda::rlsoda(y0s, ts, CsimpleSEIR_rich, pars, dllname="zikaProj", deSolve_compatible = compatible,return_time=TRUE,return_initial=TRUE,atol=1e-5,rtol=1e-5)
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
    
    if(version ==1){
        pars <- c("burnin"=burnin,"baselineProb"=baseline,pars,"mean"=mean,"var"=var,"c"=c, "tstep"=tstep)
    } else if(version ==2){
        pars <- c("burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3)
    } else if(version==3){
        pars <- c("burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6)
    } else if(version==4){
        pars <- c("burnin"=burnin,"baselineProb"=baseline,pars,"tstep"=tstep, "p1"=p1, "p2"=p2, "p3"=p3, "p4"=p4, "p5"=p5, "p6"=p6, "p7"=p7, "p8"=p8)  
    } else {
        pars <- c("sampFreq"=sampFreq,"sampPropn"=sampPropn,"mu_I"=mu_I,"sd_I"=sd_I,"mu_N"=mu_N,"sd_N"=sd_N,"probMicro"=probMicro,"baselineProb"=baseline,"burnin"=burnin,"epiStart"=epiStart,pars,"mean"=mean,"var"=var,"c"=c,"tstep"=tstep)
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


#' MCMC parameter table creation
#'
#' Sets up the parameter table for use in \code{\link{run_metropolis_MCMC}}
#' @param version the version of the model to use. See \code{\link{print_version_names}}
#' @param realDat if real data is available, adds parameters for the states given
#' @param saveFile if provided, writes the parameter table to the given filename
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, initial step sizes and whether or not the parameter should be fixed.
#' @export
#' @useDynLib zikaProj
createParTable <- function(version=1,realDat=NULL, saveFile=NULL){
    names <- c("sampFreq","sampPropn","mu_I","sd_I","mu_N","sd_N","probMicro","baselineProb","burnin","epiStart","L_M","D_EM","L_H","D_C","D_F","D_EH","D_IH","b","p_HM","p_MH","constSeed","mean","var","c","tstep","p1","p2","p3","p4","p5","p6","p7","p8")    
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
    paramTable[paramTable[,"names"]=="baselineProb",2:ncol(paramTable)] <- c(0.0002,"all",0,1,0.1,0,0,0.001)
    paramTable[paramTable[,"names"]=="burnin",2:ncol(paramTable)] <- c(0,"all",0,1000,0.1,1,0,1000)
    paramTable[paramTable[,"names"]=="epiStart",2:ncol(paramTable)] <- c(0,"all",0,700,0.1,1,0,700)
    paramTable[paramTable[,"names"]=="L_M",2:ncol(paramTable)] <- c(14,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="D_EM",2:ncol(paramTable)] <- c(10.5,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="L_H",2:ncol(paramTable)] <- c(365*70,"all",0,200*365,0.1,1,0,200*365)
    paramTable[paramTable[,"names"]=="D_C",2:ncol(paramTable)] <- c(365*18,"all",0,25*365,0.1,1,0,25*365)
    paramTable[paramTable[,"names"]=="D_F",2:ncol(paramTable)] <- c(0.75*365,"all",0,365,0.1,1,0,365)
    paramTable[paramTable[,"names"]=="D_EH",2:ncol(paramTable)] <- c(5.9,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="D_IH",2:ncol(paramTable)] <- c(5,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="b",2:ncol(paramTable)] <- c(0.25,"all",0,100,0.1,1,0,100)
    paramTable[paramTable[,"names"]=="p_HM",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="p_MH",2:ncol(paramTable)] <- c(0.5,"all",0,1,0.1,1,0,1)
    paramTable[paramTable[,"names"]=="constSeed",2:ncol(paramTable)] <- c(0,"all",0,1000,0.1,1,400,600)
    paramTable[paramTable[,"names"]=="mean",2:ncol(paramTable)] <- c(5,"all",0,100,0.1,0,8,20)
    paramTable[paramTable[,"names"]=="var",2:ncol(paramTable)] <- c(3,"all",0,100,0.1,0,1,5)
    paramTable[paramTable[,"names"]=="c",2:ncol(paramTable)] <- c(1,"all",0,100,0.1,0,0.1,0.5)
    paramTable[paramTable[,"names"]=="tstep",2:ncol(paramTable)] <- c(7,"all",0,100,0.1,1,0,7)
    paramTable[paramTable[,"names"]=="p1",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p2",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p3",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p4",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p5",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p6",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p7",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)
    paramTable[paramTable[,"names"]=="p8",2:ncol(paramTable)] <- c(0.1,"all",0,1,0.1,0,0.01,0.1)

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
    numberPars <- 18
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
        paramTable[index,] <- c("propn",1,place,0,1,0.1,1,0.01,0.2)
        index <- index + 1
        paramTable[index,] <- c("constSeed",400,place,0,2000,0.1,0,400,600)
        index <- index + 1
        paramTable[index,] <- c("mean",17,place,0,100,0.1,0,8,20)
        index <- index + 1
        paramTable[index,] <- c("var", 3, place, 0,1000,0.1,0,1,5)
        index <- index + 1
        paramTable[index,] <- c("c",0.15,place,0,100,0.1,0,0.1,0.3)
        index <- index + 1
        paramTable[index,] <- c("incPropn",0.0025,place,0,1,0.1,0,0.001,0.005)
        index <- index + 1
        paramTable[index,] <- c("baselineInc",0.0001,place,0,1,0.1,0,0.00005,0.0005)
        index <- index + 1
        paramTable[index,] <- c("p1",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p2",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p3",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p4",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p5",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p6",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p7",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
        paramTable[index,] <- c("p8",0.1,place,0,1,0.1,1,0.01,0.1)
        index <- index + 1
    }
    if(!is.null(saveFile)) write.table(paramTable,saveFile,row.names=FALSE,sep=",")
    return(paramTable)
}



#' MCMC parameter setup
#'
#' Sets up the parameter table for use in \code{\link{run_metropolis_MCMC}}
#' @param version which model version is being run (see \code{\link{print_version_names}})
#' @param realDat data frame of the real microcephaly data
#' @param sharedProb boolean indicating if all states share the same risk curve
#' @param parFile if already specified, where to find the parameter table file without state parameters
#' @param stateParFile as parFile, but with state specific parameters added
#' @return a matrix of needed settings for the MCMC algorithm. For each parameter, gives a name, lower and upper bounds, boolean for log scale, initial step sizes, log proposal and whether or not the parameter should be fixed.
#' @export
setupParTable <- function(version=1, realDat=NULL, sharedProb=FALSE, parFile = "", stateParFile = ""){
    ## Checks if given filename exists. If not, creates a fresh parameter table
    if(file.exists(parFile)){
        print(paste("Reading in: ",parFile,sep=""))
        paramTable <- read.table(parFile, header=TRUE,sep=",",stringsAsFactors=FALSE)
    }
    else paramTable <- createParTable(version,realDat)
    useNames <- names(setupParsLong(version))
    
    ## Gets the parameter names used for each version of the model
    allMicroPars <- c("mean","var","c","p1","p2","p3","p4","p5","p6","p7","p8")
    if(version==2) microPars <- c("p1","p2","p3")
    else if(version==3) microPars <- c("p1","p2","p3","p4","p5","p6")
    else if(version==4) microPars <- c("p1","p2","p3","p4","p5","p6","p7","p8")
    else microPars <- c("mean","var","c")
    notMicroPars <- allMicroPars[!(allMicroPars %in% microPars)]
    
    ## Remove all unneccessary parameters
    paramTable <- paramTable[paramTable[,"names"] %in% useNames,]
    if(!sharedProb) paramTable <- paramTable[!(paramTable[,"names"] %in% microPars),]

    ## Adds state specific parameters
    stateParTable <- NULL
    if(!is.null(realDat) | file.exists(stateParFile)){
        stateParTable <- setupStateParTable(realDat, stateParFile)
        stateParTable <- stateParTable[!(stateParTable[,"names"] %in% paramTable$names),]
        stateParTable <- stateParTable[!(stateParTable[,"names"] %in% notMicroPars),]
        if(sharedProb) stateParTable <- stateParTable[!(stateParTable[,"names"] %in% microPars),]
    }

    paramTable <- rbind(paramTable, stateParTable)
    
    paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")] <- lapply(paramTable[,c("values","lower_bounds","upper_bounds","steps","fixed","start_lower","start_upper")], FUN=as.numeric)
    
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
    peakTimes <- character(length(unique_states))
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
        peakTime <- y[which.max(y[,"incidence"]),1]
        message(peakTime)
        peakTimes[place=peakTime]
    }
    peakTimes
    
}
