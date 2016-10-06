#' MCMC diagnostic tests
#'
#' Runs some basic MCMC diagnostics on the given chain and saves a few plots. The diagnostics are:
#' \itemize{
#' \item{Gelman Diagnostics: }{Saves the Gelman plots at the given file location}
#' \item{Auto-correlation: }{Saves autocorrelation plots at the given file location}
#' }
#' @param mcmc_chains the entire MCMC chain to be tested
#' @param filename the full file path at which to save the diagnostics. _gelman.pdf will be appended, for example
#' @param param_table the parameter table
#' @param VERBOSE boolean flag for additional output. Defaulst to FALSE
#' @return returns any error messages raised during the tests
mcmc_diagnostics <- function(mcmc_chains, filename, param_table,VERBOSE=FALSE){
    errors <- NULL
    final <- NULL
    if(length(mcmc_chains) > 1){
        gelman.error <- tryCatch({
            if(VERBOSE) print("Saving Gelman diagnostics")
            gelman.filename <- paste(filename,"_gelman.pdf",sep="")
            pdf(gelman.filename)
            gelman.plot(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            if(VERBOSE) print("Gelman diagnostics:")
            gelman.diag(as.mcmc.list(mcmc_chains)[,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war1){
            if(VERBOSE) print(paste("A warnings occured in gelman diagnostics: ", war1))
            final <- war1
        }, error=function(err1){
            if(VERBOSE) print(paste("An error occured in gelman diagnostics: ", err1))
            final <- err1
        }, finally = {
            dev.off()
            errors <- c(errors, final)
        })
    }

    for(i in 1:length(mcmc_chains)){
        autocorr.error <- tryCatch({
            if(VERBOSE) print("Saving auto-correlation plot")
            autocorr.filename <- paste(filename,"_autocor_",i,".pdf",sep="")
            pdf(autocorr.filename)
            autocorr.plot(mcmc_chains[[i]][,which(param_table$fixed==0)])
            final <- NULL
        }, warnings=function(war2){
            if(VERBOSE) print(paste("A warnings occured in autocorr plot: ", war2))
            final <- war2
        }, error=function(err2){
            if(VERBOSE) print(paste("An error occured in autocorr plot diagnostics: ", err2))
            final <- err2
        }, finally = {
            dev.off()
            errors <- c(errors, final)            
        })
    }
    
    return(errors)
}



#' @useDynLib zikaProj
make_optim_r0 <- function(t_pars, pars, incDat){
    function(x){
        pars["density"] <- x[1]
        pars["constSeed"] <- x[2]
        pars["incPropn"] <- x[3]
        pars["baselineInc"] <- x[4]
        y0s <- generate_y0s(pars["N_H"],pars["density"])

        ts <-
        
        y <- solveModelSimple(t_pars,y0s,pars)
        tmpY <- y[y[,"times"] >= min(incDat[,"startDay"]) & y[,"times"] <= max(incDat[,"endDay"]),]
        N_H <- average_buckets(rowSums(tmpY[,c("I_H","S_H","E_H","R_H")]), incDat[,"buckets"])
        inc <- diff(tmpY[,"incidence"])
        inc <- colSums(matrix(inc,nrow=(1/t_pars["step"]) * 7))
                                        #inc <- average_buckets(tmpY[,"I_H"], incDat[,"buckets"])
        perCapInc <- (1-(1-(inc/N_H))*(1-pars["baselineInc"]))*pars["incPropn"]
        return(-incidence_likelihood(perCapInc, incDat[,"inc"],incDat[,"N_H"]))
    }
}


#' @useDynLib zikaProj
make_optim_micro <- function(t_pars, pars, parNames,microDat){
    function(x){
        pars[names(pars) %in% parNames] <- x

        y0s <- generate_y0s(pars["N_H"],pars["density"])
        
        y <- solveModelSimple(t_pars,y0s,pars)
        probs <- generate_micro_curve(pars)
        probM <- generate_probM(y[,"I_M"], pars["N_H"], probs, pars["b"], pars["p_MH"], pars["baselineProb"], 1)*pars["propn"]
        probM <- probM[which(y[,"times"] >= min(microDat$startDay) & y[,"times"] <= max(microDat$endDay))]
        probM <- average_buckets(probM, microDat$buckets)
        
        return(-likelihood_probM(microDat$microCeph, microDat$births, probM))
    }
}

#' @useDynLib zikaProj
optimise <- function(state, parTab, t_pars, incDat, testDat){
    incPars <- c("density","constSeed","incPropn","baselineInc")
    microcephPars <- c("baselineProb","mean","var","c","p1","p2","p3","p4","p5","p6","p7","p8")
    checkLn <- NaN
    tmpDat <- incDat[incDat$local == state,]
    tmpTab <- parTab[parTab$local %in% c("all",state),]
    tmpMicro <- testDat[testDat$local == state,]
    indices <- parTab$local %in% c("all",state)
    while(is.nan(checkLn) | is.infinite(checkLn)){
        startPars <- generate_start_pars(tmpTab)
        if(!is.null(incDat)){
            lik_r0 <- make_optim_r0(t_pars, startPars, tmpDat)
            opti1 <- optim(startPars[names(startPars) %in% incPars],lik_r0,control=list("maxit"=5000))
            startPars[names(startPars) %in% incPars] <- opti1$par
        }
        lik_micro <- make_optim_micro(t_pars, startPars, startPars[names(startPars) %in% microcephPars], tmpMicro)
        opti2 <- optim(startPars[names(startPars) %in% microcephPars], lik_micro, control=list("maxit"=5000))
        startPars[names(startPars) %in% microcephPars] <- opti2$par
        startPars <- bounds_check(startPars, parTab)
        checkLn <- posterior(t_pars, startPars, parTab$names, parTab$local,testDat$startDay,testDat$endDay,testDat$buckets,testDat$microCeph,testDat$births,testDat$local,incDat)
    }
    tmpPars <- parTab$values
    tmpPars[indices] <- startPars
    names(tmpPars) <- parTab$names
    return(tmpPars)
}



#' @useDynLib zikaProj
generate_start_pars <- function(parTab){
    startPars <- parTab$values
    names(startPars) <- parTab$names
    for(i in which(parTab$fixed==0)){
        startPars[i] <- runif(1,parTab[i,"lower_bounds"],parTab[i,"upper_bounds"])
    }
    return(startPars)
}

#' @useDynLib zikaProj
bounds_check <- function(x, parTab){
    pars <- x
    for(i in 1:length(x)){
            if(x[i] < parTab[i,"lower_bounds"]) pars[i] <- parTab[i,"lower_bounds"]
            if(x[i] > parTab[i,"upper_bounds"]) pars[i] <- parTab[i,"upper_bounds"]
    }
    return(pars)
}





#' Solve ODE model with age classes
#'
#' Given a list of parameters as generated by \code{\link{setupListPars}}, solves the ODE model and returns a named data frame that can be passed to \code{\link{plot_dynamics}}
#' @param t_pars vector of two components - the total run time in days and the time step for the ODE solver
#' @param y0s the initial population sizes for the ODE model
#' @param pars a named vector with all of the necessary parameters to solve the ODE model. See \code{\link{setupParsODE}}
#' @return a data frame of the solved ODE model
#' @seealso \code{\link{solveModelSimple}}
solveModelComplex <- function(t_pars,y0s, pars){
    ## Run time of model
    time_length <- t_pars[1]
    ## Step size
    time_by <- t_pars[2]

    ## Burn in and delayed epidemic start
    burnin <- pars["burnin"]
    epiStart <- pars["epiStart"]

    ## Package ODE pars
    pars <- pars[c("L_M","D_EM","L_H","D_C","D_F","D_EH","D_IH","b","p_HM","p_MH","constSeed")]

    ## Temporarily remove initial seeding  for burn in period
    I0 <- y0s["I_A"]
    y0s[11] <- 0
    y0s[5] <- y0s["S_A"] + I0

    ## Create time series for ode solver
    if(burnin+epiStart > time_by){
        ts1 <- seq(0, burnin+epiStart,by=time_by)
        y1 <- ode(y0s,ts1,func="derivs",parms=pars,dllname="zikaProj",initfunc="initmod",hmax=1e-4,nout=0)
        if(is.null(y1) || nrow(y1) < 1) return("Error")
        
        colnames(y1) <- c("times","S_M","E_M","I_M","S_C","S_A","S_F","E_C","E_A","E_F","I_C","I_A","I_F","R_C","R_A","R_F")
        y0s2 <- y1[nrow(y1),2:ncol(y1)]
        y1[,"times"] <- y1[,"times"] - burnin
        y1 <- y1[y1[,"times"] >= 0,]
       
    } else {
        y0s2 <- y0s
    }

    ## Restore I0 to positive value
    y0s2["I_A"] <- I0
    y0s2["S_A"] <- y0s2["S_A"] - I0
    
    ts2 <- seq(epiStart, time_length,by=time_by)
    y <- ode(y0s2,ts2,func="derivs",parms=pars,dllname="zikaProj",initfunc="initmod",hmax=1e-4,nout=0)
    if(is.null(y)) return("Error")

    if(epiStart > 0){
        y <- rbind(y1,y)
    }
    
    #y <- as.data.frame(y)
    colnames(y) <- c("times","S_M","E_M","I_M","S_C","S_A","S_F","E_C","E_A","E_F","I_C","I_A","I_F","R_C","R_A","R_F")
    return(y)
    
}
