



#' Obsolete - the R implementation of the ODE model
#'
#' Set of ODEs for the zika SEIR model to be used by deSolve
#' @param t the vector of times to be solved over
#' @param y the compartment states
#' @param pars the parameters used to solve the ODE model
#' @return set of derivatives at the given time point
#' @export
#' @useDynLib zikaProj
zika.ode <-function(t, y,pars){
  L_M <- pars[1]
  D_EM <- pars[2]
   
  L_H <- pars[3]
  D_C <- pars[4]
  D_F <- pars[5]
  D_EH <- pars[6]
  D_IH <- pars[7]
  B_H <- L_H - D_C
 
  b <- pars[8]
  P_HM <- pars[9]
  P_MH <- pars[10]
  
  t_seed <- pars[11]
  I0 <- pars[12]
  epsilon <- 0.0001
  seed <- 0
  offset <- 0
  
  if(t >= t_seed){
      offset <- pars[13]
      if(t < t_seed + epsilon){
          seed <- I0
      }
  }
  
  S_M <- y[1]
  E_M <- y[2]
  I_M <- y[3]
  S_C <- y[4]
  S_A <- y[5]
  S_F <- y[6]
  
  R_C <- y[13]
  R_A <- y[14]
  R_F <- y[15]
  
  N_H <- sum(y[4:15])
  N_M <- sum(y[1:3])
  
  lambda_M <- b*(N_M/N_H)*P_HM*(I_C + I_A + I_F + offset + seed)/N_H
  lambda_H <- b*(N_M/N_H)*P_MH*I_M/N_H
  
  dS_M <- N_M/L_M - S_M/L_M - lambda_M*S_M
  dE_M <- lambda_M*S_M - E_M/L_M - E_M/D_EM
  dI_M <- E_M/D_EM - I_M/L_M
  
  dS_C <- N_H/L_H - S_C/L_H - S_C/D_C - lambda_H*S_C
  dS_A <- S_C/D_C - S_A/L_H - S_A/B_H + S_F/D_F - lambda_H*S_A
  dS_F <- S_A/B_H - S_F/D_F - S_F/L_H - lambda_H*S_F
  
  dE_C <- lambda_H*S_C - E_C/L_H           - E_C/D_EH           - E_C/D_C
  dE_A <- lambda_H*S_A - E_A/L_H - E_A/B_H - E_A/D_EH + E_F/D_F + E_C/D_C
  dE_F <- lambda_H*S_F - E_F/L_H + E_A/B_H - E_F/D_EH - E_F/D_F
  
  dI_C <- E_C/D_EH - I_C/L_H - I_C/D_IH           - I_C/D_C
  dI_A <- E_A/D_EH - I_A/L_H - I_A/D_IH - I_A/B_H + I_C/D_C + I_F/D_F
  dI_F <- E_F/D_EH - I_F/L_H - I_F/D_IH + I_A/B_H           - I_F/D_F 
  
  dR_C <- I_C/D_IH - R_C/L_H - R_C/D_C
  dR_A <- I_A/D_IH + R_C/D_C - R_A/L_H - R_A/B_H + R_F/D_F
  dR_F <- I_F/D_IH + R_A/B_H - R_F/D_F - R_F/L_H
  
  dIf_A <- I_F/D_F
  df_B <- S_F/D_F + E_F/D_F + R_F/D_F + I_F/D_F
  dE <- E_C/D_EH + E_A/D_EH + E_F/D_EH
  
  return(list(c(dS_M,dE_M,dI_M,dS_C,dS_A,dS_F,dE_C,dE_A,dE_F,dI_C,dI_A,dI_F,dR_C,dR_A,dR_F,dIf_A,df_B,dE)))
}


test_SEIR <- function(t,y,pars){
    L_M <- pars[1]
    D_EM <- pars[2]
    
    L_H <- pars[3]
    D_C <- pars[4]
    D_F <- pars[5]
    D_EH <- pars[6]
    D_IH <- pars[7]
    B_H <- L_H - D_C
    
    b <- pars[8]
    P_HM <- pars[9]
    P_MH <- pars[10]
    
    t_seed <- pars[11]
    I0 <- pars[12]

    S_M <- y[1]
    E_M <- y[2]
    I_M <- y[3]
    S_H <- y[4]
    E_H <- y[5]
    I_H <- y[6]
    R_H <- y[7]

    N_H <- S_H + E_H + I_H + R_H
    N_M <- S_M + E_M + I_M

    lambda_M <- b*P_HM*(I_H)/N_H
    lambda_H <- b*P_MH*I_M/N_H

    dS_M <- -lambda_M*S_M - S_M/L_M + N_M/L_M
    dE_M <- lambda_M*S_M - E_M/D_EM - E_M/L_M
    dI_M <- E_M/D_EM - I_M/L_M

    dS_H <- -lambda_H*S_H - S_H/L_H
    dE_H <- lambda_H*S_H  - E_H/D_EH - E_H/L_H
    dI_H <- E_H/D_EH - I_H/D_IH - I_H/L_H
    dR_H <- I_H/D_IH - R_H/L_H

    dDeaths <- S_M/L_M + E_M/L_M +  I_M/L_M
    dBirths <- N_M/L_M
    return(list(c(dS_M, dE_M,dI_M,dS_H,dE_H,dI_H,dR_H,dDeaths,dBirths)))
    
}


#' Simulates head circumference data from an SEIR vector borne model
#'
#' Takes a list of necessary arguments to solve a system of ODEs that generates head circumferences of new born infants
#' @param allPars A list of parameter vectors. First vector should be an array of time points over which to solve the SEIR model. Second vector should be the vector of starting conditions. Third vector should be all of the parameters.
#' @return a matrix of simulated head circumferences over time.
#' @export
#' @useDynLib zikaProj
zika.sim <- function(allPars,headMeasurements=TRUE,buckets=NULL){
    ## Get length and time step for ODE solver
    y <- solveModel(allPars)
    y <- as.data.frame(y)
    y0s <- allPars[[2]]
    pars <- allPars[[3]]
    
    sampFreq <- pars["sampFreq"]
    sampPropn <- pars["sampPropn"]
    mu_I <- pars["mu_I"]
    sd_I <- pars["sd_I"]
    mu_N <- pars["mu_N"]
    sd_N <- pars["sd_N"]
    probMicro <- pars["probMicro"]
    baselineProb <- pars["baselineProb"]
    lifeExpectancy <- pars["L_H"]
 
    daysPerYear <- nrow(y)/max(y$times)
    birthsPerYear <- sum(y0s[4:6])/lifeExpectancy
    birthsPerDay <- ceiling(birthsPerYear/daysPerYear)
    if(headMeasurements){
        if(!is.null(buckets)){
            alphas_I<- calculate_alphas_buckets(as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),probMicro,buckets)
        }
        else {
            alphas_I<- calculate_alphas(as.matrix(unname(y[,c("I_F","S_F","E_F","R_F")])),probMicro,sampFreq)
        }
        alphas_N <- 1 - alphas_I
    }
    else {
        if(!is.null(buckets)){
            alphas_I <- calculate_alphas_prob_buckets(as.matrix(unname(y[,c("times","I_F","S_F","E_F","R_F")])),probMicro, baselineProb, buckets)
        }
        else {
            alphas_I <- calculate_alphas_prob_sampfreq(as.matrix(unname(y[,c("I_F","S_F","E_F","R_F")])),probMicro, baselineProb, sampFreq)
        }
    }
    N <- sampPropn*birthsPerDay*sampFreq
    index <- 1
    all <- NULL
    while(index <= length(alphas_I)){
        if(!is.null(buckets)) N <- sampPropn*(birthsPerYear*(buckets[index,"end"]-buckets[index,"start"]))
        else N <- sampPropn*(birthsPerDay*sampFreq)

        if(headMeasurements){
            components <- sample(1:2,c(alphas_I[index],alphas_N[index]),size=N,replace=TRUE)
            mus <- c(mu_I,mu_N)
            sds <- c(sd_I,sd_N)
            
            distribution <- round(rnorm(n=N,mean=mus[components],sd=sds[components]),digits=1)
            all[[index]] <- distribution
        }
        else all[[index]] <- rbind(ceiling(alphas_I[index]*N), ceiling(N-(alphas_I[index]*N)))
        index <- index + 1
    }
    tmp <- as.matrix(rbind.fill(lapply(all,function(x) {as.data.frame(t(x))})))
    if(!headMeasurements | !is.null(buckets)) colnames(tmp) <- c("microCeph","births")
    return(tmp)
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
    meanPrior <- dnorm(pars["mean"],15,5,1)
    return(meanPrior)
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



#' SEIR dynamics plot
#'
#' Plots SEIR dynamics given a data frame of solved ODEs
#' @param y The data frame or matrix of times and population sizes
#' @param N_H human population size
#' @param N_M mosquito population size
#' @param file.name optional filename at which to save the plot. Must be a PNG. If NULL, does not open the png device.
#' @export
#' @useDynLib zikaProj
plot_dynamics <- function(y, N_H, N_M, file.name = NULL){
    y <- as.data.frame(y)
    n <- ncol(y)
    cols <- c("times","S_M","E_M","I_M","S_C","S_A","S_F","E_C","E_A","E_F","I_C","I_A","I_F","R_C","R_A","R_F")
    y <- y[,1:length(cols)]
    colnames(y) <- cols

    if(!is.null(file.name)){
        png(file.name)
    }
    par(mfrow=c(2,2))
    plot(y$S_M~y$times,col="green",ylim=c(0,N_M),type='l',main="Mosquito Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_M~y$times,col="red")
    lines(y$I_M~y$times,col="blue")

    plot(y$S_C~y$times,col="green",ylim=c(0,0.3*N_H),type='l',main="Children Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_C~y$times,col="red")
    lines(y$I_C~y$times,col="blue")
    lines(y$R_C~y$times,col="purple")

    plot(y$S_A~y$times,col="green",ylim=c(0,0.8*N_H),type='l',main="Adult Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_A~y$times,col="red")
    lines(y$I_A~y$times,col="blue")
    lines(y$R_A~y$times,col="purple")

    plot(y$S_F~y$times,col="green",ylim=c(0,0.004*N_H),type='l',main="First Trimester Dynamics",xlab="Time",ylab="Incidence")
    lines(y$E_F~y$times,col="red")
    lines(y$I_F~y$times,col="blue")
    lines(y$R_F~y$times,col="purple")
    if(!is.null(file.name)){
        dev.off()
    }
    
}

#' Head circumference heatmap
#'
#' Given a matrix or data frame of head sizes over time (rows represent sampling times), plots a heatmap showing distribution and mean head sizes over time.
#' @param dat matrix of head count data. Rows represent sampling times and columns represent individual measurements.
#' @return A ggplot object with the heatmap of head sizes over time. White line shows mean head size.
#' @export
plotDataHeatMap <- function(dat){
tmp <- createCounts(dat)
meanDat <- tmp[[2]]
tmp <- tmp[[1]]

plot <- ggplot(tmp) + geom_raster(aes(x=Day,y=Size,fill=Proportion),interpolate=FALSE) +
    geom_line(data=meanDat,aes(y=y,x=x),linetype=2,colour="white",size=1)+
    scale_fill_gradientn(colours=c("darkblue","red")) +
    scale_y_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Size),by=1),limits=c(19,50),labels=seq(0,max(tmp$Size),by=1))+
    scale_x_continuous(expand=c(0,0),breaks=seq(0,max(tmp$Day),by=max(tmp$Day)/15),labels=round(seq(0,max(tmp$Day),by=max(tmp$Day)/15),digits=0))+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor=element_blank(),
        panel.background = element_blank(),
        text=element_text(size=16,colour="gray20"),
        axis.line=element_line(colour="gray20"),
        axis.line.x = element_line(colour = "gray20"),
        axis.line.y=element_line(colour="gray20")
    )
return(plot)
}



#' Auxiliary setup
#'
#' Given the run name, sets up all of the parameters needed to produce all plots
#' @param runName the name of the runto be plotted. This should match the filename of the MCMC chains
#' @param runDir the full directory path where the MCMC chains are located
#' @param allDatFile local filename of actual data OPTIONAL
#' @param nchain number of chains run
#' @param mcmcPars named vector (burnin, adaptive, thin) to format the chain correctly
#' @param version version of the model used
#' @return a list of all the parameters needed for plotting
#' @export
#' @useDynLib zikaProj
plot_setup <- function(runName, runDir="~/net/home/zika/outputs", allDatFile = "~/net/home/zika/allDat20.04.16.csv",nchains=3, mcmcPars=c("burnin"=50000,"adaptive"=100000,"thin"=50),version=3){
    curWD <- getwd()
    setwd(runDir)
    realDat <- read.csv(allDatFile)
    parTab <- setupParTable(version=version,realDat=realDat)
    pars <- setupListPars(duration = 3000,N_H = 9000000,N_M=27000000,version=version)

    ## Get filenames of chains
    filenames <- NULL
    for(i in 1:nchains) filenames[[i]] <- paste(sprintf("%s_%d", runName, i),"B_chain.csv",sep="_")

    priors <- NULL
    peakTimes <- NULL
    
    if(runName == "simulation_single_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco")
        realDat <- testDat
    } else if(runName == "simulation_single_priors"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco")
        peakTimes <- data.frame("start"=400,"end"=500,"local"=as.character(unique(states[states != "all"])))
    } else if(runName == "simulation_multi_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
    } else if(runName=="simulation_multi_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_A_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia")
    } else if(runName=="real_A_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_A_propn"){
        states <- c("all","pernambuco","bahia")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
        parTab[parTab$names=="propn","fixed"] <- 0
    } else if(runName=="real_B_unconstrained"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
    } else if(runName=="real_B_prior"){
        parTab[parTab$names=="propn","fixed"] <- 1
        states <- c("all","pernambuco","bahia","saopaulo")
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName=="real_B_propn"){
        states <- c("all","pernambuco","bahia","saopaulo")
        parTab[parTab$names=="propn","fixed"] <- 0
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    } else if(runName == "real_3_all"){
        states <- c("all","pernambuco", "bahia","saopaulo", "paraiba", "maranhao", "ceara","sergipe","riodejaneiro","piaui","riograndedonorte","minasgerais", "matogrosso", "alagoas", "para", "acre", "goias", "espiritosanto","tocantins")
        parTab[parTab$names=="propn","fixed"] <- 0
        peakTimes <- data.frame("start"=460,"end"=560,"local"=as.character(unique(states[states != "all"])))
    }
    
    parTab <- parTab[parTab$local %in% states,]
    testDat <- generate_multiple_data(pars[[1]],parTab,NULL)
    if(grepl("simulation",runName)) realDat <- testDat
    places <- states[states != "all"]
    
    if(!is.null(peakTimes)) peakTimes$local <- as.character(peakTimes$local)

    realDat <- realDat[realDat$local %in% places,colnames(realDat) %in% colnames(testDat)]
    
    parTab[parTab$names=="b","fixed"] <- 1
    parTab[parTab$names=="mean","fixed"] <- 0
    parTab[parTab$names=="var","fixed"] <- 0
    parTab[parTab$names=="scale","fixed"] <- 0

    allChains <- NULL
    for(i in 1:length(filenames)){
        tmpChain <- read.csv(filenames[[i]])
        tmpChain <- tmpChain[tmpChain[,1] > (mcmcPars["adaptive"]+mcmcPars["burnin"]),]
        tmpChain <- cbind(tmpChain, i)
        allChains[[i]] <- tmpChain
    }
    setwd(curWD)
    return(list("chains"=allChains,"data"=realDat,"parTab"=parTab,"states"=places,"pars"=pars))
}



