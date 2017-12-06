#' BIC calculation
#' 
#' Given an MCMC chain melted by location and subset location, calculates the BIC for the given chain
#' @param chain the MCMC chain to be tested
#' @param location the subset location to take. Note that 'local' must be a column to subset rows by.
#' @param parTab the parameter table used for this chain
#' @param allDat the microcephaly data
#' @param incDat the incidence data
#' @return a single BIC value
#' @export
calculate_BIC <- function(chain, location, parTab, allDat, incDat){
  tmpDat <- allDat[allDat$local==location,]
  tmpInc <- incDat[allDat$local==location,]
  n <- nrow(tmpDat) + nrow(tmpInc)
  
  tmpChain <- chain[chain$location==location,]
  maxLik <- -max(tmpChain$lnlike)
  B <- length(parTab[parTab$fixed==0,"values"])*log(n)
  return(2*maxLik + B)
}

#' Deviance
#'
#' Calculates deviance of a vector with a given likelihood function
#' @param x a vector of parameters used to solve the likelihood function
#' @param likelihood a likelihood function taking just a vector of parameters
#' @return a single deviance value
#' @export
calc_deviance  <- function(x, likelihood){
    return(-2*(likelihood(x)))
}

#' Posterior mean
#'
#' Calculates posterior mean of a given MCMC chain
#' @param chain the MCMC chain with a lnlike column
#' @return the posterior mean
#' @export
calc_post_mean <- function(chain){
    return(-2*mean(chain$lnlike))
}

#' Posterior variance
#'
#' Calculates posterior variance from a given MCMC chain
#' @param chain the MCMC chain with a lnlike column
#' @return the posterior variance
#' @export
calc_post_var <- function(chain){
    meanBar <- calc_post_mean(chain)
    tmp <- 0
    for(i in 1:nrow(chain)){
        tmp <- tmp + (-2*chain[i,"lnlike"] - meanBar)^2
    }
  varBar <- (1/(nrow(chain)-1)) * sum(tmp)
  return(varBar)
}

#' DIC
#'
#' Calculates DIC from a given MCMC chain. Optionally can look at subset by location.
#' @param chain the MCMC chain with a lnlike column and all columns (number of columns will be used
#' @param location optionally, subset the chain by a given location
#' @return a single DIC value
#' @export
calculate_DIC <- function(chain,location=NULL){
    tmpChain <- chain
    if(!is.null(location)) tmpChain <- chain[chain$location==location,]
    DIC <- calc_post_mean(tmpChain) + calc_post_var(tmpChain)/2
    return(DIC)
}

#' Microcephaly risk range
#'
#' Getting microcephaly risk curve range from a chain for a given location
#' @param chain the MCMC chain
#' @param location the given location to subset by, if present.
#' @param runs number of samples to take
#' @param limit the risk limit to consider as being at risk
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a list with three vectors. 1: first day of significant risk, 2: last day of significant risk; 3: the number of days spent at risk
get_microceph_range <- function(chain,location=NULL,runs, limit=0.001,scale=1){
    ## If there is a location column and a location is specified, subset chain by this
    if(!is.null(chain$location) & !is.null(location)){
        chain <- chain[chain$location==location,]
    } 

    ## Get random samples from chain
    samples <- sample(nrow(chain),runs)
    index <- 1
    allProbs <- NULL
    lower_lim <- NULL
    upper_lim <- NULL
    range_lim <- NULL
    
    for(i in samples){
        ## Get pars from each row
        tmpPars <- get_index_pars(chain,i)
        tmpPars["tstep"] <- 1

        ## Generate microcephaly curve from this
        probs <- generate_micro_curve(tmpPars)/scale

        ## Get those days above the limit
        tmp <- which(probs > limit)
        if(length(tmp) > 0){
            ## Note that we are 0 indexed
            lower_lim[index] <- tmp[1] - 1
            upper_lim[index] <- tmp[length(tmp)] - 1
            range_lim[index] <- upper_lim[index] - lower_lim[index]
            index <- index+1
        }
    }
    return(list("lower_lim"=lower_lim,"upper_lim"=upper_lim,"range"=range_lim))
}

#' Lower microcephaly risk bound
#'
#' For a given MCMC chain, calculates the first day at risk
#' @param chain the MCMC chain to run over
#' @param limit the probability of microcephaly considered at risk
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a vector of first risk days
#' @export
#' @seealso \code{\link{get_microceph_range}}
get_lower_micro_bound <- function(chain,limit=0.001,scale=1){
  lower <- numeric(nrow(chain))
  for(i in 1:nrow(chain)){
    tmpPars <- as.numeric(chain[i,])
    names(tmpPars) <- colnames(chain)
    tmpPars["tstep"] <- 1
    probs <- generate_micro_curve(tmpPars)/scale
    tmp <- which(probs > limit)
    if(length(tmp) > 0){
      lower[i] <- tmp[1] -1
    }
    else lower[i] <- 0
  }
  return(lower)
}

#' Upper microbound
#'
#' For a given MCMC chain, calculates the last day at risk
#' @param chain the MCMC chain to run over
#' @param limit the probability of microcephaly considered at risk
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a vector of last risk days
#' @export
#' @seealso \code{\link{get_microceph_range}}
get_upper_micro_bound <- function(chain,limit=0.001,scale=1){
  lower <- numeric(nrow(chain))
  for(i in 1:nrow(chain)){
    tmpPars <- as.numeric(chain[i,])
    names(tmpPars) <- colnames(chain)
    tmpPars["tstep"] <- 1
    probs <- generate_micro_curve(tmpPars)/scale
    tmp <- which(probs > limit)
    if(length(tmp) > 0){
      lower[i] <- tmp[length(tmp)] -1
    }
    else lower[i] <- NA
  }
  return(lower)
}

#' Max micro risk
#'
#' For a given MCMC chain, calculates the maximum risk of microcephaly given infection
#' @param chain the MCMC chain to run over
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a vector of estimated maximum risk from the MCMC chain
#' @export
get_max_micro <- function(chain,scale){
  maxPs <- numeric(nrow(chain))
  for(i in 1:nrow(chain)){
      tmpPars <- get_index_pars(chain,i)
      tmpPars["tstep"] <- 1
      probs <- generate_micro_curve(tmpPars)
      maxPs[i] <- max(probs)
  }
  return(maxPs)
}

#' Max micro risk day
#'
#' For a given MCMC chain, calculates the day of maximum risk of microcephaly given infection
#' @param chain the MCMC chain to run over
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a vector of estimated maximum risk days from the MCMC chain
#' @export
get_max_micro_day <- function(chain,scale=1){
    maxWs <- numeric(nrow(chain))
    for(i in 1:nrow(chain)){
        tmpPars <- get_index_pars(chain,i)
        tmpPars["tstep"] <- 1
        probs <- generate_micro_curve(tmpPars)/scale
        maxWs[i] <- which.max(probs) - 1
    }
    return(maxWs)
}

#' Load MCMC chains
#'
#' From a working directory with MCMC chains (files ending chain.csv), reads these in and puts them into an MCMC list
#' @param location the path to the directory to read from
#' @param asList if TRUE, returns the MCMC chains as a list
#' @param convertMCMC if TRUE, converts the chains to mcmc objects (coda)
#' @param unfixed if TRUE, only includes columns for estimated parameters (not fixed parameters)
#' @param thin thins the read in MCMC chains by the specified amount (useful for speed when testing)
#' @param burnin discards the burn in period
#' @param unique_names if TRUE, gives column names as if read in by read.csv
#' @return the list or data frame of MCMC chains
#' @export
load_mcmc_chains <- function(location="",asList=FALSE, convertMCMC=FALSE,unfixed=TRUE, thin=1, burnin=100000,unique_names=TRUE){
    chains <- Sys.glob(file.path(location,"*_multivariate_chain.csv"))
    if(length(chains) < 1){
        message("Error - no chains found")
        return(NULL)
    }
    
    ## Read in the MCMC chains with fread for speed
    read_chains <- lapply(chains,data.table::fread,data.table=FALSE)
    
    ## Thin and remove burn in
    read_chains <- lapply(read_chains, function(x) x[seq(1,nrow(x),by=thin),])
    read_chains <- lapply(read_chains,function(x) x[x$sampno > burnin,])
    
    ## Use names of a read.csv chain
    if(unique_names){
        tmpChain <- read.csv(chains[1])
        for(i in 1:length(read_chains)){
            colnames(read_chains[[i]]) <- colnames(tmpChain)
        }
    }
    for(i in 1:length(read_chains)){
        read_chains[[i]]$chain <- i
    }
    ## Get the estimated parameters only
    if(unfixed){
        fixed <- read_inipars(location)$fixed
        read_chains <- lapply(read_chains, function(x) x[,c(which(fixed==0)+1,ncol(x))])
    }

    ## Convert to MCMC
    if(convertMCMC){
        read_chains <- lapply(read_chains,coda::as.mcmc)
        if(asList) read_chains <- coda::as.mcmc.list(read_chains)
    } else {
        if(!asList) read_chains <- do.call("rbind",read_chains)
    }
    return(read_chains)
}

#' Read partab
#'
#' Read in the initial parameter file, which is the same as the parameter table produced from zikaInfer
#' @param location the path to the directory to read from
#' @return the full parameter table
#' @export
read_inipars <- function(location=""){
    pars <- Sys.glob(file.path(location,"*inipars.csv"))
    if(length(pars) > 0){
        read_pars <- data.table::fread(pars[1],data.table=FALSE)
    } else {
        read_pars <- NULL
    }
    return(read_pars)
}

#' Attack rate simeq
#' 
#' Simeq cost function for calculating attack rate from R0
#' @param par the current estimated attack rate
#' @param R0 the value of R0 to be tested
#' @return the difference between the two sides of the final size calculation
#' @export
simeq <- function(par,R0){
  A <- par[1]
  f1 <- A - (1-exp(-R0*A))
  f1
}

#' Gamma mode
#'
#' Calculates the mode of a gamma distribution given the mean and variance
#' @param mean the gamma mean
#' @param var the gamma variance
#' @return the gamma mode
#' @export
calculate_gamma_mode <- function(mean, var){
  theta <- var/mean
  k <- mean/theta
  mode <- (k-1)*theta
}

#' Attack rate
#'
#' Calculates attack rate using final size equation and R0
#' @param R0 the value of R0 to test
#' @return the estimated attack rate
#' @export
calculate_AR <- function(r0){
    attacks <- sapply(r0, function(x) nleqslv::nleqslv(runif(1,0.5,0.8),simeq,R0=x)$x)
}


#' Post calculations
#'
#' Adds some parameter estimates for peripheral microcephaly risk parameters eg. mode. Note that this only applies to version 1 of the model. If not version 1, will return lots of zeros
#' @param chain the full MCMC chain to do post analysis on
#' @param version the version of the model run here
#' @param microceph_limit the probability of microcephaly given infection to use as the cut off
#' @param scale if a non 100 percent reporting proportion, need to divide by assumed proportion
#' @return an MCMC chain of the added analyses
#' @export
extra_microceph_calculations <- function(chain,version=1,microceph_limit=0.001,scale=1){
    if(version == 1) {
        ## Get break days as trimesters
        break1 <- 14*7
        break2 <- 14*7+ break1
        break3 <- break2 + 12*7

        
        ## Calculate microcephaly stats for each sample     
        allProbs <- matrix(nrow=nrow(chain),ncol=280)
        
        for(i in 1:nrow(chain)){
            pars <- get_index_pars(chain,i)
            probs <- allProbs[i,] <- generate_micro_curve(pars)/scale
        }      
        tr1 <- rowMeans(allProbs[,1:break1])
        tr2 <- rowMeans(allProbs[,break1:break2])
        tr3 <- rowMeans(allProbs[,break2:break3])
        
        mode <- calculate_gamma_mode(chain$mean,chain$var)
        maxWeek <- get_max_micro_day(chain,scale)
        maxP <- get_max_micro(chain,scale)
        lower <- get_lower_micro_bound(chain,microceph_limit,scale)
        upper <- get_upper_micro_bound(chain,microceph_limit,scale)
        
        range <- upper - lower
        
    } else {
        mode <- 0
        maxWeek <- 0
        maxP <- 0
        tr1 <- 0
        tr2 <- 0
        tr3 <- 0
        lower <- upper <- range <- 0
    }
    chain <- cbind(mode=mode,maxWeek=maxWeek,lower=lower,upper=upper,range=range,mode=mode,max=maxP,tr1=tr1,tr2=tr2,tr3=tr3)
    return(chain)
}

#' Combine MCMC chains
#'
#' Given a list of MCMC chains, rbinds them all to give one big chain
#' @param chains the list of chains
#' @return a single data frame
#' @export
combine_chains <- function(chains){
    return(do.call("rbind",chains))
}

#' Summarise MCMC chain
#'
#' Gives the summary statistics and quantiles of all columns from an MCMC chain
#' @param chain the MCMC chain
#' @return a data frame of MCMC summary stats
#' @export
summarise_chain <- function(chain){
  tmp <- summary(as.mcmc(chain))
  return(cbind(tmp$statistics,tmp$quantiles))
}

#' Isolate location from chain
#'
#' Given the name of a location and its corresponding index, returns only those parameters relevant to that location from the MCMC chain
#' @param local the name of the location
#' @param i the index of the location used for indexing the right column (for example, first location in the table is 0, second is 1 etc...)
#' @param chains the MCMC chain to extract from
#' @param parTab the parameter table used to generate this MCMC chain. This is needed to ensure that parameter names are correct
#' @param extraNames vector of additional names in the MCMC chain
#' @return an MCMC chain for just this location, with R0 and attack rate calculated
#' @export
isolate_location <- function(local, i, chains, parTab,extraNames=NULL){
    ## Create a vector of the parameter names associated with this location
    location_pars <- parTab[parTab$local==local,"names"]
    if(i >= 1) location_pars <- paste(location_pars,".",i,sep="")
    locationNames <- c(location_pars,parTab[parTab$local=="all","names"],extraNames)

    ## Create vector of non-indexed names to name the final chain
    blankNames <- parTab[parTab$local==local,"names"]
    blankNames <- c(blankNames,parTab[parTab$local=="all","names"],extraNames)
    
    ## Get only those parameters related to this location and rename the chain
    tmpChain <- chains[,locationNames]
    colnames(tmpChain) <- blankNames

    ## Calculate RO and attack rate
    r0 <- r0.vector(tmpChain)
    attackRate <- calculate_AR(r0)

    ## Create a new chain with only location specific parameters
    newChain <- cbind(tmpChain, R0=r0,AR=attackRate,location=local)
    return(newChain)
}


#' Cumulative true
#'
#' Gets number of consecutive values in a vector that are TRUE
#' @param x the vector to be tested
#' @return a vector counting number of successive TRUEs
#' @export
cumul_true <- function(x)  {
    rl <- rle(x)
    len <- rl$lengths
    v <- rl$values
    cumLen <- cumsum(len)
    z <- x
    ## replace the 0 at the end of each zero-block in z by the
    ## negative of the length of the preceding 1-block....
    iDrops <- c(0, diff(v)) < 0
    z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
    ## ... to ensure that the cumsum below does the right thing.
    ## We zap the cumsum with x so only the cumsums for the 1-blocks survive:
    x*cumsum(z)
}


#' Median and quantiles
#'
#' Calculates the median and 95% quantiles from a vector
#' @param x the vector
#' @return vector with 95% quantiles and median
#' @export
median.quantile <- function(x){
  out <- quantile(x, probs = c(0.025,0.5,0.975))
  names(out) <- c("ymin","y","ymax")
  return(out)
}
