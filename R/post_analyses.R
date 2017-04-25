#' BIC calculation
#' 
#' Given an MCMC chain melted by state and subset state, calculates the BIC for the given chain
#' @param chain the MCMC chain to be tested
#' @param state the subset state to take. Note that 'local' must be a column to subset rows by.
#' @param parTab the parameter table used for this chain
#' @param allDat the microcephaly data
#' @param incDat the incidence data
#' @return a single BIC value
#' @export
calculate_BIC <- function(chain, state, parTab, allDat, incDat){
  tmpDat <- allDat[allDat$local==state,]
  tmpInc <- incDat[allDat$local==state,]
  n <- nrow(tmpDat) + nrow(tmpInc)
  
  tmpChain <- chain[chain$state==state,]
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
#' Calculates DIC from a given MCMC chain. Optionally can look at subset by state.
#' @param chain the MCMC chain with a lnlike column and all columns (number of columns will be used
#' @param state optionally, subset the chain by a given state
#' @return a single DIC value
#' @export
calculate_DIC <- function(chain,state=NULL){
    tmpChain <- chain
    if(!is.null(state)) tmpChain <- chain[chain$state==state,]
    
    DIC <- calc_post_mean(tmpChain) + calc_post_var(tmpChain)/2
    return(DIC)
}

#' Microcephaly risk range
#'
#' Getting microcephaly risk curve range from a chain for a given state
#' @param chain the MCMC chain
#' @param state the given state to subset by, if present.
#' @param runs number of samples to take
#' @param limit the risk limit to consider as being at risk
#' @param scale scales microcephaly risk curve by assumed reporting proportion
#' @return a list with three vectors. 1: first day of significant risk, 2: last day of significant risk; 3: the number of days spent at risk
get_microceph_range <- function(chain,state=NULL,runs, limit=0.001,scale=1){
    ## If there is a state column and a state is specified, subset chain by this
    if(!is.null(chain$state) & !is.null(state)){
        chain <- chain[chain$state==state,]
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
    return(list(lower_lim,upper_lim,range_lim))
}

#' Lower microbound
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
#' From a working directory with MCMC chains (files ending _slice_chain.csv), reads these in and puts them into an MCMC list
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
    chains <- Sys.glob(file.path(location,"*_slice_chain.csv"))
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
#' Read in the initial parameter file, which is the same as the parameter table produced from zikaProj
#' @param location the path to the directory to read from
#' @return the full parameter table
#' @export
read_inipars <- function(location=""){
    pars <- Sys.glob("*inipars.csv")
  read_pars <- fread(pars[1],data.table=FALSE)       
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


#' State conversion
#'
#' Contains a vector of the brazilian state names to convert lower-case simplified names to actual strings
#' @param name the vector or scalar of lower case state names
#' @return the same vector but with actual strings
#' @export
convert_name_to_state <- function(name){
  country_names <- c("pernambuco"="Pernambuco","amapa"="Amapá","amazonas" ="Amazonas","distritofederal" = "Distrito Federal",
                     "bahia"="Bahia","saopaulo"="São Paulo","paraiba"="Paraíba","maranhao"="Maranhão","ceara"="Ceará",
                     "sergipe"="Sergipe","riodejaneiro"="Rio de Janeiro","piaui"="Piauí","riograndedonorte"="Rio Grande Norte",
                     "minasgerais"="Minas Gerais", "matogrosso"="Mato Grosso","alagoas"="Alagoas","para"="Pará","acre"="Acre",
                     "espiritosanto"="Espírito Santo","goias"="Goiás","tocantins"="Tocantins","matogrossodosul"="Mato Grosso do Sul",
                     "matogrossdosul"="Mato Grosso do Sul","parana"="Paraná","riograndedosul"="Rio Grande do Sul","rondonia"="Rondônia",
                     "roraima"="Roraima","santacatarina"="Santa Catarina", "colombia"="Colombia","northeast"="Northeast Brazil")
  return(country_names[name])
}

#' Post calculations
#'
#' Adds some parameter estimates for peripheral microcephaly risk parameters eg. mode. Note that this only applies to version 1 of the model; if not version 1, will return lots of zeros
#' @param chain the full MCMC chain to do post analysis on
#' @param version the version of the model run here
#' @param microceph_limit the probability of microcephaly given infection to use as the cut off
#' @param scale if a non 100% reporting proportion, need to divide by assumed proportion
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


#' Get states lowercase
#'
#' Gives a vector of all states in Brazil. Names are actual names, values are lower case
#' @return vector of all states
#' @export
get_states <- function(){
    states <- c(
    "Acre"="acre",
    "Alagoas"="alagoas",
    "Amapá"="amapa",
    "Amazonas"="amazonas",
    "Bahia"="bahia",
    "Ceará"="ceara",
    "Distrito Federal"="distritofederal",
    "Espírito Santo"="espiritosanto",      
    "Goiás"="goias",
    "Maranhão"="maranhao",
    "Mato Grosso do Sul"="matogrossdosul",
    "Mato Grosso"="matogrosso",
    "Minas Gerais"="minasgerais",
    "Pará"="para",
    "Paraíba"="paraiba",
    "Paraná"="parana",
    "Pernambuco"="pernambuco",
    "Piauí"="piaui",
    "Rio de Janeiro"="riodejaneiro",
    "Rio Grande do Norte"="riograndedonorte",
    "Rio Grande do Sul"="riograndedosul",
    "Rondônia"="rondonia",
    "Roraima"="roraima",
    "São Paulo"="saopaulo",
    "Santa Catarina"="santacatarina",
    "Sergipe"="sergipe",
    "Tocantins"="tocantins",
    "Colombia"="colombia",
    "Northeast Brazil"="northeast"
  )
  return(states)
}


#' Get states uppercase
#'
#' Gives a vector of all states in Brazil. Names are lowercase names, values are actual names
#' @return vector of all states
#' @export
states_to_title <- function(){
  states <- c(
      "acre"="Acre",
      "alagoas"="Alagoas",
      "amapa"="Amapá",
      "amazonas"="Amazonas",
      "bahia"="Bahia",
      "ceara"="Ceará",
      "distritofederal"="Distrito Federal",
      "espiritosanto"= "Espírito Santo",
      "goias"="Goiás",
      "maranhao"="Maranhão",
      "matogrossdosul"="Mato Grosso do Sul",
      "matogrosso"="Mato Grosso",
      "minasgerais"="Minas Gerais",
      "para"="Pará",
      "paraiba"="Paraíba",
      "parana"="Paraná",
      "pernambuco"="Pernambuco",
      "piaui"="Piauí",
      "riodejaneiro"= "Rio de Janeiro",
      "riograndedonorte"="Rio Grande do Norte",
      "riograndedosul"= "Rio Grande do Sul",
      "rondonia"= "Rondônia",
      "roraima"="Roraima",
      "saopaulo"="São Paulo",
      "santacatarina"="Santa Catarina",
      "sergipe"="Sergipe",
      "tocantins"="Tocantins",
      "colombia"="Colombia",
      "northeast"="Northeast Brazil"
  )
  return(states)
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

#' Isolate state from chain
#'
#' Given the name of a state and its corresponding index, returns only those parameters relevant to that state from the MCMC chain
#' @param local the name of the state
#' @param i the index of the state used for indexing the right column (for example, first state in the table is 0, second is 1 etc...)
#' @param chains the MCMC chain to extract from
#' @param parTab the parameter table used to generate this MCMC chain. This is needed to ensure that parameter names are correct
#' @param extraNames vector of additional names in the MCMC chain
#' @return an MCMC chain for just this state, with R0 and attack rate calculated
#' @export
isolate_state <- function(local, i, chains, parTab,extraNames=NULL){
    ## Create a vector of the parameter names associated with this state
    state_pars <- parTab[parTab$local==local,"names"]
    if(i >= 1) state_pars <- paste(state_pars,".",i,sep="")
    stateNames <- c(state_pars,parTab[parTab$local=="all","names"],extraNames)

    ## Create vector of non-indexed names to name the final chain
    blankNames <- parTab[parTab$local==local,"names"]
    blankNames <- c(blankNames,parTab[parTab$local=="all","names"],extraNames)
    
    ## Get only those parameters related to this state and rename the chain
    tmpChain <- chains[,stateNames]
    colnames(tmpChain) <- blankNames

    ## Calculate RO and attack rate
    r0 <- r0.vector(tmpChain)
    attackRate <- calculate_AR(r0)

    ## Create a new chain with only state specific parameters
    newChain <- cbind(tmpChain, R0=r0,AR=attackRate,state=local)
    return(newChain)
}

#' Convert name to factor
#'
#' Function to convert state strings to factors (with correct accents etc)
#' @param states the vector of lower case state names to convert
#' @return the converted states with proper naming and factored
#' @export
convert_name_to_state_factor <- function(states){
    order <- get_correct_order(normalised=TRUE)
    country_names <- c("pernambuco"="Pernambuco","amapa"="Amapá",
                       "amazonas" ="Amazonas","distritofederal" = "Distrito Federal","bahia"="Bahia",
                       "saopaulo"="São Paulo","paraiba"="Paraíba","maranhao"="Maranhão","ceara"="Ceará",
                       "sergipe"="Sergipe","riodejaneiro"="Rio de Janeiro","piaui"="Piauí",
                       "riograndedonorte"="Rio Grande Norte","minasgerais"="Minas Gerais",
                       "matogrosso"="Mato Grosso","alagoas"="Alagoas","para"="Pará","acre"="Acre",
                     "espiritosanto"="Espírito Santo","goias"="Goiás","tocantins"="Tocantins","matogrossodosul"="Mato Grosso do Sul","matogrossdosul"="Mato Grosso do Sul","parana"="Paraná","riograndedosul"="Rio Grande do Sul","rondonia"="Rondônia","roraima"="Roraima","santacatarina"="Santa Catarina")
  states <- factor(country_names[states],levels=order)
  return(states)
}

#' Subset epidemic states
#'
#' Gets those states that fulfill the epidemic criteria of 2sd above mean since July for 3 months
#' @param dat the data frame of data to look through (microcephaly or incidence)
#' @param micro if TRUE, looks at microcephaly data. Incidence data otherwise
#' @param lim in days, the start of the period of consideration
#' @return a list with a vector of epidemic and non-epidemic states
#' @export
get_epidemic_states <- function(dat,micro=TRUE,lim=575){
    if(micro){
        include <- NULL
        exclude <- NULL
        for(local in unique(dat$local)){
            tmpDat <- dat[dat$local==local,c("startDay","microCeph")]
            tmp <- tmpDat[tmpDat$startDay < lim,2]
            mean <- mean(tmp)
            sd <- sd(tmp)
            a <- tmpDat[tmpDat$startDay >= lim,2] >= (mean + 2*sd)
            if(any(cumul_true(a) >= 3)){
                include <- c(include, local)
            } else {
                exclude <- c(exclude, local)
            }
        }
    } else  {
        include <- NULL
        exclude <- NULL
        for(local in unique(dat$local)){
            tmpDat <- dat[dat$local==local,c("startDay","inc")]
            tmp <- tmpDat[tmpDat$startDay < lim,2]
            mean <- mean(tmp)
            sd <- sd(tmp)
            a <- tmpDat[tmpDat$startDay >= lim,2] >= (mean + 2*sd)
            if(any(cumul_true(a) >= 3)){
                include <- c(include, local)
            } else {
                exclude <- c(exclude, local)
            }
        }
    }
    return(list("include"=include,"exclude"=exclude))
}

#' Read all data
#'
#' Function to read in microcephaly and incidence data. Converts strings to factors, and orders them by microcephaly cases
#' @param microCephFile full filepath to microcephaly data
#' @param incDatFile full filepath ot incidence data
#' @param epidemicStates if TRUE, gets only epidemic states. Otherwise gives all states
#' @param order if TRUE, orders the states by per capita microcephaly rate
#' @param convert if TRUE, converts strings to factors
#' @return a list of microcephaly and incidence data
#' @export
read_all_data <- function(microCephFile="~/net/home/zika/data/microcephaly_data.csv",
                          incDatFile="~/net/home/zika/data/inc_data.csv",
                          epidemicStates=FALSE, order=TRUE, convert=TRUE){
    ## Read in microcephaly data
    allDat <- read.csv(microCephFile,stringsAsFactors=FALSE)

    ## Read in incidence data
    incDat <- read.csv(incDatFile,stringsAsFactors=FALSE)

    if(epidemicStates){
        states <- get_epidemic_states()
        allDat <- allDat[allDat$local %in% states$include,]
        incDat <- incDat[incDat$local %in% states$include,]
    }

    ## Convert to characters and factors
    if(convert){
        allDat$local <- as.character(allDat$local)
        allDat$local <- convert_name_to_state_factor(allDat$local)
        incDat$local <- as.character(incDat$local)
        incDat$local <- convert_name_to_state_factor(incDat$local)
    }
    if(order){
        order <- get_correct_order(named=convert,normalised=TRUE)
        allDat$local <- factor(allDat$local,levels=order)
        incDat$local <- factor(incDat$local,levels=order)
    }
    
    return(list(microDat=allDat,incDat=incDat))
}

#' Cumulative true
#'
#' Gets number of consecutive values in a vector that are TRUE
#' @param x the vector to be tested
#' @return a vector counting number of successive TRUEs
#' @export
cumul_true <- function(x)  {
  #x <- !x
  rl <- rle(x)
  len <- rl$lengths
  v <- rl$values
  cumLen <- cumsum(len)
  z <- x
  # replace the 0 at the end of each zero-block in z by the
  # negative of the length of the preceding 1-block....
  iDrops <- c(0, diff(v)) < 0
  z[ cumLen[ iDrops ] ] <- -len[ c(iDrops[-1],FALSE) ]
  # ... to ensure that the cumsum below does the right thing.
  # We zap the cumsum with x so only the cumsums for the 1-blocks survive:
  x*cumsum(z)
}

#' State order
#'
#' Gets the correct state order for all plots, sorted by number of microcephaly cases
#' @param birth_dat the data frame of data
#' @param named if TRUE, converts the state names to full names
#' @param normalised if TRUE, orders by per birth microcephaly rate
#' @return a vector of ordered state names
#' @export
get_correct_order <- function(named=TRUE,normalised=FALSE){
    birth_dat <- read.csv("~/net/home/zika/data/microcephaly_data.csv",stringsAsFactors=FALSE)
    if(named) birth_dat$local <- convert_name_to_state(birth_dat$local)
    if(normalised){
        order <- names(sort(sapply(unique(birth_dat$local),function(x) sum(birth_dat[birth_dat$local==x,"microCeph"]/birth_dat[birth_dat$local==x,"births"])),decreasing = TRUE))
    } else {
        order <- names(sort(sapply(unique(birth_dat$local),function(x) sum(birth_dat[birth_dat$local==x,"microCeph"])),decreasing = TRUE))
    }
    return(order)
}

## Converts the simplified state to its full name as a character
convert_name_to_state <- function(name){
  country_names <- c("colombia"="Colombia","pernambuco"="Pernambuco","amapa"="Amapá","amazonas" ="Amazonas",
                     "distritofederal" = "Distrito Federal","bahia"="Bahia","saopaulo"="São Paulo",
                     "paraiba"="Paraíba","maranhao"="Maranhão","ceara"="Ceará","sergipe"="Sergipe",
                     "riodejaneiro"="Rio de Janeiro","piaui"="Piauí","riograndedonorte"="Rio Grande Norte",
                     "minasgerais"="Minas Gerais", "matogrosso"="Mato Grosso","alagoas"="Alagoas",
                     "para"="Pará","acre"="Acre","espiritosanto"="Espírito Santo","goias"="Goiás",
                     "tocantins"="Tocantins","matogrossodosul"="Mato Grosso do Sul",
                     "matogrossdosul"="Mato Grosso do Sul","parana"="Paraná","riograndedosul"="Rio Grande do Sul",
                     "rondonia"="Rondônia","roraima"="Roraima","santacatarina"="Santa Catarina")
  return(country_names[name])
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
