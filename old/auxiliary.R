#' Create count table
#'
#' Given a matrix of data points (rows = times, columns = individuals), returns the summarised count data
#' @param dat the matrix to be counted
#' @return a data frame of counted data
#' @export
#' @useDynLib zikaProj
createCounts <- function(dat){
    all <- NULL
    means <- NULL
    all <- apply(dat,1,table)
    all <- lapply(all,as.data.frame)
    all <- lapply(seq_along(all),function(x) cbind(all[[x]],x))
    all <- do.call("rbind",all)
    means <- apply(dat,1,mean)
    
    all[,1] <- as.numeric(as.character(all[,1]))
    colnames(all) <- c("Size","Proportion","Day")
    meanDat <- data.frame(x=seq(1,max(all$Day),by=1),y=means)
    return(list(all,meanDat))
}

#' Sets up R environment
#' @export
#' @useDynLib zikaProj
setupEnv <- function(){
    system("R CMD SHLIB mymod.c")
    dyn.load("mymod.so")
}


#' Un cumulative sum
#'
#' Given a vector representing a cumulative sum, returns a vector of the same length with new "incidence" for each time point
#' @param x the vector to be "uncumsummed"
#' @return the un cumulative summed vector
#' @export
#' @useDynLib zikaProj
un_cumsum <- function(x){
  tmp <- numeric(length(x)-1)
  i <- length(x)
  while(i > 1){
   tmp[i-1] <- x[i] - x[i-1]
   i <- i - 1
   }
  return(tmp)
}

#' Final size calculation
#'
#' Cost function to be used for the final epidemic size calculation. Using a given R0, if the argument of par[1] is the correct attack rate then this function will return 0. Try: nleqslv(A0, simeq, R0=3)
#' @param par array of length one containing the attack rate
#' @param R0 the R0 used for solving
#' @return difference between the LHS and RHS of the final size equation
#' @export
#' @useDynLib zikaProj
simeq <- function(par,R0){
    A <- par[1]
    f1 <- A - (1-exp(-R0*A))
    f1
}

calc.alphas <- function(y, sampFreq,probMicro){
    daysPerYear <- nrow(y)/max(y$times)
    i <- 1 + sampFreq
    index <- 1
    all <- NULL
    birthsPerYear <- sum(y0s[4:6])/sampFreq
    birthsPerDay <- ceiling(birthsPerYear/daysPerYear)
    alphas <- numeric(ceiling(max(y$times)*daysPerYear/sampFreq))
    while(i <= nrow(y)){
        is <- y[(i-sampFreq):i, "If"]
        ns <- rowSums(y[(i-sampFreq):i,c("Sf","Ef","If","Rf")])
        propns <- mean(na.omit(is/ns))
        alphas[index] <- probMicro * propns
        index <- index + 1
        i <- i + sampFreq
    }
    return(alphas)
}
