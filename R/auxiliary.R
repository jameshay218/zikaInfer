#' Create count table
#'
#' Given a matrix of data points (rows = times, columns = individuals), returns the summarised count data
#' @param dat the matrix to be counted
#' @return a data frame of counted data
#' @export
#' @useDynLib zikaProj
createCounts <- function(dat){
    all <- NULL
    for(i in 1:nrow(dat)){
        tmp <- as.data.frame(table(dat[i,]))
        tmp[,2] <- tmp[,2]/sum(tmp[,2])
        tmp <- cbind(tmp, i)
        all <- rbind(all, tmp)
    }
    return(all)
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
