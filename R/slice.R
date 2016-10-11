#' Slice MCMC sampler
#' 
#' Typically uses a univariate slice sampler for proposal update. Performs MCMC on a given target function
#' @param target pointer to the posterior (target) function, that should take only one argument - the vector of model parameters
#' @param x_init initial paramter values for the target function
#' @param fixed a vector of 1/0 (TRUE/FALSE) matching the target function, indicating which parameters should be updated
#' @param nsteps number of iterations
#' @param w vector of widths for the slice sampler
#' @param lower optional vector of lower parameter bounds
#' @param upper options vector of upper parameter points
#' @param print_every how often should output be printed
#' @return a dataframe of the MCMC chain
#' @export
#' @useDynLib
mcmc_simple <- function(target, x_init, fixed, nsteps, w, lower=-Inf, upper=-Inf, print_every=1,save_file=NULL){
    npar <- length(x_init)
    hist_pars <- matrix(NA, ncol=npar, nrow=nsteps)
    hist_prob <- rep(NA, nsteps)
    colnames(hist_pars) <- names(x_init)

    lower <- recycle(lower, npar)
    upper <- recycle(upper, npar)
    w     <- recycle(w,     npar)

    check_bounds(lower, upper, x_init)
    
    y_init <- target(x_init)

    if (!is.finite(y_init)) {
        stop("Starting point must have finite probability")
    }

    we_should_print <- make_every_so_often(print_every)
    if(!is.null(save_file)) write.table(t(c(i=0,x_init,lnlik=y_init)),file=save_file,row.names=FALSE,col.names=FALSE,sep=",",append=FALSE)
    
    last_save <- 1
    for (i in seq_len(nsteps)) {
        tmp <- sampler_slice(target, x_init, fixed, y_init, w, lower, upper)
        x_init <- hist_pars[i,] <- tmp[[1]]
        y_init <- hist_prob[i]  <- tmp[[2]]
        if (we_should_print()) {
            message(sprintf("%d: {%s} -> %2.5f", i,
                            paste(sprintf("%2.4f", x_init), collapse=", "),
                            y_init))
            if(!is.null(save_file)){
                write.table(data.frame(i=seq(last_save,i,by=1),hist_pars[last_save:i,],lnlik=hist_prob[last_save:i]),file=save_file,row.names=FALSE,col.names=FALSE,sep=",",append=TRUE)
                last_save <- i+1
            }
        }
    }

    data.frame(i=seq_along(hist_prob), hist_pars, hist_prob)
}

#' Slice sampler
#'
#' Multivariate slice sampler for given target function
#' @param lik target function
#' @param x_init intial parameter values for the target function
#' @param fixed vector of bools to inidicate which parameters should be updated
#' @param y_init current target function value
#' @param w vector of widths for the slice sampler
#' @param lower optional vector of lower parameter bounds
#' @param upper options vector of upper parameter points
#' @return a list. First component is updated parameters, second is corresponding target function evaluation
#' @export
#' @useDynLib
sampler_slice <- function(lik, x_init, fixed, y_init, w, lower, upper) {
  for (i in fixed) {
    xy <- slice_1d(make_unipar(lik, x_init, i),
                   x_init[i], y_init, w[i], lower[i], upper[i])
    x_init[i] <- xy[1]
    y_init    <- xy[2]
  }

  list(x_init, y_init)
}

#' Univariate slice sampler
#'
#' Univariate slice sampler for given target function
#' @param f target function
#' @param x_init intial parameter value for the target function
#' @param y_init current target function value
#' @param w width for the slice sampler
#' @param lower lower parameter bounds
#' @param upper upper parameter points
#' @return a single sampled value
slice_1d <- function(f, x_init, y_init, w, lower, upper) {
  z <- y_init - rexp(1)
  r <- slice_isolate(f, x_init, y_init, z, w, lower, upper)
  slice_sample(f, x_init, z, r)
}

#' Slice sampler
#'
#' Finds actual bounds for the slice sampler
slice_isolate <- function(f, x_init, y_init, z, w, lower, upper) {
  u <- runif(1) * w
  L <- x_init - u
  R <- x_init + (w-u)

  while (L > lower && f(L) > z) {
    L <- L - w
  }
  while (R < upper && f(R) > z) {
    R <- R + w
  }

  c(max(L, lower), min(R, upper))
}

#' Working function for slice sampler
slice_sample <- function(f, x_init, z, r) {
  r0 <- r[1]
  r1 <- r[2]

  repeat {
    xs <- runif(1, r0, r1)
    ys <- f(xs)
    if (ys > z) {
      break
    }
    if (xs < x_init) {
      r0 <- xs
    } else {
      r1 <- xs
    }
  }
  c(xs, ys)
}

## utilities:
make_unipar <- function(f, x, i) {
  force(f)
  force(x)
  force(i)
  function(z) {
    x[i] <- z
    f(x)
  }
}

#' Checks/forces given vector to desired length
recycle <- function(x, length, name=deparse(substitute(x))) {
  if (length(x) == 1) {
    rep(x, length)
  } else if (length(x) == length) {
    x
  } else {
    stop(sprintf("'%s' of incorrect length", name))
  }
}

#' Checks that given value is within bounds
check_bounds <- function(lower, upper, x0=NULL) {
  if (!is.null(x0) && (any(x0 < lower) || any(x0 > upper))) {
    stop("Starting parameter falls outside of problems bounds")
  }
  if (any(lower >= upper)) {
    stop("'upper' must be strictly greater than 'lower'")
  }
}

#' Function that evaluations to true every nth time it is called
make_every_so_often <- function(iterations=1) {
  if (iterations == 1L) {
    function() {
      TRUE
    }
  } else if (iterations > 0) {
    i <- 0L # counter, will be updated
    iterations <- as.integer(iterations)
    function() {
      i <<- i + 1L
      i %% iterations == 0L
    }
  } else {
    function() {
      FALSE
    }
  }
}

#' Basic MCMC function
mcmc <- function(target, x_init, nsteps, w,
                 lower=-Inf, upper=Inf, print_every=1) {
  npar <- length(x_init)
  hist_pars <- matrix(NA, ncol=npar, nrow=nsteps)
  hist_prob <- rep(NA, nsteps)

  colnames(hist_pars) <- names(x_init)

  lower <- recycle(lower, npar)
  upper <- recycle(upper, npar)
  w     <- recycle(w,     npar)

  check_bounds(lower, upper, x_init)

  y_init <- target(x_init)
  if (!is.finite(y_init)) {
    stop("Starting point must have finite probability")
  }

  we_should_print <- make_every_so_often(print_every)

  for (i in seq_len(nsteps)) {
    tmp <- sampler_slice(target, x_init, y_init, w, lower, upper)
    x_init <- hist_pars[i,] <- tmp[[1]]
    y_init <- hist_prob[i]  <- tmp[[2]]

    if (we_should_print()) {
      message(sprintf("%d: {%s} -> %2.5f", i,
                      paste(sprintf("%2.4f", x_init), collapse=", "),
                      y_init))
    }
  }

  data.frame(i=seq_along(hist_prob), hist_pars, hist_prob)
}

#' Simplify posterior function
#'
#' Simplifies the Zika model posterior function to only require the parameter values
#' @param t_pars time vector over which to solve the ODE model
#' @param values ODE parameters
#' @param names names for the ODE parameters
#' @param local names of all the states that are being considered
#' @param startDays vector of the start day for each bucket of data
#' @param endDays vector of the end day for each bucket of data
#' @param buckets bucket sizes
#' @param microCeph vector of microcephaly incidence for each bucket
#' @param births vector of total briths for each bucket
#' @param data_locals corresponding vector of state names for the microcephaly data
#' @param incDat defaults to NULL. If provided, uses the incidence data to improve parameter estimates.
#' @param actualPeakTimes defaults to NULL. Data frame of priors belief on incidence peak times. Provide a row for each state and three columns - the name of the state, the lower bound and the upper bound on peak time. Puts a uniform prior on this range.
#' @param allPriors defaults to NULL. Arguments for parameter priors, if desired.
#' @return a single value for the posterior
#' @export
#' @useDynLib zikaProj
make_posterior <- function(t_pars, values, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat=NULL,actualPeakTimes=NULL,allPriors=NULL, ...){
    counter <- 1L
    function(x){
        counter <<- counter + 1L
        return(posterior_complex_buckets(t_pars, x, names, local, startDays, endDays, buckets, microCeph, births, data_locals, incDat, actualPeakTimes,allPriors))
    }
}

