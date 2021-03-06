% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posteriors.R
\name{create_posterior}
\alias{create_posterior}
\title{Simplify posterior call}
\usage{
create_posterior(parTab, data, PRIOR_FUNC, ...)
}
\arguments{
\item{parTab}{the parameter table. See \code{\link{exampleParTab}}}

\item{data}{the microcephaly data. See \code{\link{exampleMicroDat}}}

\item{PRIOR_FUNC}{pointer to a function to calculate prior values. Can take "pars" and "..." as arguments, where "pars" is the model parameters.}

\item{incDat}{should also include incidence data, though this can be left as NULL. See \code{\link{exampleIncDat}}.}

\item{peakTimes}{include data on epidemic peak time ranges. Can be left as NULL. See \code{\link{examplePeakTimes}}}

\item{ts}{time vector over which to solve the ODE model}

\item{version}{usually just leave this as "binomial", indicating assumed binomially distributed errors. I've added a "forecast" version which expects parameters to do with the lack of second wave.}
}
\value{
a single value for the posterior
}
\description{
Simplifies the call to the posterior function using closures. This function is set up to work with the lazymcmc package, found at github.com/jameshay218/lazymcmc
}

