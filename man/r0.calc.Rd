% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model.R
\name{r0.calc}
\alias{r0.calc}
\title{R0 calculation}
\usage{
r0.calc(params)
}
\arguments{
\item{params}{Vector of parameters matching those returned by \code{\link{setupListPars}}}
}
\value{
A single value for R0
}
\description{
Calculates the R0 of the SEIR model given a vector of parameters. R0 defined as number of expected human cases given introduction of 1 infected human into a totally naive population of humans and mosquitoes.
}
\seealso{
\code{\link{b.calc}}
}

