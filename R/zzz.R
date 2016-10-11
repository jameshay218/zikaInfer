CsimpleSEIR_rich <- NULL
#' @useDynLib zikaProj
#' @importFrom Rcpp evalCpp
#' @importFrom deSolve ode
.onLoad <- function(...) {
    CsimpleSEIR_rich <<- getNativeSymbolInfo("simpleSEIR_rich", "zikaProj")
}

#' @import gridExtra
#' @import gtable
#' @import data.table
#' @import plyr
#' @import lattice
#' @import reldist
#' @import coda
NULL
