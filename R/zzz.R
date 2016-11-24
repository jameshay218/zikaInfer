CsimpleSEIR_rich <- NULL
#' @useDynLib zikaProj
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @import grDevices
#' @import graphics
#' @import ggplot2
#' @import gridExtra
#' @import gtable
.onLoad <- function(...) {
    CsimpleSEIR_rich <<- getNativeSymbolInfo("simpleSEIR_rich", "zikaProj")
}
