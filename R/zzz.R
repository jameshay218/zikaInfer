C_SEIR_model_rlsoda <- NULL
C_SEIR_model_lsoda <- NULL
C_initmodSEIR <- NULL
#' @useDynLib zikaProj, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import stats
#' @import utils
#' @import grDevices
#' @import graphics
#' @import ggplot2
#' @import gridExtra
#' @import gtable
.onLoad <- function(...) {
  C_SEIR_model_rlsoda<<- getNativeSymbolInfo("SEIR_model_rlsoda", PACKAGE = "zikaProj")
  C_SEIR_model_lsoda <<- getNativeSymbolInfo("SEIR_model_lsoda", PACKAGE = "zikaProj")
  C_initmodSEIR <<- getNativeSymbolInfo("initmodSEIR", PACKAGE = "zikaProj")
  
}
