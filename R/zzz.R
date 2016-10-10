CsimpleSEIR_rich <- NULL
.onLoad <- function(...) {
    CsimpleSEIR_rich <<- getNativeSymbolInfo("simpleSEIR_rich", "zikaProj")
}
