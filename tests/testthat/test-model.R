context("model")

test_that("integrators agree", {
  ts <- seq(0,3003,by=1)
  pars <- setupParsLong()
  pars <- c(pars, "L_H"=70*365, "N_H"=1000000,"density"=3,"constSeed"=100)
  y0s <- generate_y0s(pars["N_H"], pars["density"])

  lsoda_y <- solveModelSimple_lsoda(ts, y0s, pars, TRUE)
  class(lsoda_y) <- "matrix"
  at <- attributes(lsoda_y)
  attributes(lsoda_y) <- at[names(at) %in% c("class","dim","dimnames")]
  
  rlsoda_y <- solveModelSimple_rlsoda(ts,y0s,pars,TRUE)
  expect_equal(rlsoda_y,lsoda_y, tolerance=1e-5)
})



speed_test <- function(){
    library(inline)
    cpp_if_src <- '
  Rcpp::NumericVector xa(a);
  int n_xa = xa.size();
  for(int i=0; i < n_xa; i++) {
    if(xa[i]<0) xa[i] = 0;
  }
  return xa;
'
    cpp_if <- cxxfunction(signature(a="numeric"), cpp_if_src, plugin="Rcpp")
    
    ts <- seq(0,3003,by=1)
    pars <- setupParsLong()
    pars <- c(pars, "L_H"=70*365, "N_H"=1000000,"density"=3,"constSeed"=100)
    y0s <- generate_y0s(pars["N_H"], pars["density"])

    a <- function(){
        lsoda_y <- solveModelSimple_lsoda(ts, y0s, pars, TRUE)
    }
    b <- function(){
        rlsoda_y <- solveModelSimple_rlsoda(ts,y0s,pars,TRUE)
        rlsoda_y[,"I_H"] <- (abs(rlsoda_y[,"I_H"])  + rlsoda_y[,"I_H"])/2
    }

    c <- function(){
        rlsoda_y <- solveModelSimple_rlsoda(ts,y0s,pars,TRUE)
        rlsoda_y[,"I_H"] <- cpp_if(rlsoda_y[,"I_H"])
    }
    d <- function(){
        rlsoda_y <- solveModelSimple_rlsoda(ts,y0s,pars,TRUE)
        rlsoda_y[rlsoda_y[,"I_H"] < 0,"I_H"] <- 0
    }
    res <- microbenchmark::microbenchmark(a(),b(),c(),d())

}
