#include <Rcpp.h>
using namespace Rcpp;

//' Converts to unit scale
//'
//' @param x the double to be converted
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a unit scale
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
double toUnitScale(double x, double min, double max){
  return((x-min)/(max-min));
}

//' Converts to linear scale
//'
//' @param x the double to be converted back to linear scale
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a linear scale
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}

//' Converts to linear scale
//'
//' @param x the double to be converted back to linear scale
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a linear scale
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
double likelihood(NumericMatrix dat, NumericMatrix alphas, NumericVector mus, NumericVector sds){
  double lnlik = 0;
  double tmp = 0;
  for(int i = 0; i < dat.nrow(); ++i){
    for(int j = 0; j < dat.ncol(); ++j){
        if(!NumericVector::is_na(dat(i,j))){
	tmp = 0;
	for(int q = 0; q < alphas.ncol();++q){
	  tmp += alphas(i,q)*(pow(2.0*M_PI*pow(sds(q),2),-0.5) * exp(-0.5*pow(dat(i,j)-mus(q),2) / pow(sds(q),2)));
	}
	lnlik += log(tmp);
	 }
    }
  }
  return(lnlik);
}

//' Converts to linear scale
//'
//' @param x the double to be converted back to linear scale
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a linear scale
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
NumericVector calculate_alphas(NumericMatrix y, double probMicro, double sampFreq){
  int i = 0 + sampFreq;
  int index = 0;
  NumericVector alphas(y.nrow()/sampFreq);
  NumericVector tmpPropn(sampFreq);
  Rcpp::Rcout << "Here" << std::endl;
  while(i < y.nrow()){
    Rcpp::Rcout << i << std::endl;
    for(int j =0;j <= sampFreq;++j){
      tmpPropn[j] = y(i-j,0)/(y(i-j,0)+y(i-j,1)+y(i-j,2)+y(i-j,3));
    }
    alphas[index++] = probMicro*Rcpp::mean(tmpPropn);
    i += sampFreq;
  }
  return(alphas);
}

//[[Rcpp::export]]
double proposal_function(double current, double step){
  double update;
  double move;
  update = R::rnorm(current, step*step);
  return(update);
  
}

//[[Rcpp::export]]
double scaletuning2(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  return(step);
}

