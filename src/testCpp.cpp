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

//' Calculate gaussian mixture model likelihood
//'
//' Given a matrix of data points (columns represent individuals, rows represent sampling times), computes the log likelihood using the gaussian mixture model with known alphas, mus and sds.
//' @param dat the matrix of data
//' @param alphas matrix of of alphas. 2 columns (1 for both distributions) with number of rows matching number of rows in dat
//' @param mus vector of mus for the 2 distributions 
//' @param sds vector of sds for the 2 distributions
//' @return a single log likelihood
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


//' Calculate likelihood for threshold data
//'
//' Given a matrix of threshold data (ie. microcephaly or not) for different sampling times, gives a log likelihood with the given parameter values
//' @param dat the matrix of data. 2 columns for counts of positive and negative, and N rows for each sampling time
//' @param alphas matrix of of alphas. 2 columns (1 for both distributions) with number of rows matching number of rows in dat
//' @param mus vector of mus for the 2 distributions 
//' @param sds vector of sds for the 2 distributions
//' @param threshold the threshold value below which a data point is classified as having microcephaly
//' @return a single log likelihood
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
double likelihood_threshold(NumericMatrix dat, NumericMatrix alphas, NumericVector mus, NumericVector sds, double threshold){
  double lnlik = 0;
  double p = 0;
  for(int i =0; i <dat.nrow(); ++i){
    p = alphas(i,0)*R::pnorm(threshold,mus[0],sds[0],1,0) + alphas(i,1)*R::pnorm(threshold,mus[1],sds[1],1,0);
    lnlik += R::dbinom(dat(i,0),dat(i,1),p,1);
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
  while(i < y.nrow()){
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

