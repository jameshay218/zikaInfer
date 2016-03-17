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
//' @param dat the matrix of data. 2 columns for counts of positive and total births, and N rows for each sampling time
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
  for(int i =0; i < alphas.nrow(); ++i){
    p = alphas(i,0)*R::pnorm(threshold,mus[0],sds[0],1,0) + alphas(i,1)*R::pnorm(threshold,mus[1],sds[1],1,0);
    //    Rcpp::Rcout << "Proportion: " << dat(i,0)/dat(i,1) << "   P: " << p << "   Likelihood: " << R::dbinom(dat(i,0),dat(i,1),p,1) << std::endl;
    lnlik += R::dbinom(dat(i,0),dat(i,1),p,1);
  }
  return(lnlik);
}

//' Calculate likelihood for probability data
//'
//' Given a matrix of threshold data (ie. microcephaly or not) for different sampling times, gives a log likelihood with the given parameter values
//' @param dat the matrix of data. 2 columns for counts of positive and total births, and N rows for each sampling time
//' @param alphas vector of probability of giving birth to an infant with microcephaly over time
//' @return a single log likelihood
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
double likelihood_prob(NumericMatrix dat, NumericVector alphas){
  double lnlik = 0;
  double p = 0;
  for(int i =0; i < alphas.length(); ++i){
    lnlik += R::dbinom(dat(i,0),dat(i,1),alphas[i],1);
  }
  return(lnlik);
}

//' Calculates probabilities of microcephaly over time for given time bucket sizes
//'
//' @param y the matrix of pregnant adult counts. First column should be times
//' @param probMicro probability of developing microcephaly given infection
//' @param baselineProb the baseline probability of giving birth to an infant with microcephaly
//' @param times matrix of times to create alphas over. First column is start of bucket, last column is end of bucket
//' @return the vector of alphas
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
NumericVector calculate_alphas_prob_buckets(NumericMatrix y, double probMicro, double baselineProb, NumericMatrix times){
  int i = 0;
  int index = 0;
  double start = 0;
  double end = 0;
  double tmp = 0;
  int j = 0;
  NumericVector alphas(times.nrow());

  // Get to start of y matrix that has relevant information
  while(y(i,0) < start) i++;

  // Go through all of the y matrix
  while(i < y.nrow() && index < alphas.size()){
    tmp = 0;
    j = 0;

    // Get upper and lower bounds of time
    start = times(index,0);
    end = times(index,1);
    //    Rcpp::Rcout << start << " " << end << std::endl;

    // Increase y index until at start of bucket
    while(y(i,0) < start) i++;

    // Go through y until end of bucket, storing incidence over this time
    while(y(i,0) < end && i < y.nrow()){
      tmp += y(i,1)/(y(i,1)+y(i,2)+y(i,3)+y(i,4));
      j++;
      i++;
    }

    // Take average incidence over this time to be alpha for this bucket
    alphas[index] = baselineProb + probMicro*tmp/j;
    if(alphas[index] > 1) alphas[index] = 1;
    // Increase index of bucket/alpha
    index++;
  }
  return(alphas);
}

//' Calculates probability of microcephaly over time with regular sampling intervals
//'
//' @param y the matrix of pregnant adult counts. First column should be times
//' @param probMicro probability of developing microcephaly given infection
//' @param baselineProb the baseline probability of giving birth to an infant with microcephaly
//' @param sampFreq the frequency of sampling
//' @return the value converted to a linear scale
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
NumericVector calculate_alphas_prob_sampfreq(NumericMatrix y, double probMicro, double baselineProb, int sampFreq){
  int i = 0;
  int index = 0;
  int size = y.nrow()/sampFreq;
  double tmp = 0;
  int j = 0;
  NumericVector alphas(size);
  while(i < y.nrow() && index < alphas.size()){
    tmp = 0;
    j = 0;
    while(j < sampFreq && i < y.nrow()){
      tmp += y(i,0)/(y(i,0)+y(i,1)+y(i,2)+y(i,3));
      j++;
      i++;
    }
    alphas[index] = baselineProb + probMicro*tmp/j;
    if(alphas[index] > 1) alphas[index] = 1;
    index++;
  }
  return(alphas);
}




//[[Rcpp::export]]
NumericVector p_test(NumericMatrix dat, NumericMatrix alphas, NumericVector mus, NumericVector sds, double threshold){
  double lnlik = 0;
  double p = 0;
  NumericVector ps(dat.nrow());
  for(int i =0; i < alphas.nrow(); ++i){
    p = alphas(i,0)*R::pnorm(threshold,mus[0],sds[0],1,0) + alphas(i,1)*R::pnorm(threshold,mus[1],sds[1],1,0);
    ps[i] = p;
  }
  return(ps);
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
NumericVector calculate_alphas(NumericMatrix y, double probMicro, int sampFreq){
  int i = 0;
  int index = 0;
  int size = y.nrow()/sampFreq;
  double tmp = 0;
  int j = 0;
  //Rcpp::Rcout << size << std::endl;
  NumericVector alphas(size);
  while(i < y.nrow() && index < alphas.size()){
    tmp = 0;
    j = 0;
    while(j < sampFreq && i < y.nrow()){
      tmp += y(i,0)/(y(i,0)+y(i,1)+y(i,2)+y(i,3));
      j++;
      i++;
    }
    alphas[index] = probMicro*tmp/j;
    index++;
    //    if(j > 0) alphas[index++] = probMicro*tmp/j;
  }
  // if(index < (alphas.size())) Rcpp::Rcout << "Error - invalid index" << std::endl;
  //Rcpp::Rcout << index << std::endl;
  return(alphas);
}



//' Calculates alphas for given time bucket sizes
//'
//' @param y the matrix of pregnant adult counts. First column should be times
//' @param probMicro probability of developing microcephaly given infection
//' @param times matrix of times to create alphas over. First column is start of bucket, last column is end of bucket
//' @return the vector of alphas
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
NumericVector calculate_alphas_buckets(NumericMatrix y, double probMicro, NumericMatrix times){
  int i = 0;
  int index = 0;
  double start = 0;
  double end = 0;
  double tmp = 0;
  int j = 0;
  NumericVector alphas(times.nrow());

  // Get to start of y matrix that has relevant information
  while(y(i,0) < start) i++;

  // Go through all of the y matrix
  while(i < y.nrow() && index < alphas.size()){
    tmp = 0;
    j = 0;

    // Get upper and lower bounds of time
    start = times(index,0);
    end = times(index,1);
    //    Rcpp::Rcout << start << " " << end << std::endl;

    // Increase y index until at start of bucket
    while(y(i,0) < start) i++;

    // Go through y until end of bucket, storing incidence over this time
    while(y(i,0) < end && i < y.nrow()){
      tmp += y(i,1)/(y(i,1)+y(i,2)+y(i,3)+y(i,4));
      j++;
      i++;
    }

    // Take average incidence over this time to be alpha for this bucket
    alphas[index] = probMicro*tmp/j;
    //    Rcpp::Rcout << alphas[index] << std::endl;
    // Increase index of bucket/alpha
    index++;
  }
  return(alphas);
}



//[[Rcpp::export]]
double scaletuning2(double step, double popt, double pcur){
  if(pcur>=1) pcur = 0.99;
  if(pcur<=0) pcur = 0.01;
  step = (step*R::qnorm(popt/2,0,1,1,0))/R::qnorm(pcur/2,0,1,1,0);
  return(step);
}


//[[Rcpp::export]]
double proposal_function(double current, double lower, double upper, double step){
  double update;
  double move;
  double new1;
  new1 = toUnitScale(current,lower,upper);
  
  do {
    new1 = toUnitScale(current,lower,upper);
    update = R::rnorm(0, 1);
    new1 = new1 + update*step;
  } while(new1 > 1 || new1 < 0);
  
  new1 = fromUnitScale(new1,lower,upper);
  
  return(new1);
}
