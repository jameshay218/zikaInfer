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

//' FOI calculation
//'
//' Calculates the force of infection over time for an SEIR model
//' @param IM a numeric vector of number of infected mosquitoes over time
//' @param NH the constant human population size
//' @param b the per vector per day bite rate
//' @param pMH the probability of transmission upon bite
//' @param the time step over which to bucket the force of infection. Best left to 1 for 1 day
//' @return the vector of FOI
//' @export
//' @useDynLib zikaProj 
//[[Rcpp::export]]
NumericVector generate_foi(NumericVector IM, double NH, double b, double pMH, double tstep){
  NumericVector foi(IM.size()/tstep);
  int max = foi.size();
  int index = 0;
  double tmp;
  for(int i = 0; i < max; ++i){
    tmp = 0;
    for(int j = 0; j < tstep; ++j){
      tmp += IM[index++]*b*pMH/NH;
    }
    foi[i] = tmp/tstep;
  }
  return(foi);
}

//' Cumulative susceptible probability
//' 
//' Calculates the cumulative probability of remaining susceptible over the course of the epidemic
//' @param foi the time varying force of infection as generated by \link{\code{generate_foi}}
//' @param tstep the time step for the buckets. Best left to 1.
//' @return the vector of cumulative escape probabilities
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
NumericVector generate_riskS(NumericVector foi, double tstep){
  NumericVector riskS(foi.size());
  int max = riskS.size();
  riskS[0] = exp(-foi[0]*tstep);
  for(int i = 1; i < max; ++i){
    riskS[i] = riskS[i-1]*exp(-foi[i]*tstep);
  }
  return(riskS);
}
//' Risk of infection
//' 
//' Calculates the relative risk of infection over time
//' @param foi the time varying force of infection as generated by \link{\code{generate_foi}}
//' @param riskS the cumulative probability of remaining susceptible over the course of the epidemic as generated by \link{\code{generate_riskS}}
//' @param tstep the time step for the buckets. Best left to 1.
//' @return the vector of relative infection risks
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
NumericVector generate_riskI(NumericVector foi, NumericVector riskS, double tstep){
  NumericVector riskI(foi.size());
  int max = riskI.size();
  riskI[0] = 1-exp(-foi[0]*tstep);
  for(int i = 1; i < max; ++i){
    riskI[i] = riskS[i-1]*(1-exp(-foi[i]*tstep));
  }
  return(riskI);
}


//' Time varying microcephaly risk AUX
//'
//' Calculates the probability of a birth having microcephaly at a given time
//' @param riskI the relative risk of infection over time as generated by \link{\code{generate_riskI}}
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation. Time steps should be the same as riskI
//' @param bp the baseline probability of microcephaly
//' @return the vector of microcephaly probabilities
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
NumericVector generate_probM_aux(NumericVector riskI, NumericVector probM, double bp){
  NumericVector allProbs(riskI.size());
  double tmp = 0;
  int max = riskI.size();
  int testDay = 0;
  int minDay = 0;
  int wow = 0;
  int week = 0;
  for(int i = 0; i < max; ++i){
    tmp = 0;
    testDay = i;
    minDay = testDay-probM.size();
    if(minDay < 0) minDay = 0;
    week = 0;
    for(int j = minDay; j < testDay; ++j){
      tmp += riskI[j]*probM[j-(testDay-probM.size())];
    }
    allProbs[i] = 1 - (1-tmp)*(1-bp);
  }
  return(allProbs);
}

//' Time varying microcephaly risk
//'
//' Given the epidemic model parameters and number of infected mosquitoes over time, calculates the probability of a birth having microcephaly at a given time
//' @param IM a numeric vector of number of infected mosquitoes over time
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation
//' @param NH the constant human population size
//' @param b the per vector per day bite rate
//' @param pMH the probability of transmission upon bite
//' @param bp the baseline probability of microcephaly
//' @param tstep the time step for the buckets. Best left to 1.
//' @return the vector of microcephaly probabilities
//' @export
//' @useDynLib zikaProj
//' @seealso \link{\code{generate_probM_aux}}
//[[Rcpp::export]]
NumericVector generate_probM(NumericVector IM, NumericVector probM, double NH, double b, double pMH, double bp, double tstep){
  NumericVector foi = generate_foi(IM, NH, b, pMH, tstep);
  NumericVector riskS = generate_riskS(foi, tstep);
  NumericVector riskI = generate_riskI(foi, riskS, tstep);
  return(generate_probM_aux(riskI, probM, bp));
}

//' Likelihood function for time-varying microcephaly
//'
//' Calculates the likelihood of observing a vector of microcephaly births given the total number of births and microcephaly probabilities. Note that all vectors must be equal lengths.
//' @param microBirths the vector of observed microcephaly cases over time
//' @param allBirths the corresponding total number of births
//' @param probM the corresponding vector of microcephaly probabilities as calculated by \link{\code{generate_probM}}.
//' @return a single likelihood value
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
double likelihood_probM(NumericVector microBirths, NumericVector allBirths, NumericVector probM){
  double lnlik = 0;
  int max = probM.length();
  for(int i = 0; i < max; ++i){
    lnlik += R::dbinom(microBirths[i],allBirths[i],probM[i],1);
  }
  return(lnlik);
}

//' Likelihood function for time-varying microcephaly
//'
//' Calculates the likelihood of observing a vector of microcephaly births given the total number of births and microcephaly probabilities. Takes model generated/used parameters.
//' @param microBirths the vector of observed microcephaly cases over time
//' @param allBirths the corresponding total number of births
//' @param IM a numeric vector of number of infected mosquitoes over time
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation
//' @param NH the constant human population size
//' @param b the per vector per day bite rate
//' @param pMH the probability of transmission upon bite
//' @param bp the baseline probability of microcephaly
//' @param tstep the time step for the buckets. Best left to 1.
//' @return a single likelihood value
//' @export
//' @useDynLib zikaProj
//' @seealso \link{\code{likelihood_probM}}
//[[Rcpp::export]]
double likelihood_probM_all(NumericVector microBirths, NumericVector allBirths, NumericVector IM, NumericVector probM, double NH, double b, double pHM, double bp, double tstep){
  NumericVector allP = generate_probM(IM, probM, NH, b, pHM, bp, tstep);
  return(likelihood_probM(microBirths, allBirths, allP));
}

//' Averages a vector based on bucket sizes
//'
//' Given a vector (a) and another vector of bucket sizes, returns the averaged vector (a)
//' @param a the vector to be bucketed and averaged
//' @param buckets the vector of bucket sizes to average a over
//' @return the vector of averaged a
//' @export
//' @useDynLib zikaProj
//[[Rcpp::export]]
NumericVector average_buckets(NumericVector a, NumericVector buckets){
  NumericVector results(buckets.size());
  int index = 0;
  for(int i = 0; i < buckets.size(); ++i){
    results[i] = 0;
    for(int j = 0; j < buckets[i]; ++j){
      results[i] += a[index++];
    }
    results[i] = results[i]/(double)buckets[i];
  }
  return(results);
}
