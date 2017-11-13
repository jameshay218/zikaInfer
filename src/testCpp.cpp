#include <Rcpp.h>
using namespace Rcpp;

//' Converts to unit scale
//'
//' @param x the double to be converted
//' @param min the minimum value on the linear scale
//' @param max the maximum value on the linear scale
//' @return the value converted to a unit scale
//' @export
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
//[[Rcpp::export]]
double fromUnitScale(double x, double min, double max){
  return(min + (max-min)*x);
}

//' FOI calculation
//'
//' Calculates the force of infection over time for an SEIR model
//' @param IM a numeric vector of number of infected mosquitoes over time
//' @param NH the constant human population size
//' @param b the per vector per day bite rate
//' @param pMH the probability of transmission upon bite
//' @param tstep the time step over which to bucket the force of infection. Best left to 1 for 1 day
//' @return the vector of FOI
//' @export
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
//' @param foi the time varying force of infection as generated by generate_foi
//' @param tstep the time step for the buckets. Best left to 1.
//' @return the vector of cumulative escape probabilities
//' @export
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
//' @param foi the time varying force of infection as generated by generate_foi
//' @param riskS the cumulative probability of remaining susceptible over the course of the epidemic as generated by generate_riskS
//' @param tstep the time step for the buckets. Best left to 1.
//' @return the vector of relative infection risks
//' @export
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
//' @param riskI the relative risk of infection over time as generated by generate_riskI
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation. Time steps should be the same as riskI
//' @param bp the baseline probability of microcephaly
//' @return the vector of microcephaly probabilities
//' @export
//[[Rcpp::export]]
NumericVector generate_probM_aux(NumericVector riskI, NumericVector probM, double bp){
  NumericVector allProbs(riskI.size());
  double tmp = 0;
  int max = riskI.size();
  int testDay = 0;
  int minDay = 0;

  //double * c_riskI = REAL(riskI), c_probM = REAL(probM), c_allProbs = REAL(allProbs);
  int probM_size = probM.size();

  for(int i = 0; i < max; ++i){
    tmp = 0;
    testDay = i;
    minDay = testDay - probM_size;
    if(minDay < 0) minDay = 0;
    for(int j = minDay; j < testDay; ++j){
      tmp += riskI[j]*probM[j-(testDay-probM_size)];
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
//[[Rcpp::export]]
NumericVector generate_probM(NumericVector IM, double NH, NumericVector probM, double b, double pMH, double bp, double tstep){
  NumericVector foi = generate_foi(IM, NH, b, pMH, tstep);
  NumericVector riskS = generate_riskS(foi, tstep);
  NumericVector riskI = generate_riskI(foi, riskS, tstep);
  return(generate_probM_aux(riskI, probM, bp));
}

//' Likelihood function for time-varying microcephaly
//'
//' Calculates the likelihood of observing a vector of microcephaly births given the total number of births and microcephaly probabilities. Note that all vectors must be equal lengths.. Assuming binomial distribution.
//' @param microBirths the vector of observed microcephaly cases over time
//' @param allBirths the corresponding total number of births
//' @param probM the corresponding vector of microcephaly probabilities as calculated by generate_probM.
//' @return a single likelihood value
//' @export
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
//' Calculates the likelihood of observing a vector of microcephaly births given the total number of births and microcephaly probabilities. Note that all vectors must be equal lengths. Assuming normal distribution.
//' @param microBirths the vector of observed microcephaly cases over time
//' @param allBirths the corresponding total number of births
//' @param probM the corresponding vector of microcephaly probabilities as calculated by generate_probM.
//' @param lik_sd assumed standard deviation of the likelihood function
//' @return a single likelihood value
//' @export
//[[Rcpp::export]]
double likelihood_probM_norm(NumericVector microBirths, NumericVector allBirths, NumericVector probM, double lik_sd){
  double lnlik = 0;
  int max = probM.length();
  for(int i = 0; i < max; ++i){
    //Rcpp::Rcout << microBirths[i]/allBirths[i] << " " << probM[i] << std::endl;
    lnlik += R::dnorm(microBirths[i],probM[i]*allBirths[i],lik_sd,1);
  }
  return(lnlik);
}


//' Averages a vector based on bucket sizes
//'
//' Given a vector (a) and another vector of bucket sizes, returns the averaged vector (a)
//' @param a the vector to be bucketed and averaged
//' @param buckets the vector of bucket sizes to average a over
//' @return the vector of averaged a
//' @export
//[[Rcpp::export]]
NumericVector average_buckets(NumericVector a, NumericVector buckets){
  NumericVector results(buckets.size());
  int index = 0;
  int j = 0;
  for(int i = 0; i < buckets.size(); ++i){
    results[i] = 0;
    j = 0;
    for(j = 0; (j < buckets[i]) & (index < a.size()); ++j){
      results[i] += a[index++];
    }
    results[i] = results[i]/j;//(double)buckets[i];
  }
  return(results);
}

//' Sums a vector based on bucket sizes
//'
//' Given a vector (a) and another vector of bucket sizes, returns the summed vector (a)
//' @param a the vector to be bucketed
//' @param buckets the vector of bucket sizes to sum a over
//' @return the vector of summed a
//' @export
//[[Rcpp::export]]
NumericVector sum_buckets(NumericVector a, NumericVector buckets){
  NumericVector results(buckets.size());
  int index = 0;
  for(int i = 0; i < buckets.size(); ++i){
    results[i] = 0;
    for(int j = 0; (j < buckets[i]) & (index < a.size()); ++j){
      results[i] += a[index++];
    }
  }
  return(results);
}



//' Time varying microcephaly risk not aborted
//'
//' Calculates the probability of a birth having microcephaly at a given time
//' @param riskI the relative risk of infection over time as generated by generate_riskI
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation. Time steps should be the same as riskI
//' @param bp the baseline probability of microcephaly
//' @param abortion_rate probability of a women having an abortion given that she was infected
//' @param birth_reduction reduced probability of getting infected, synonymous with probability of avoiding pregnancy
//' @param switch_t integer for index after which behaviour changes
//' @param abortion_time integer for last allowable abortion day/week in gestational time
//' @param ABORTED if TRUE, then return the number of aborted cases. If FALSE, return the number of not aborted cases.
//' @return the vector of microcephaly probabilities
//' @export
//[[Rcpp::export]]
NumericVector generate_probM_forecast(NumericVector riskI, NumericVector probM, double bp, double abortion_rate, double birth_reduction, int switch_t, int abortion_time, bool ABORTED){
  NumericVector allProbs(riskI.size());
  double tmp = 0;
  double tmp1 = 0;
  double tmp_bp = bp;
  double day_prob  = 0;
  double aborted_prop = 1 - abortion_rate;
  double aborted_prop_before = 0;
  int max = riskI.size();
  int testDay = 0;
  int minDay = 0;

  if(ABORTED){
    aborted_prop = abortion_rate;
    aborted_prop_before = 1;
  }
  int probM_size = probM.size();

  // For each day
  for(int i = 0; i < max; ++i){
    tmp = 0;
    day_prob = 0;
    testDay = i;
    minDay = testDay - probM_size;
    if(minDay < 0) minDay = 0;
    
    // Check all preceding 40 days.
    for(int j = minDay; j < testDay; ++j){
      // Prob that didn't get baseline microceph up to now
      tmp_bp = pow(1-bp,j-minDay);
      /* If after switch time and before abortion date cutoff, 
	 prob is probability of getting infection (modified 
	 for additional avoidance), multiplied by the probability
	 of getting microcephaly given infection in that gestational week
	 multiplied by the probability of not getting an abortion */
      if((j >= switch_t) & (j - (testDay - probM_size)) < abortion_time){
	tmp1 = riskI[j]*(1-birth_reduction)*probM[j-(testDay-probM_size)];
	// Prob of getting ZIKV associated microceph and didn't get baseline up to now
	// Or got baseline microceph this day
	// Or both happened on this day
	day_prob = tmp1*tmp_bp + tmp_bp*bp + tmp1*(1-tmp_bp)*bp*tmp_bp;
	day_prob = aborted_prop*day_prob;
	
	/* If after switch time but not within abortion time, just modified infection risk */
      } else if(j >= switch_t){
	
	tmp1 = riskI[j]*(1-birth_reduction)*probM[j-(testDay-probM_size)];
	day_prob = tmp1*tmp_bp + tmp_bp*bp + tmp1*(1-tmp_bp)*bp*tmp_bp;
	/* Otherwise, just normal risk of getting infected multiplied by risk of
	   microcephaly given infection in that gestational week */
      } else {
	tmp1 = riskI[j]*probM[j-(testDay-probM_size)];
	day_prob = tmp1*tmp_bp + tmp_bp*bp + tmp1*(1-tmp_bp)*bp*tmp_bp;
	day_prob = day_prob*(1-aborted_prop_before);
      }
      tmp += day_prob;
    }
    allProbs[i] = tmp;
  }
  return(allProbs);
}



//' Time varying microcephaly risk AUX
//'
//' Calculates the probability of a birth having microcephaly at a given time
//' @param riskI the relative risk of infection over time as generated by generate_riskI
//' @param probM the time varying risk of developing microcephaly given infection over the course of gestation. Time steps should be the same as riskI
//' @param bp the baseline probability of microcephaly
//' @param abortion_rate probability of a women getting an abortion given that she was infected
//' @param birth_reduction reduced probability of getting infected, synonymous with probability of avoiding pregnancy
//' @param switch_t integer for index after which behaviour changes
//' @param abortion_time integer for last allowable abortion day/week
//' @return the vector of microcephaly probabilities
//' @export
//[[Rcpp::export]]
NumericVector generate_probM_forecast_OLD(NumericVector riskI, NumericVector probM, double bp, double abortion_rate, double birth_reduction, int switch_t, int abortion_time){
  NumericVector allProbs(riskI.size());
  double tmp = 0;
  double day_prob  = 0;
  int max = riskI.size();
  int testDay = 0;
  int minDay = 0;

  //double * c_riskI = REAL(riskI), c_probM = REAL(probM), c_allProbs = REAL(allProbs);
  int probM_size = probM.size();

  for(int i = 0; i < max; ++i){
    tmp = 0;
    day_prob = 0;
    testDay = i;
    minDay = testDay - probM_size;
    if(minDay < 0) minDay = 0;
    for(int j = minDay; j < testDay; ++j){
      //Rcpp::Rcout << "Testing pregnancy on day: " << j << std::endl;
      /* If after switch time and before abortion date cutoff, 
	 prob is probability of getting infection (modified 
	 for additional avoidance), multiplied by the probability
	 of getting microcephaly given infection in that gestational week
	 multiplied by the probability of not getting an abortion */
      if((j >= switch_t) & (j - (testDay - probM_size)) < abortion_time){
	day_prob = riskI[j]*(1-birth_reduction)*(1 - abortion_rate)*probM[j-(testDay-probM_size)];
	/* If after switch time but not within abortion time, just modified infection risk */
      } else if(j >= switch_t){
	day_prob = riskI[j]*(1-birth_reduction)*probM[j-(testDay-probM_size)];
	/* Otherwise, just normal risk of getting infected multiplied by risk of
	   microcephaly given infection in that gestational week */
      } else {
	day_prob = riskI[j]*probM[j-(testDay-probM_size)];
      }
      tmp += day_prob;
    }
    allProbs[i] = 1 - (1-tmp)*(1-bp);
  }
  return(allProbs);
}


