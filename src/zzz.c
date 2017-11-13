#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>

SEXP _zikaProj_toUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP);
SEXP _zikaProj_fromUnitScale(SEXP xSEXP, SEXP minSEXP, SEXP maxSEXP);
SEXP _zikaProj_generate_foi(SEXP IMSEXP, SEXP NHSEXP, SEXP bSEXP, SEXP pMHSEXP, SEXP tstepSEXP);
SEXP _zikaProj_generate_riskS(SEXP foiSEXP, SEXP tstepSEXP);
SEXP _zikaProj_generate_riskI(SEXP foiSEXP, SEXP riskSSEXP, SEXP tstepSEXP);
SEXP _zikaProj_generate_probM_aux(SEXP riskISEXP, SEXP probMSEXP, SEXP bpSEXP);
SEXP _zikaProj_generate_probM(SEXP IMSEXP, SEXP NHSEXP, SEXP probMSEXP, SEXP bSEXP, SEXP pMHSEXP, SEXP bpSEXP, SEXP tstepSEXP);
SEXP _zikaProj_likelihood_probM(SEXP microBirthsSEXP, SEXP allBirthsSEXP, SEXP probMSEXP);
SEXP _zikaProj_likelihood_probM_norm(SEXP microBirthsSEXP, SEXP allBirthsSEXP, SEXP probMSEXP, SEXP lik_sdSEXP);
SEXP _zikaProj_average_buckets(SEXP aSEXP, SEXP bucketsSEXP);
SEXP _zikaProj_sum_buckets(SEXP aSEXP, SEXP bucketsSEXP);
SEXP _zikaProj_generate_probM_forecast(SEXP riskISEXP, SEXP probMSEXP, SEXP bpSEXP, SEXP abortion_rateSEXP, SEXP birth_reductionSEXP, SEXP switch_tSEXP, SEXP abortion_timeSEXP, SEXP ABORTEDSEXP);
SEXP _zikaProj_generate_probM_forecast_OLD(SEXP riskISEXP, SEXP probMSEXP, SEXP bpSEXP, SEXP abortion_rateSEXP, SEXP birth_reductionSEXP, SEXP switch_tSEXP, SEXP abortion_timeSEXP);


static const R_CMethodDef cMethods[] = {
  // Disabled until the rlsoda registration issue is fixed
  //{"C_SEIR_model_lsoda", (DL_FUNC) &SEIR_model_lsoda, 4, NULL},
  //{"C_SEIR_model_rlsoda", (DL_FUNC) &SEIR_model_rlsoda, 4, NULL},
  {NULL, NULL, 0, NULL}  
};


 
static const R_CallMethodDef CallEntries[] = {
  {"_zikaProj_toUnitScale", (DL_FUNC) &_zikaProj_toUnitScale, 3},
  {"_zikaProj_fromUnitScale", (DL_FUNC) &_zikaProj_fromUnitScale, 3},
  {"_zikaProj_generate_foi", (DL_FUNC) &_zikaProj_generate_foi, 5},
  {"_zikaProj_generate_riskS", (DL_FUNC) &_zikaProj_generate_riskS, 2},
  {"_zikaProj_generate_riskI", (DL_FUNC) &_zikaProj_generate_riskI, 3},
  {"_zikaProj_generate_probM_aux", (DL_FUNC) &_zikaProj_generate_probM_aux, 3},
  {"_zikaProj_generate_probM", (DL_FUNC) &_zikaProj_generate_probM, 7},
  {"_zikaProj_likelihood_probM", (DL_FUNC) &_zikaProj_likelihood_probM, 3},
  {"_zikaProj_likelihood_probM_norm", (DL_FUNC) &_zikaProj_likelihood_probM_norm, 4},
  {"_zikaProj_average_buckets", (DL_FUNC) &_zikaProj_average_buckets, 2},
  {"_zikaProj_sum_buckets", (DL_FUNC) &_zikaProj_sum_buckets, 2},
  {"_zikaProj_generate_probM_forecast", (DL_FUNC) &_zikaProj_generate_probM_forecast, 8},
  {"_zikaProj_generate_probM_forecast_OLD", (DL_FUNC) &_zikaProj_generate_probM_forecast_OLD, 7},
  {NULL, NULL, 0}
};

void R_init_zikaProj(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, CallEntries, NULL, NULL);
  // Rich needs to work out how to get rlsoda to behave with registered routines, so this is disabled for now.
  //R_useDynamicSymbols(dll, FALSE);
  //R_forceSymbols(dll, TRUE);
}
