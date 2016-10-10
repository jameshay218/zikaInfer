#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "rlsoda.h"

// Registration:
static const R_CallMethodDef call_methods[] = {
  {"Clsoda",  (DL_FUNC) &r_lsoda,  16},
  {NULL,      NULL,                 0}
};
void R_init_rlsoda(DllInfo *info) {
  R_registerRoutines(info, NULL, call_methods, NULL, NULL);
}
