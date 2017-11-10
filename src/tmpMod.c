int SIR_model(double t, double *y, double *dydt, void *data){
  double *parms = (double*) data;
  double beta = parms[0];
  double gamma = parms[1];

  double S = y[0];
  double I = y[1];
  double R = y[2];

  dydt[0] = -beta*y[0]*y[1];
  dydt[1] = beta*y[0]*y[1] - gamma*y[1];
  dydt[2] = gamma*y[2];
  return 0;
}
