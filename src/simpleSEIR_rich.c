int simpleSEIR_rich(double t, double *y, double *ydot, void *data)
{
  double *parms = (double*)data;
  double L_M = parms[0], L_H = parms[1], D_EM = parms[2], D_EH =  parms[3], D_IH = parms[4], 
    b = parms[5], p_HM = parms[6], p_MH = parms[7], seeding = parms[8];

  double S_M = y[0]; // Susceptible mosquitoes
  double E_M = y[1]; // Exposed mosquitoes
  double I_M = y[2]; // Infected mosquitoes

  
  double S_H = y[3]; // Susceptible humans
  double E_H = y[4]; // Exposed humans
  double I_H = y[5]; // Infected humans
  double R_H = y[6]; // Recovered humans

  double N_H = S_H + E_H + I_H + R_H;
  double N_M = S_M + E_M + I_M;
    
  double lambda_M = b*p_HM*(I_H)/N_H; // + 1.0/(7.0*N_H);
  double lambda_H = (b*p_MH*I_M/N_H);// + 1.0/(7.0*N_H);

  
  // Changes in mosquito population
  if(t < seeding){
    S_H = S_H + I_H;
    I_H = 0;
    ydot[0] = N_M/L_M - y[0]/L_M;
    ydot[1] = - E_M/L_M;
    ydot[2] = - I_M/L_M;
    
    ydot[3] = N_H/L_H - S_H/L_H;
    ydot[4] = - E_H/L_H;
    ydot[5] = - I_H/L_H;
    ydot[6] = - R_H/L_H;
    ydot[7] = 0;

  } else {
    ydot[0] = N_M/L_M - y[0]/L_M - lambda_M*y[0];
    ydot[1] = lambda_M*y[0] - E_M/L_M - E_M/D_EM;
    ydot[2] = E_M/D_EM - I_M/L_M;
    
    ydot[3] = N_H/L_H - S_H/L_H - lambda_H*S_H;
    ydot[4] = lambda_H*S_H - E_H/L_H - E_H/D_EH;
    ydot[5] = E_H/D_EH - I_H/L_H - I_H/D_IH;
    ydot[6] = I_H/D_IH - R_H/L_H;

    ydot[7] = E_H/D_EH - E_H/L_H;
  }
  return 0;
}
