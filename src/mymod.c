#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static double parms[13];
#define L_M parms[0]
#define D_EM parms[1]
#define L_H parms[2]
#define D_C parms[3]
#define D_F parms[4]
#define D_EH parms[5]
#define D_IH parms[6]
#define b parms[7]
#define P_HM parms[8]
#define P_MH parms[9]
#define t_seed parms[10]
#define I0 parms[11]
#define offset1 parms[12]

/* initializer */
void initmod(void (* odeparms)(int *, double *))
{
  int N=13;
  odeparms(&N, parms);
}
/* Derivatives and 1 output variable */
void derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
  double epsilon = 0.0001; // Constant to check for seeding event
  double B_H = L_H - D_C; // Rate of entering the first trimester class (ie. rate of getting pregnant)

  double S_M = y[0]; // Susceptible mosquitoes
  double S_C = y[3]; // Susceptible children 
  double S_A = y[4]; // Susceptible adults
  double S_F = y[5]; // Susceptible first trimester adults
 
  double E_M = 0; // Exposed mosquitoes
  double I_M = 0; // Infected mosquitoes
  
  double E_C = 0; // Exposed children
  double I_C = 0; // Infected children
  double E_A = 0; // Exposed adults
  double I_A = 0; // Infected adults
  double E_F = 0; // Exposed first trimester adults
  double I_F = 0; // Infected first trimester adults
  double R_C = 0; // Recovered children
  double R_A = 0; // Recovered adults
  double R_F = 0; // Recovered first trimester adults
  double offset = 0;

  E_M = y[1];
  I_M = y[2];
  E_C = y[6];
  E_A = y[7];
  E_F = y[8];
  I_C = y[9];
  I_A = y[10];
  I_F = y[11];
  R_C = y[12];
  R_A = y[13];
  R_F = y[14];
  double seed = 0;

  // Checks if passed the seeding time of the epidemic. If yes, add constant seeding to FOI.
  if(*t >= t_seed){
    offset = offset1;
    // Just after seeding time, give a pulse of initial infecteds
    if(*t < t_seed + epsilon){
      seed = I0;
    }
  }

  // Calculate total population sizes
  double N_H = S_C + S_A + S_F + E_C + I_C + E_A + I_A + E_F + I_F + R_C + R_A + R_F;
  double N_M = S_M + E_M + I_M;

  // Force of infection from humans to mosquitoes. Note that this comes from infected humans plus the initial seed and constant influx of infecteds
  double lambda_M = b*P_HM*(I_C + I_A + I_F + offset + seed)/N_H;
  // FOI from mosquitoes to humans
  double lambda_H = b*P_MH*(I_M)/N_H;

  
  // Changes in mosquito population
  ydot[0] = N_M/L_M - y[0]/L_M - lambda_M*y[0];
  ydot[1] = lambda_M*y[0] - E_M/L_M - E_M/D_EM;
  ydot[2] = E_M/D_EM - I_M/L_M;
  
  // Change in susceptible humans
  ydot[3] = N_H/L_H - y[3]/L_H - y[3]/D_C - lambda_H*y[3];
  ydot[4] = y[3]/D_C - y[4]/L_H - y[4]/B_H + y[5]/D_F - lambda_H*y[4];
  ydot[5] = y[4]/B_H - y[5]/D_F - y[5]/L_H - lambda_H*y[5];
  
  // Change of exposed humans
  ydot[6] = lambda_H*y[3] - E_C/L_H           - E_C/D_EH           - E_C/D_C;
  ydot[7] = lambda_H*y[4] - E_A/L_H - E_A/B_H - E_A/D_EH + E_F/D_F + E_C/D_C;
  ydot[8] = lambda_H*y[5] - E_F/L_H + E_A/B_H - E_F/D_EH - E_F/D_F;
  
  // Change of infected humans
  ydot[9] =  E_C/D_EH - I_C/L_H - I_C/D_IH           - I_C/D_C;
  ydot[10] = E_A/D_EH - I_A/L_H - I_A/D_IH - I_A/B_H + I_C/D_C + I_F/D_F;
  ydot[11] =  E_F/D_EH - I_F/L_H - I_F/D_IH + I_A/B_H           - I_F/D_F ;
  
  // Change of recovered humans
  ydot[12] = I_C/D_IH - R_C/L_H - R_C/D_C;
  ydot[13] = I_A/D_IH + R_C/D_C - R_A/L_H - R_A/B_H + R_F/D_F;
  ydot[14] = I_F/D_IH + R_A/B_H - R_F/D_F - R_F/L_H;
  
  // Saving movements out of pregnancy class for post processing
  ydot[15] = I_F/D_F;
  ydot[16] = y[5]/D_F + E_F/D_F + R_F/D_F + I_F/D_F;
 
	 
  //  ydot[17] = lambda_H*y[3] + lambda_H*y[4] + lambda_H*y[5];
  // Save human incidence
  ydot[17] = E_C/D_EH + E_A/D_EH + E_F/D_EH;

  yout[0] = I_F/D_F;
  yout[1] = I_F/D_IH;
  yout[2] = y[5]/D_F + E_F/D_F + R_F/D_F + I_F/D_F;
  yout[3] = E_C/D_EH + E_A/D_EH + E_F/D_EH;
 


}
