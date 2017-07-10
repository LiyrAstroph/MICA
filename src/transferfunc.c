#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void transfer_function(double *theta)
{
  FILE *ftran, *ftran_comp;
  const int nt=ntau;
  int i, j;
  double tau, phi, phi1, phi2, phii, phierr, err1, err2, err3, dlnphii, w, fk;

  ftran = fopen("data/transfer.txt", "w");
  ftran_comp = fopen("data/transfer_comp.txt", "w");
  
  w = theta[2];
  for(i=0; i<nt; i++)
  {
    TF_tau[i] = (tau_lim_up-tau_lim_low)/(nt-1.0) * i + tau_lim_low;
  }
#ifdef TOPHAT
  for(i=0; i<nt; i++)
  {
    tau = TF_tau[i];
    phi = 0.0;
    phierr = 0.0;
    phi1 = phi2 = 0.0;
    fprintf(ftran, "%f ", tau);
    for(j=0; j<nc; j++)
    {
      if( tau >= grid_tau[j] - w &&  tau < grid_tau[j] + w )
      {
        phii = exp(theta[3+j])/(2.0*w);
        phi += phii;
        err1 = theta_best_var[(3+j)*2];
        dlnphii = err1;
        phi1 += exp(log(phii) - dlnphii);

        err1 = theta_best_var[(3+j)*2+1];
        dlnphii = err1;
        phi2 += exp(log(phii) + dlnphii);
        phierr += pow(phii * dlnphii, 2.0);
      }
    }
    fprintf(ftran, "%e %e %e\n", phi, phi1, phi2);
    TF[i] = phi;
  }
#elif defined JAVELIN
  for(i=0; i<nt; i++)
  {
    tau = TF_tau[i];
    phi = 0.0;
    phierr = 0.0;
    phi1 = phi2 = 0.0;
    fprintf(ftran, "%f ", tau);
    fk = theta[3];
    if( tau >= theta[4] - w &&  tau < theta[4] + w )
    {
      phii = fk/(2.0*w);
      phi += phii;
      err1 = theta_best_var[(3+j)*2];
      dlnphii = err1;
      phi1 += exp(log(phii) - dlnphii);
      err1 = theta_best_var[(3+j)*2+1];
      dlnphii = err1;
      phi2 += exp(log(phii) + dlnphii);
      phierr += pow(phii * dlnphii, 2.0);
    }
    fprintf(ftran, "%e %e %e\n", phi, phi1, phi2);
    TF[i] = phi;
  }
#else
  for(i=0; i<nt; i++)
  {
    tau = TF_tau[i];
    phi = 0.0;
    phierr = 0.0;
    phi1 = phi2 = 0.0;
    fprintf(ftran, "%f ", tau);
    fprintf(ftran_comp, "%f ", tau);
    for(j=0; j<nc; j++)
    {
      fk = exp(theta[3+j]);

      phii = fk/sqrt(2.0*PI)/w * exp( - 0.5* pow( (tau - grid_tau[j])/w, 2.0 ));
      phi += phii;

      err1 = theta_best_var[(3+j)*2];
      err2 = theta_best_var[2*2]/w;
      err3 = pow( (tau - grid_tau[j])/w, 2.0 ) * err2;
      dlnphii = sqrt(pow(err1, 2.0) + pow(err2, 2.0) + pow(err3, 2.0));
      phi1 += exp(log(phii) - dlnphii);

      //phi1 += exp(theta[3+j]-theta_best_var[(3+j)*2])/sqrt(2.0*PI)/(w )   
      //* exp( - 0.5* pow( (tau - grid_tau[j])/(w), 2.0 ));
     // phi2 += exp(theta[3+j]+theta_best_var[(3+j)*2+1])/sqrt(2.0*PI)/(w)   
      //* exp( - 0.5* pow( (tau - grid_tau[j])/(w), 2.0 ));

      err1 = theta_best_var[(3+j)*2+1];
      err2 = theta_best_var[2*2+1]/w;
      err3 = pow( (tau - grid_tau[j])/w, 2.0 ) * err2;
      dlnphii = sqrt(pow(err1, 2.0) + pow(err2, 2.0) + pow(err3, 2.0));
      phi2 += exp(log(phii) + dlnphii);
      phierr += pow(phii * dlnphii, 2.0);

      //printf("%e %e %e %e %e\n", phi2, exp(log(phii) + dlnphii), err1, err2, err3);
      fprintf(ftran_comp, "%e ", phii);
    }
    fprintf(ftran, "%e %e %e\n", phi, phi1, phi2);
    fprintf(ftran_comp, "\n");
    TF[i] = phi;
  }

  double tau_mean, tau_mean_err, tauk, fac1, fac2;

  fac1 = 0.0;
  fac2 = 0.0;

  err1=err2=0.0;
  for(j=0; j<nc; j++)
  {
    tauk = grid_tau[j];
    fk = exp(theta[3+j]);
    fac1 += tauk * fk;
    fac2 += fk;

    err1 += tauk*tauk * pow(fk*0.5*(theta_best_var[(3+j)*2] + theta_best_var[(3+j)*2 + 1]), 2.0);
    err2 += pow(fk*0.5*(theta_best_var[(3+j)*2] + theta_best_var[(3+j)*2 + 1]), 2.0);
  }

  tau_mean = fac1 / fac2;
  tau_mean_err = tau_mean * sqrt(err1/fac1/fac1 + err2/fac2/fac2);
  printf("tau_mean: %f %f\n", tau_mean, tau_mean_err);
  fprintf(fp_results, "tau_mean: %f %f\n", tau_mean, tau_mean_err);
#endif
  fclose(ftran);
  fclose(ftran_comp);
}
