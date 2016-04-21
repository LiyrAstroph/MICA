#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

int mcmc_sampling(char *fname_mcmc, double (* prob_fun)(double *theta))
{
  const int n_cov_update=10000;
  double prob, prob_new, scale, ratio, u;
  double theta[ntheta], theta_new[ntheta];
  double var_cov_mat[ntheta*ntheta], Pmatrix[ntheta*ntheta], Prandvec[ntheta], Py[ntheta];
  double *theta_mcmc;
  int istep, iaccpt, i, j, info, flag;

  FILE *fmcmc;
  fmcmc = fopen(fname_mcmc, "w");
  theta_mcmc = array_malloc(ntheta*n_cov_update);

/* 
 * sampling in logrithm space.
 */
  memcpy(theta, theta_input, ntheta*sizeof(double));
  
  for(i=0; i<ntheta; i++)
  {
    var_cov_mat[i*ntheta+i] = sigma_input[i] * sigma_input[i];
    for(j=0; j<i; j++)
      var_cov_mat[i*ntheta+j] = var_cov_mat[j*ntheta + i] = 0.0;
  }
  memcpy(Pmatrix, var_cov_mat, ntheta*ntheta*sizeof(double));
  Chol_decomp_U(Pmatrix, ntheta, &info);

  scale = 1.0;

  prob = prob_fun(theta);
  istep = 0;
  iaccpt = 0;
  while(istep < nmcmc)
  {
    flag = 1;
    while(flag)
    {
      for(i=0; i<ntheta; i++)
      {
        Prandvec[i]  = gsl_ran_gaussian(gsl_r, 1.0);
      }
      multiply_matvec(Pmatrix, Prandvec, ntheta, Py);
      flag = 0;
      for(i=0; i<ntheta; i++)
      {
        if(theta_fixed[i]==1)
        {
          theta_new[i] = theta[i];
        }
        else
        {
          theta_new[i] = theta[i] + scale * Py[i];
          if(theta_new[i]<theta_range[i][0] || theta_new[i]>theta_range[i][1])
          {
            flag = 1;
            break;
          }
        }
      }
    }
       
    prob_new = prob_fun(theta_new);
    
    ratio = prob_new - prob;
    if(ratio > 0.0)
    {
      memcpy(theta, theta_new, ntheta*sizeof(double));
      prob = prob_new;
      iaccpt++;
    }
    else
    {
      u = gsl_rng_uniform(gsl_r);
      if(log(u) < ratio)
      { 
        memcpy(theta, theta_new, ntheta*sizeof(double));
        prob = prob_new;
        iaccpt++;
      }
    }

    for(i=0; i<ntheta; i++)
      theta_mcmc[i*n_cov_update + istep%n_cov_update] = theta[i];
    if(istep%100==0)printf("%d\n", istep);

    fprintf(fmcmc, "%d", istep);
    for(i=0;i<ntheta;i++)
    {
      fprintf(fmcmc,"\t%f", theta[i]);
    }
    fprintf(fmcmc, "\t%e\n", prob);

    istep++;
    
    if(istep%n_cov_update == 0)
    {
      get_cov_matrix(theta_mcmc, n_cov_update, ntheta);
      memcpy(var_cov_mat, cov_matrix, ntheta*ntheta*sizeof(double));
      display_mat(var_cov_mat, ntheta, ntheta);
      memcpy(Pmatrix, var_cov_mat, ntheta*ntheta*sizeof(double));
      Chol_decomp_U(Pmatrix, ntheta, &info);
    }
  }

  printf("Accept rate: %f\n", 1.0*iaccpt/istep);

  fclose(fmcmc);
  free(theta_mcmc);
  return 0;
}