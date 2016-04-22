#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>

#include "allvars.h"
#include "proto.h"

void get_cov_matrix(double *theta, int nstep, int ntheta)
{
  int i, j;
  double cov;
  
  for(i=0; i<ntheta; i++)
  {
    cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[i*nstep], 1, nstep);
    cov = cov>1.0e-6 ? cov : 1.0e-6;
    cov_matrix[i*ntheta + i] = cov;
    for(j=0; j<i; j++)
    {
      cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[j*nstep], 1, nstep);
      cov = cov>0.0001?cov:0.0;
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = cov;
    }
  }
}
void get_cov_matrix_diag(double *theta, int nstep, int ntheta)
{
  int i, j;
  double cov;
  
  for(i=0; i<ntheta; i++)
  {
    cov = gsl_stats_covariance(&theta[i*nstep], 1, &theta[i*nstep], 1, nstep);
    cov = cov>1.0e-6 ? cov : 1.0e-6;
    cov_matrix[i*ntheta + i] = cov;
    for(j=0; j<i; j++)
    {
      cov_matrix[i*ntheta + j] = cov_matrix[j*ntheta + i] = 0.0;
    }
  }
}

void mcmc_stats(char * fname)
{
  FILE *fp;
  double **theta, tmp;
  int i, nstep, istep;
  char buf[1000], *pstr, buf1[100];
  
  theta = malloc(ntheta*sizeof(double));
  for(i=0; i<ntheta; i++)
  {
    theta[i] = malloc((nmcmc)*sizeof(double));
  }

  fp = fopen(fname, "r");
  if(fp==NULL)
  {
    strcpy(str_error_exit, fname);
    error_exit(2);
  }

  istep = 0;
  nstep = 0;
  while(!feof(fp) && nstep<=nmcmc)
  {
    fscanf(fp, "%d", &nstep);
    for(i=0; i<ntheta; i++)
    {
      fscanf(fp, "%lf", &theta[i][istep]);
    }
    fscanf(fp,"%lf\n", &tmp);
    istep++;
  }
  printf("nstep: %d\n", istep);
  for(i=0; i<ntheta; i++)
  {
    gsl_sort(&theta[i][nbuilt], 1, istep-nbuilt);
    //printf("ss:%f %f %f\n", theta[i*n_mcmc + 40000], theta[i*n_mcmc+ 40100], theta[i*n_mcmc+ 40200]);
    theta_best[i] = gsl_stats_mean(&theta[i][nbuilt], 1, istep-nbuilt);
    theta_best_var[i*2] = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt], 1, istep-nbuilt, 0.1585);
    theta_best_var[i*2+1] = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt], 1, istep-nbuilt, 1.0-0.1585);

    theta_best_var[i*2] = theta_best[i] - theta_best_var[i*2];
    theta_best_var[i*2+1] -= theta_best[i];

    printf("%f %f %f\n", theta_best[i], theta_best_var[i*2], theta_best_var[i*2+1]);
  }

  for(i = 0; i<ntheta; i++)
  {
    free(theta[i]);
  }
  free(theta);

  fclose(fp);
  return;
}