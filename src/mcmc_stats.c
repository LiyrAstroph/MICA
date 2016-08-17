#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>


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
  double **theta, tmp, err_max1, err_max2, theta_up, theta_low;
  int i, nstep, istep;

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
  while(!feof(fp) && nstep < nmcmc)
  {
    fscanf(fp, "%d", &istep);
    for(i=0; i<ntheta; i++)
    {
      fscanf(fp, "%lf", &theta[i][nstep]);
    }
    fscanf(fp,"%lf\n", &tmp);
    nstep++;
  }
  printf("nstep: %d\n", nstep);
  if(nstep <= nbuilt)
  {
    strcpy(str_error_exit, "mcmc.txt");
    error_exit(10);
  }

  printf("Sts: ID    par   err1   err2\n");
  for(i=0; i<ntheta; i++)
  {
    gsl_sort(&theta[i][nbuilt-1], 1, nstep-nbuilt);
    //printf("ss:%f %f %f\n", theta[i*n_mcmc + 40000], theta[i*n_mcmc+ 40100], theta[i*n_mcmc+ 40200]);
    theta_best[i] = gsl_stats_mean(&theta[i][nbuilt-1], 1, nstep-nbuilt);
    //theta_best[i] = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt-1], 1, nstep-nbuilt, 0.500);
    theta_low = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt-1], 1, nstep-nbuilt, 0.1585);
    theta_up = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt-1], 1, nstep-nbuilt, 1.0-0.1585);

    err_max1 = (theta_best[i] - theta_range[i][0])/2.35;
    err_max2 = (theta_range[i][1] - theta_best[i])/2.35;


    theta_best_var[i*2] = max(min( theta_best[i] - theta_low, err_max1 ), 0.0);
    theta_best_var[i*2+1] = max( min( theta_up - theta_best[i], err_max2 ), 0.0);


    printf("Sts: %d %f %f %f\n", i, theta_best[i], theta_best_var[i*2], theta_best_var[i*2+1]);

  }
  
  fclose(fp);
  return;
}