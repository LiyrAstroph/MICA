#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_histogram.h>
#include <mpfit.h>


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

/*
 * For mpfit
 */
struct vars_struct
{
  double *x, *y, *ey;
};
int fitfunc(int m, int n, double *p, double *dy, double **devc, void *vars);


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

  const int  nh=30;
  int j;
  double hmax, hmin;
  gsl_histogram *h[ntheta];
  for(i=0; i<ntheta; i++)
  {
    h[i] = gsl_histogram_alloc(nh);
    hmax = gsl_stats_max(&theta[i][nbuilt], 1, istep-nbuilt);
    hmin = gsl_stats_min(&theta[i][nbuilt], 1, istep-nbuilt);
    gsl_histogram_set_ranges_uniform(h[i], hmin, hmax);

    for(j = nbuilt; j<nmcmc; j++)
    {
      gsl_histogram_increment(h[i], theta[i][j]);
    }
  }

  
/*  struct vars_struct v;
  mp_result result;
  mp_config config;
  mp_par pars[3];
  double p[3], perr[3];
  int status;

  config.maxiter = 1000;
  
  for(i=0; i<ntheta; i++)
  {
    hmin = h[i]->range[0];
    hmax = h[i]->range[nh]; 

    v.x = h[i]->range;
    v.y = h[i]->bin;
    memset(&result, 0, sizeof(result));
    result.xerror = perr;
    memset(&pars[0], 0, sizeof(pars));

    pars[0].limited[0] =  1;
    pars[0].limits[0] = 0.0;

    pars[1].limited[0] =  1;
    pars[1].limits[0] = hmin;
    pars[1].limited[1] =  1;
    pars[1].limits[1] = hmax;

    pars[2].limited[0] =  1;
    pars[2].limits[0] = 0.0;
    //pars[2].limited[1] =  1;
    //pars[2].limits[1] = (hmax - hmin);

    p[1] = gsl_histogram_mean(h[i]);
    p[2] = gsl_histogram_sigma(h[i]);
    p[0] = gsl_histogram_max_val(h[i]);


    status = mpfit(fitfunc, nh, 3, p, pars, &config, &v, &result);

    printf("Fit:%f %f %f\n", p[0], p[1], p[2]);
    theta_best[i] = p[1];
    theta_best_var[i*2] = theta_best_var[i*2+1] = p[2];
  } */

  for(i = 0; i<ntheta; i++)
  {
    free(theta[i]);
  }
  free(theta);

  fclose(fp);
  return;
}

int fitfunc(int m, int n, double *p, double *dy, double **devc, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey, f;

  x = v->x;
  y = v->y;
  //ey = v->ey;

  for (i=0; i<m; i++) {
    f = p[0] * exp(-0.5*pow(x[i] - p[1], 2.0)/p[2]/p[2]);
    dy[i] = (y[i] - f)/1.0;
  }
  return 0;
}