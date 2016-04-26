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


const int  nh=30;
void mcmc_stats(char * fname)
{
  FILE *fp;
  double **theta, tmp, err_max1, err_max2;
  int i, j, nstep, istep;
  char buf[1000], *pstr, buf1[100];
  double pars[3];

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

  for(i=0; i<ntheta; i++)
  {
    gsl_sort(&theta[i][nbuilt-1], 1, nmcmc-nbuilt);
    //printf("ss:%f %f %f\n", theta[i*n_mcmc + 40000], theta[i*n_mcmc+ 40100], theta[i*n_mcmc+ 40200]);
    theta_best[i] = gsl_stats_mean(&theta[i][nbuilt-1], 1, nmcmc-nbuilt);
    theta_best_var[i*2] = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt-1], 1, nmcmc-nbuilt, 0.1585);
    theta_best_var[i*2+1] = gsl_stats_quantile_from_sorted_data(&theta[i][nbuilt-1], 1, nmcmc-nbuilt, 1.0-0.1585);

    err_max1 = (theta_best[i] - theta_range[i][0])/2.35;
    err_max2 = (theta_range[i][1] - theta_best[i])/2.35;

    theta_best_var[i*2] = theta_best[i] - theta_best_var[i*2];
    theta_best_var[i*2+1] -= theta_best[i];

    theta_best_var[i*2] = min( theta_best_var[i*2], err_max1 );
    theta_best_var[i*2+1] = min( theta_best_var[i*2+1], err_max2 );


    printf("Sts: %d %f %f %f\n", i, theta_best[i], theta_best_var[i*2], theta_best_var[i*2+1]);

    /*par_fit(&theta[i][nbuilt-1], nmcmc-nbuilt, pars);

    theta_best[i] = pars[1];

    err_max1 = (theta_best[i] - theta_range[i][0])/2.35;
    err_max2 = (theta_range[i][1] - theta_best[i])/2.35;

    theta_best_var[i*2] = min( pars[2], err_max1);
    theta_best_var[i*2+1] = min(pars[2], err_max2);
    printf("Fit: %d %f %f %f %f\n", i, theta_best[i], theta_best_var[i*2], theta_best_var[i*2+1], err_max1);*/

  }
  
  fclose(fp);
  return;
}

/*
 * For mpfit
 */
struct vars_struct
{
  double *x, *y;// *ey;
};

int par_fit(double *theta, int n, double *pfit)
{
  int i, j;
  double hmax, hmin;
  gsl_histogram *hist = gsl_histogram_alloc(nh);
  
  struct vars_struct v;
  mp_result result;
  mp_config fitconfig;
  mp_par pars_set[3];
  double pars[3], perr[3], x[nh], y[nh];
  int status;

// use Gaussian to fit the distribution
  gsl_stats_minmax(&hmin, &hmax, theta, 1, n);
  gsl_histogram_set_ranges_uniform(hist, hmin, hmax);
  for(j = 0; j<n; j++)
  {
    gsl_histogram_increment(hist, theta[j]);
  }
    
  memcpy(x, hist->range, nh*sizeof(double));
  memcpy(y, hist->bin, nh*sizeof(double));
  v.x = x;
  v.y = y;
  memset(&fitconfig, 0, sizeof(fitconfig));
  memset(&result, 0, sizeof(result));
  result.xerror = perr;
  memset(&pars_set[0], 0, sizeof(pars_set));
    
  fitconfig.maxiter = 100;
    
  pars_set[0].limited[0] =  1;
  pars_set[0].limits[0] = 1.0e-5;

  pars_set[1].limited[0] =  1;
  pars_set[1].limits[0] = hmin;
  pars_set[1].limited[1] =  1;
  pars_set[1].limits[1] = hmax;

  pars_set[2].limited[0] =  1;
  pars_set[2].limits[0] = 1.0e-5;
    //pars[2].limited[1] =  1;
    //pars[2].limits[1] = (hmax - hmin);

  pars[1] = gsl_histogram_mean(hist);
  pars[2] = gsl_histogram_sigma(hist);
  pars[0] = gsl_histogram_max_val(hist);

  status = mpfit(&fitfunc, nh, 3, pars, pars_set, &fitconfig, (void *)&v, &result);
  
  memcpy(pfit, pars, 3*sizeof(double));
  
  gsl_histogram_free(hist);
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