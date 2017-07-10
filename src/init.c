
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "allvars.h"
#include "proto.h"

void init()
{
  int i;
  double dT, T1, T2;

  strcpy(fname_results, "data/results.txt");
  fp_results = fopen(fname_results, "w");

#ifdef JAVELIN
  fprintf(fp_results, "JAVELIN-single top-hat transfer function.\n");
  printf("JAVELIN.\n");
  nc_lim_low = nc_lim_up = 1;
#elif defined TOPHAT
  fprintf(fp_results, "multiple top-hat transfer function.\n");
  printf("TOPHATs.\n");
#else
  fprintf(fp_results, "multiple Gaussian transfer function.\n");
  printf("GAUSSIANs.\n");
#endif

  memory_alloc_data();
  read_data();

  nall = ncon + nline;
  nall_data = ncon_data + nline_data;
  memory_alloc();

  scale_con = scale_line = 1.0;
  scale_light_curves();

  memcpy(Fall_data, Fcon_data, ncon_data*sizeof(double));
  memcpy(Fall_data+ncon_data, Fline_data, nline_data*sizeof(double));
  
  len_con = Tcon_data[ncon_data-1] - Tcon_data[0];
  len_line = Tline_data[nline_data-1] - Tline_data[0];

  cad_con = 0.0;
  for(i=0; i<ncon_data-1; i++)
  {
    cad_con += Tcon_data[i+1] - Tcon_data[i];
  }
  cad_con /= (ncon_data-1);

  cad_line = 0.0;
  for(i=0; i<nline_data-1; i++)
  {
    cad_line += Tline_data[i+1] - Tline_data[i];
  }
  cad_line /= (nline_data-1);

  printf("len: %f %f\n", len_con, len_line);
  printf("cad: %f %f\n", cad_con, cad_line);

  T1 = Tcon_data[0];
  T2 = Tcon_data[ncon_data-1];

  dT = T2-T1;
  T1 -= 0.1*dT;
  T2 += 0.1*dT;
  dT = (T2-T1)/(ncon-1.0);
  for(i=0; i<ncon; i++)
    Tcon[i] = T1 + dT *i;

  T1 = Tline_data[0];
  T2 = Tline_data[nline_data-1];

  dT = T2-T1;
  T1 -= 0.1*dT;
  T2 += 0.1*dT;
  dT = (T2-T1)/(nline-1.0);
  for(i=0; i<nline; i++)
    Tline[i] = T1 + dT *i;
}

void scale_light_curves()
{
  int i;
  double mean, norm;
  
  mean = 0.0;
  norm = 0.0;
  for(i=0; i<ncon_data; i++)
  {
    mean += Fcon_data[i]/(Fcerrs_data[i] * Fcerrs_data[i]);
    norm += 1.0/(Fcerrs_data[i] * Fcerrs_data[i]);
  }
  mean /= norm;
  scale_con = mean;

  for(i=0; i<ncon_data;i++)
  {
    Fcon_data[i] /= mean;
    Fcerrs_data[i] /= mean;
  }

  mean = 0.0;
  norm = 0.0;
  for(i=0; i<nline_data; i++)
  {
    mean += Fline_data[i]/(Flerrs_data[i]*Flerrs_data[i]);
    norm += 1.0/(Flerrs_data[i]*Flerrs_data[i]);
  }
  mean /= norm;
  scale_line = mean;
  for(i=0; i<nline_data;i++)
  {
    Fline_data[i] /= mean;
    Flerrs_data[i] /= mean;
  }
  
  printf("scale: %e %e\n", scale_con, scale_line);
  fprintf(fp_results, "scale: %e %e\n", scale_con, scale_line);
}

/*
 *
 */
void memory_alloc_data()
{
  Tcon_data = malloc(ndata_max*sizeof(double));
  if(Tcon_data==NULL)
  {
    strcpy(str_error_exit, "Tcon_data");
    error_exit(7);
  }

  Fcon_data = malloc(ndata_max*sizeof(double));
  if(Fcon_data==NULL)
  {
    strcpy(str_error_exit, "Fcon_data");
    error_exit(7);
  }

  Fcerrs_data = malloc(ndata_max*sizeof(double));
  if(Fcerrs_data==NULL)
  {
    strcpy(str_error_exit, "Fcerrs_data");
    error_exit(7);
  }

  Tline_data = malloc(ndata_max*sizeof(double));
  if(Tline_data==NULL)
  {
    strcpy(str_error_exit, "Tline_data");
    error_exit(7);
  }

  Fline_data = malloc(ndata_max*sizeof(double));
  if(Fline_data==NULL)
  {
    strcpy(str_error_exit, "Fline_data");
    error_exit(7);
  }

  Flerrs_data = malloc(ndata_max*sizeof(double));
  if(Flerrs_data==NULL)
  {
    strcpy(str_error_exit, "Flerrs_data");
    error_exit(7);
  }
}
void memory_alloc()
{
  workspace_ipiv = malloc(ndata_max*sizeof(double));
  workspace = malloc(10*ndata_max*sizeof(double) + 10);
  
  Fall_data = array_malloc(nall_data);

  Smat = array_malloc(nall_data*nall_data);
  Nmat = array_malloc(nall_data*nall_data);
  Cmat = array_malloc(nall_data*nall_data);
  ICmat = array_malloc(nall_data*nall_data);
  ICvmat = array_malloc(nall_data*nall_data);

  USmat = array_malloc(nall*nall_data);
  ASmat = array_malloc(nall*nall);

  Tmat1 = array_malloc(nall_data*nall_data);
  Tmat2 = array_malloc(nall_data*nall_data);
  Tmat3 = array_malloc(nall_data*nall_data);
  Tmat4 = array_malloc(nall_data*nall_data);

  Tcon = array_malloc(ncon);
  Fcon = array_malloc(ncon);
  Fcerrs = array_malloc(ncon);
  Tline = array_malloc(nline);
  Fline = array_malloc(nline);
  Flerrs = array_malloc(nline);

  theta_range = matrix_malloc(ntheta_max, 2);
  theta_best = array_malloc(ntheta_max);
  theta_best_var = array_malloc(ntheta_max*2);
  theta_input = array_malloc(ntheta_max);
  sigma_input = array_malloc(ntheta_max);
  theta_fixed = malloc(ntheta_max*sizeof(int));

  theta_best_con = array_malloc(ntheta_max);
  theta_best_var_con = array_malloc(ntheta_max*2);

  cov_matrix = array_malloc(ntheta_max * ntheta_max);

  grid_tau = array_malloc(nc_lim_up);
  TF_tau = array_malloc(ntau);
  TF = array_malloc(ntau);

  gsl_T = gsl_rng_default;
  gsl_r = gsl_rng_alloc (gsl_T);
  gsl_rng_set(gsl_r, time(NULL));
}
