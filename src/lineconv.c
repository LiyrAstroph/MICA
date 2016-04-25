#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_interp.h>

#include "allvars.h"
#include "proto.h"

void line_convolution()
{
   
  FILE *fp;
  int info, i, j, nq;
  double flux, err, tline, taup, tcon, fcon, fcon_err, dtau; 
  double *yrec_line;
  double sigma, taud;
  double *Larr, *ybuf, *Cq, *ICq, *ave, *yave;

  gsl_interp_accel *gsl_acc_conv, *gsl_acc_error_conv;
  gsl_interp  *gsl_linear_conv, *gsl_linear_error_conv;

  gsl_acc_conv = gsl_interp_accel_alloc();
  gsl_acc_error_conv = gsl_interp_accel_alloc();

  gsl_linear_conv = gsl_interp_alloc(gsl_interp_linear, ncon);
  gsl_linear_error_conv = gsl_interp_alloc(gsl_interp_linear, ncon);

  gsl_interp_init(gsl_linear_conv, Tcon, Fcon, ncon);
  gsl_interp_init(gsl_linear_error_conv, Tcon, Fcerrs, ncon);

  yrec_line = array_malloc(nline_data);
  
// ave
  taud = exp(theta_best[1]);
  sigma = exp(theta_best[0]) * sqrt(taud/2.0);
  
  nq = 2*(flag_detrend + 1);
  Larr = workspace;
  ybuf = Larr + nq*nall_data;
  Cq = ybuf + nall_data;
  ICq = Cq + nq*nq;
  ave = ICq + nq*nq;
  yave = ave + nall_data;

  set_covar_mat(theta_best);

  for(i=0; i<nall_data*nall_data; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }

  for(i=0; i<ncon_data; i++)
  {
    if(flag_detrend==0)
    {
      Larr[i*nq+0] = 1.0;
      Larr[i*nq+1] = 0.0;
    }
    else
    {
      Larr[i*nq + 0] = 1.0;
      Larr[i*nq + 1] = Tcon_data[i];
      Larr[i*nq + 2] = 0.0;
      Larr[i*nq + 3] = 0.0;
    }
  }
  for(i=0; i<nline_data; i++)
  {
    if(flag_detrend==0)
    {
      Larr[(i+ncon_data)*nq+0] = 0.0;
      Larr[(i+ncon_data)*nq+1] = 1.0;
    }
    else
    {
      Larr[(i+ncon_data)*nq + 0] = 0.0;
      Larr[(i+ncon_data)*nq + 1] = 0.0;
      Larr[(i+ncon_data)*nq + 2] = 1.0;
      Larr[(i+ncon_data)*nq + 3] = Tline_data[i];
    }
  }

  memcpy(ICmat, Cmat, nall_data*nall_data*sizeof(double));
  inverse_mat(ICmat, nall_data, &info);

  multiply_mat_MN(ICmat, Larr, Tmat1, nall_data, nq, nall_data);
  multiply_mat_MN_transposeA(Larr, Tmat1, ICq, nq, nq, nall_data);

  memcpy(Cq, ICq, nq*nq*sizeof(double));
  inverse_mat(Cq, nq, &info);
  
  multiply_mat_MN_transposeA(Larr, ICmat, Tmat1, nq, nall_data, nall_data);
  multiply_mat_MN(Cq, Tmat1, Tmat2, nq, nall_data, nq);
  multiply_mat_MN(Tmat2, Fall_data, ave, nq, 1, nall_data);
  multiply_mat_MN(Larr, ave, yave, nall_data, 1, nq);
  printf("%f %f\n", ave[0], ave[1]);


  dtau = TF_tau[1] - TF_tau[0];
  fp = fopen("data/line_conv.txt", "w");

  for(i=0; i<nline_data; i++)
  {
    flux = 0.0;
    err = 0.0;
    tline = Tline_data[i]; 
    for(j=0; j<ntau; j++)
    {
      taup = TF_tau[j];
      tcon = tline - taup;
      if(tcon >=Tcon[0] && tcon <= Tcon[ncon-1])
      {
        fcon = gsl_interp_eval(gsl_linear_conv, Tcon, Fcon, tcon, gsl_acc_conv);
        fcon_err = gsl_interp_eval(gsl_linear_error_conv, Tcon, Fcerrs, tcon, gsl_acc_error_conv);
        flux += fcon * TF[j];
        err += pow(fcon_err * TF[j] * dtau, 2);
      }
    }
    flux *= dtau;
    err = sqrt(err);
    yrec_line[i] = flux + yave[i+ncon_data];
    fprintf(fp, "%f  %e %e\n", tline, yrec_line[i]*scale_line, err*scale_line);
  }
  fclose(fp);

  gsl_interp_free(gsl_linear_conv);
  gsl_interp_free(gsl_linear_error_conv);
  gsl_interp_accel_free(gsl_acc_conv);
  gsl_interp_accel_free(gsl_acc_error_conv);
}