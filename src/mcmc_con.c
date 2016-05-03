#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

void mcmc_con_run()
{
  char fname_mcmc[100];
  strcpy(fname_mcmc, "data/mcmc_con.txt");

  mcmc_con_init();
  if(flag_mcmc==1)
  {
    mcmc_sampling(fname_mcmc, &probability_con);
    mcmc_stats(fname_mcmc);
  }
  else
  {
    read_input();
  }
  
  reconstruct_con();

  memcpy(theta_best_con, theta_best, ntheta*sizeof(double));
  memcpy(theta_best_var_con, theta_best_var, ntheta*2*sizeof(double));
}

void mcmc_con_init()
{
  int i;

  ntheta = 2;

  i=0;
  theta_range[i][0] = log(1.0e-6);
  theta_range[i++][1] = log(10.0);
  theta_range[i][0] = log(1.0e-2);
  theta_range[i++][1] = log(1.0e5);

  i = 0;
  theta_input[i++] = log(0.1);
  theta_input[i++] = log(100.0);

  i = 0;
  sigma_input[i++] = 0.1;
  sigma_input[i++] = 0.5;

  i = 0;
  theta_fixed[i++] = 0;
  theta_fixed[i++] = 0;

}
/*
 * probability for continuum mcmc.
 */
double probability_con(double *theta)//double sigmahat, double taud, double alpha)
{
  double sigmahat, taud, alpha;
  double prob, prior, lndet, lndet_ICq, sigma;
  double *ybuf, *Larr, *Cq, *ICq, *ave, *yave, *ysub;
  int i, nq, info;
  
  nq = (flag_detrend + 1);
  ybuf = workspace;
  Larr = workspace+ncon_data;
  Cq = Larr + nq*ncon_data;
  ICq = Cq + nq*nq;
  ave = ICq + nq*nq;
  yave = ave + nq;
  ysub = yave + ncon_data;
  
  alpha = 1.0;
  sigmahat = exp(theta[0]);
  taud = exp(theta[1]);

  sigma = sigmahat * sqrt(taud / 2.0);

  set_covar_mat_con(sigma, taud, alpha);
  
  for(i=0;i<ncon_data*ncon_data; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }
  memcpy(ICmat, Cmat, ncon_data*ncon_data*sizeof(double));
  inverse_mat(ICmat, ncon_data, &info);
  
  for(i=0;i<ncon_data;i++)
  {
    if(flag_detrend==0)
    {
      Larr[i]=1.0;
    }
    else
    {
      Larr[i*nq + 0] = 1.0;
      Larr[i*nq + 1] = Tcon_data[i];
    }
  }
  //multiply_matvec(ICmat, Larr, ncon_data, ybuf);
  multiply_mat_MN(ICmat, Larr, Tmat1, ncon_data, nq, ncon_data);
  multiply_mat_MN_transposeA(Larr, Tmat1, ICq, nq, nq, ncon_data);
  memcpy(Cq, ICq, nq*nq*sizeof(double));
  inverse_mat(Cq, nq, &info);

  multiply_mat_MN_transposeA(Larr, ICmat, Tmat1, nq, ncon_data, ncon_data);
  multiply_mat_MN(Cq, Tmat1, Tmat2, nq, ncon_data, nq);
  multiply_mat_MN(Tmat2, Fcon_data, ave, nq, 1, ncon_data);
  multiply_mat_MN(Larr, ave, yave, ncon_data, 1, nq);

//  lambda = cblas_ddot(ncon_data, Larr, 1, ybuf, 1);
//  multiply_matvec(ICmat, Fcon_data, ncon_data, ybuf);
//  ave_con = cblas_ddot(ncon_data, Larr, 1, ybuf, 1);
//  ave_con /=lambda;
  
  for(i=0;i<ncon_data;i++)ysub[i] = Fcon_data[i] - yave[i];
  multiply_matvec(ICmat, ysub, ncon_data, ybuf);
  prob = -0.5 * cblas_ddot(ncon_data, ysub, 1, ybuf, 1);

  lndet = lndet_mat(Cmat, ncon_data, &info);
  lndet_ICq = lndet_mat(ICq, nq, &info);
  prob = prob - 0.5*lndet - 0.5*lndet_ICq;

/* penalize on larger tau than the length of continuum */ 
  prior = 0.0; 
  prior += - theta[0];
  if(theta[1] > log(cad_con) )
  {
    prior = ( log(cad_con) - theta[1]);
  }
  else
  {
    prior +=  - (log(cad_con) - theta[1]);
  }
  return prob + prior;  
}

/* calculate covariance matrix */
void set_covar_mat_con(double sigma, double tau, double alpha)
{
  double t1, t2, nerr;
  int i, j;
  
  for(i=0; i<ncon_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<=i; j++)
    {

      t2 = Tcon_data[j];
      Smat[i*ncon_data+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha));
      Smat[j*ncon_data+i] = Smat[i*ncon_data+j];

      Nmat[i*ncon_data+j] = Nmat[j*ncon_data+i] = 0.0;
    }
    nerr = Fcerrs_data[i];
    Nmat[i*ncon_data+i] = nerr*nerr;
  }
  return;
}

void set_covar_Umat_con(double sigma, double tau, double alpha)
{
  double t1, t2;
  int i, j;
 
  for(i=0; i<ncon; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<ncon_data; j++)
    {
      t2 = Tcon_data[j];
      USmat[i*ncon_data+j] = sigma*sigma * exp (- pow (fabs(t1-t2) / tau, alpha) );
    }
  }
  return;
}

/* reconstruction */
void reconstruct_con()
{
  double sigma, tau, alpha, sigmahat;
  double *ybuf, *Larr, *Larr_rec, *yave, *yave_rec, *Cq, *ICq, *ysub, *ave;
  double *Tmp1, *Tmp2, *Tmp3, *Tmp4;

  int i, nq, info, nmax;

  nq = 1+flag_detrend;

  yave = workspace;
  ysub = yave + ncon_data;
  Larr_rec = ysub + ncon_data;    //Larr_rec is a nq*ncon matrix
  yave_rec = Larr_rec + nq*ncon;  //yave_rec is ncon matrix
  ybuf = yave_rec + ncon;
  Larr = ybuf + ncon_data;   // Larr is a nq*n matrix
  Cq = Larr + nq*ncon_data;
  ICq = Cq + nq*nq;
  ave = ICq + nq*nq;

  nmax = max(ncon, ncon_data);
  Tmp1 = array_malloc(nmax*nmax);     // temporary matrixes
  Tmp2 = array_malloc(nmax*nmax);
  Tmp3 = array_malloc(nmax*nmax);
  Tmp4 = array_malloc(nmax*nmax);

  sigmahat = exp(theta_best[0]);
  tau = exp(theta_best[1]);
  sigma = sigmahat * sqrt(tau/2.0);
  alpha = 1.0;

  set_covar_mat_con(sigma, tau, alpha);
  set_covar_Umat_con(sigma, tau, alpha);

  for(i=0;i<ncon_data*ncon_data; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }
  memcpy(ICmat, Cmat, ncon_data*ncon_data*sizeof(double));
  inverse_mat(ICmat, ncon_data, &info);

  /* the best estimate for average */
  for(i=0;i<ncon_data;i++)
  {
    if(flag_detrend==0)
    {
      Larr[i]=1.0;
    }
    else
    {
      Larr[i*nq + 0] = 1.0;
      Larr[i*nq + 1] = Tcon_data[i];
    }
  }

  for(i=0;i<ncon;i++)
  {
    if(flag_detrend==0)
    {
      Larr_rec[i]=1.0;
    }
    else
    {
      Larr_rec[i*nq + 0] = 1.0;
      Larr_rec[i*nq + 1] = Tcon[i];
    }
  }
  
// calculate the covariance matrix Cq for q.
// Cq = (L^T x C^-1 x L)^-1
  multiply_mat_MN(ICmat, Larr, Tmat1, ncon_data, nq, ncon_data);
  multiply_mat_MN_transposeA(Larr, Tmat1, ICq, nq, nq, ncon_data);
  memcpy(Cq, ICq, nq*nq*sizeof(double));
  inverse_mat(Cq, nq, &info);
// calculate the best estimate of q.
// q = Cq x L^T x C^-1 x y 
  multiply_mat_MN_transposeA(Larr, ICmat, Tmat1, nq, ncon_data, ncon_data);
  multiply_mat_MN(Cq, Tmat1, Tmat2, nq, ncon_data, nq);
  multiply_mat_MN(Tmat2, Fcon_data, ave, nq, 1, ncon_data);
// subtract the linear trend   
  multiply_mat_MN(Larr, ave, yave, ncon_data, 1, nq); 
  for(i=0;i<ncon_data;i++)ysub[i] = Fcon_data[i] - yave[i];
// calculate the best estimate for s.
// s = S x C^-1 x (y - Lxq)    
  multiply_matvec(ICmat, ysub, ncon_data, ybuf);
  multiply_matvec_MN(USmat, ncon, ncon_data, ybuf, Fcon);
  multiply_mat_MN(Larr_rec, ave, yave_rec, ncon, 1, nq);

  for(i=0; i<ncon; i++)Fcon[i] = Fcon[i] + yave_rec[i];

  multiply_mat_MN(USmat, ICmat, Tmp1, ncon, ncon_data, ncon_data);
  multiply_mat_MN_transposeB(Tmp1, USmat, Tmp2, ncon, ncon, ncon_data);
  multiply_mat_MN(Tmp1, Larr, Tmp3, ncon, nq, ncon_data);
  for(i=0; i<ncon*nq; i++)Tmp3[i] -= Larr_rec[i];
  multiply_mat_MN(Tmp3, Cq, Tmp1, ncon, nq, nq);
  multiply_mat_MN_transposeB(Tmp1, Tmp3, Tmp4, ncon, ncon, nq);

  for(i=0; i<ncon; i++)
  {
    Fcerrs[i] = sqrt(sigma*sigma - Tmp2[i*ncon+i] + Tmp4[i*ncon+i]);
  }

  FILE *fp;
  fp = fopen("data/scon.txt", "w");
  for(i=0; i<ncon; i++)
  {
    fprintf(fp, "%f %f %f\n", Tcon[i], Fcon[i]*scale_con, Fcerrs[i]*scale_con);
  }
  fclose(fp);
  free(Tmp1);
  free(Tmp2);
  free(Tmp3);
  free(Tmp4);
}