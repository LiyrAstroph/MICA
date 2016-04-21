#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>

#include "allvars.h"
#include "proto.h"


double probability_con(double sigmahat, double taud, double alpha)
{
  double prob, lambda, ave_con, lndet, lndet_ICq, sigma;
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

  return prob;  
}

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