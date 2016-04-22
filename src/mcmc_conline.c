#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_randist.h>

#include "allvars.h"
#include "proto.h"

void mcmc_conline_run()
{
  char fname_mcmc[100];
  strcpy(fname_mcmc, "data/mcmc.txt");

  mcmc_conline_init();
//  mcmc_sampling(fname_mcmc, &probability_conline);
  mcmc_stats(fname_mcmc);
  reconstruct_conline();
  transfer_function(theta_best);
  aicc();
  line_convolution();
}

void mcmc_conline_init()
{
  int i, j;

  
#ifdef JAVELINE  // only one-tophat is used in JAVELINE
  nc = 1;
  ntheta = 2 + 3 + 1;

  theta_range[0][0] = log(1.0e-6);
  theta_range[0][1] = log(10.0);

  theta_range[1][0] = log(1.0);
  theta_range[1][1] = log(1.0e5);

  theta_range[2][0] = (tau_lim_up - tau_lim_low)/1000.0;
  theta_range[2][1] = (tau_lim_up - tau_lim_low)*10.0;
  
  theta_range[3][0] = log(1.0e-10);
  theta_range[3][1] = log(100.0);

  theta_range[4][0] = tau_lim_low;
  theta_range[4][1] = tau_lim_up;

  theta_range[5][0] = log(1.0e-5);
  theta_range[5][1] = log(1.0e6);

  i = 0;
  sigma_input[i++] = 0.01; // sigma
  sigma_input[i++] = 0.01; // taud
  sigma_input[i++] = 0.01; // width
  sigma_input[i++] = 0.1;  // fk
  sigma_input[i++] = 0.1;  // tauc
  sigma_input[i++] = 0.1;

  theta_input[0] = theta_best[0];
  theta_input[1] = theta_best[1];
  theta_input[2] = (tau_lim_up-tau_lim_low)/5.0;
  theta_input[4] = (tau_lim_up-tau_lim_low)/10.0;
  theta_input[3] = 2.0/theta_input[4];
  theta_input[5] = log(10.0);
  
#else

  ntheta = nc + 3 + 1;
// grid of time lag
  for(i=0; i<nc; i++)
  {
    grid_tau[i] = (tau_lim_up-tau_lim_low)/(nc-1.0) * i;
  }

// limit range for parameters
  i=0;
  theta_range[i][0] = log(1.0e-6);
  theta_range[i++][1] = log(10.0);

  theta_range[i][0] = log(1.0);
  theta_range[i++][1] = log(1.0e4);

  theta_range[i][0] = (tau_lim_up - tau_lim_low)/(nc-1.0)/10.0;
  theta_range[i++][1] = (tau_lim_up - tau_lim_low)/(nc-1.0);
  for(j=0; j<nc; j++)
  {
    theta_range[i][0] = log(1.0e-10);
    theta_range[i++][1] = log(100.0);
  }
  
  theta_range[i][0] = log(1.0e-5);
  theta_range[i++][1] = log(1.0e6);

/* input step sizes for mcmc sampling */
  i = 0;
  sigma_input[i++] = 0.01; // sigma
  sigma_input[i++] = 0.01; // taud
  sigma_input[i++] = 0.01; // width
  for(j=0; j<nc; j++)
  {
    sigma_input[i++] = 0.1;  // fk
  } 
  sigma_input[i++] = 0.1;    // systematic error

/* input values for parameters */
  theta_input[0] = theta_best[0];
  theta_input[1] = theta_best[1];
  theta_input[2] = (tau_lim_up-tau_lim_low)/(nc-1.0)/2.0;
  for(i=0; i<nc; i++)
  {
    theta_input[3+i] = log(1.0/nc);
  }
  theta_input[3 + nc] = log(10.0);

#endif


/* set if the parameters are fixed */
  for(i=0; i<ntheta; i++)
  {
    theta_fixed[i] = 0;
  }

/* for TOPHAT case, the width is fixed */
#ifdef TOPHAT  
  theta_fixed[2] = 1;
#endif

}

/* theta represets sigma, tuad, and the parameters for transfer function */
double probability_conline(double *theta)
{
  double prob, prior, lndet_C, lndet_ICq, sigma, taud;
  double *Larr, *ybuf, *Cq, *ICq, *ave, *ysub, *yrec, *yrec_err, *yave;
  int i, nq, info, sign_C, sign_Cq;
  
  taud = exp(theta[1]);
  sigma = exp(theta[0]) * sqrt(taud/2.0);

  nq = (flag_detrend + 1)*2;
  Larr = workspace;
  ybuf = Larr + nq*nall_data;
  Cq = ybuf + nall_data;
  ICq = Cq + nq*nq;
  ave = ICq + nq*nq;

  ysub = ave + nq;
  yrec = ysub + nall_data;
  yrec_err = yrec + nall_data;
  yave = yrec_err + nall_data;

  set_covar_mat(theta);

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

// cal q
  memcpy(Tmat1, Cmat, nall_data*nall_data*sizeof(double));
  memcpy(Tmat2, Larr, nall_data*nq*sizeof(double));
  //printf("FFFF\n");
  //display_mat(Tmat1, 1, nall_data);
  multiply_mat_MN_inverseA(Tmat1, Tmat2, nall_data, nq); // Tmat2 = C^-1 * L
  
  multiply_mat_MN_transposeA(Larr, Tmat2, ICq, nq, nq, nall_data); // ICq = L^T*C^-1*L
  multiply_mat_MN_transposeA(Tmat2, Fall_data, ave, nq, 1, nall_data); // ave = L^T*C^-1*y
  memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ave, nq, 1); // (L^T*C^-1*L)^-1 * L^T*C^-1*y

  multiply_mat_MN(Larr, ave, yave, nall_data, 1, nq);
  for(i=0; i<nall_data; i++)ysub[i] = Fall_data[i] - yave[i];
  memcpy(Tmat1, Cmat, nall_data*nall_data*sizeof(double));
  memcpy(ybuf, ysub, nall_data*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ybuf, nall_data, 1); 

  prob = -0.5*cblas_ddot(nall_data, ysub, 1, ybuf, 1) / (sigma*sigma);
  //printf("%f\n", prob);
  if(prob > 0.0 )  // check if prob is positive
  {
    prob = -1.0e10;
    printf("prob >0!\n");
    return prob;
  }

  //memcpy(Tmat1, Cmat, n_data*nall_data*sizeof(double));
  lndet_C = lndet_mat3(Cmat, nall_data, &info, &sign_C) + 2.0*nall_data * log(sigma);
  if(info!=0|| sign_C==-1)
  {
    prob = -1.0e10;
    printf("lndet_C %f %d!\n", lndet_C, sign_C);
    return prob;
  }

  //memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  lndet_ICq = lndet_mat3(ICq, nq, &info, &sign_Cq) - 2.0*nq*log(sigma);
  if(info!=0 || sign_Cq==-1 )
  {
    prob = -1.0e10;
    printf("lndet_ICq!\n");
    return prob;
  }
  
  prob -= 0.5*(lndet_C + lndet_ICq);

  prior = -0.5*pow(theta[0] - theta_best_con[0], 2.0)/pow(max(theta_best_var_con[0*2], theta_best_var_con[0*2+1]), 2.0)
          -0.5*pow(theta[1] - theta_best_con[1], 2.0)/pow(max(theta_best_var_con[1*2], theta_best_var_con[1*2+1]), 2.0);

  prob += prior;
  return prob;
}

double aicc()
{
  double prob, aic, aicc;
  int k, n;

  k =  nc + 4;
  n = nall_data;

  prob = probability_conline_aicc(theta_best);

  aic = 2.0*k - 2.0*prob;

  aicc = aic + 2.0*k*(k+1.0)/(n - k - 1.0);
  //aicc = -2.0*prob + k * log(n); 
  printf("aicc: %d %f\n", nc, aicc);
  return aicc;
}

double probability_conline_aicc(double *theta)
{
  double prob,  lndet_C, lndet_ICq, sigma, taud;
  double *Larr, *ybuf, *Cq, *ICq, *ave, *ysub, *yrec, *yrec_err, *yave;
  int i, nq, info, sign_C, sign_Cq;
  
  taud = exp(theta[1]);
  sigma = exp(theta[0]) * sqrt(taud/2.0);

  nq = (flag_detrend + 1)*2;
  Larr = workspace;
  ybuf = Larr + nq*nall_data;
  Cq = ybuf + nall_data;
  ICq = Cq + nq*nq;
  ave = ICq + nq*nq;

  ysub = ave + nq;
  yrec = ysub + nall_data;
  yrec_err = yrec + nall_data;
  yave = yrec_err + nall_data;

  set_covar_mat(theta);

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

// cal q
  memcpy(Tmat1, Cmat, nall_data*nall_data*sizeof(double));
  memcpy(Tmat2, Larr, nall_data*nq*sizeof(double));
  //printf("FFFF\n");
  //display_mat(Tmat1, 1, nall_data);
  multiply_mat_MN_inverseA(Tmat1, Tmat2, nall_data, nq); // Tmat2 = C^-1 * L
  
  multiply_mat_MN_transposeA(Larr, Tmat2, ICq, nq, nq, nall_data); // ICq = L^T*C^-1*L
  multiply_mat_MN_transposeA(Tmat2, Fall_data, ave, nq, 1, nall_data); // ave = L^T*C^-1*y
  memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ave, nq, 1); // (L^T*C^-1*L)^-1 * L^T*C^-1*y

  multiply_mat_MN(Larr, ave, yave, nall_data, 1, nq);
  for(i=0; i<nall_data; i++)ysub[i] = Fall_data[i] - yave[i];
  memcpy(Tmat1, Cmat, nall_data*nall_data*sizeof(double));
  memcpy(ybuf, ysub, nall_data*sizeof(double));
  multiply_mat_MN_inverseA(Tmat1, ybuf, nall_data, 1); 

  prob = -0.5*cblas_ddot(nall_data, ysub, 1, ybuf, 1) / (sigma*sigma);
  //printf("%f\n", prob);
  if(prob > 0.0 )  // check if prob is positive
  {
    prob = -1.0e10;
    printf("prob >0!\n");
    return prob;
  }

  //memcpy(Tmat1, Cmat, n_data*nall_data*sizeof(double));
  lndet_C = lndet_mat3(Cmat, nall_data, &info, &sign_C) + 2.0*nall_data * log(sigma);
  if(info!=0|| sign_C==-1)
  {
    prob = -1.0e10;
    printf("lndet_C %f %d!\n", lndet_C, sign_C);
    return prob;
  }

  //memcpy(Tmat1, ICq, nq*nq*sizeof(double));
  lndet_ICq = lndet_mat3(ICq, nq, &info, &sign_Cq) - 2.0*nq*log(sigma);
  if(info!=0 || sign_Cq==-1 )
  {
    prob = -1.0e10;
    printf("lndet_ICq!\n");
    return prob;
  }
  
  prob -= 0.5*(lndet_C + lndet_ICq);
  
  return prob;
}

void set_covar_mat(double *theta)
{
  double t1, t2, nerr, syserr;
  double taud, sigma;
  int i, j;

  //alpha = 1.0;
  taud = exp(theta[1]);
  sigma = exp(theta[0]) * sqrt(taud/2.0);
  syserr = exp(theta[ntheta-1]);

// first con-con
  for(i=0; i<ncon_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<i; j++)
    {
      t2 = Tcon_data[j];
      Smat[i*nall_data+j] = Smat[j*nall_data+i] = exp( - fabs(t2-t1)/taud );
      Nmat[i*nall_data+j] = Nmat[j*nall_data+i] = 0.0;
    }
    nerr = Fcerrs_data[i];
    Nmat[i*nall_data+i] = (nerr*nerr + syserr*syserr)/(sigma*sigma);
    Smat[i*nall_data+i] = 1.0;
  }  

// then con-line and line-con
  for(i=0; i<ncon_data; i++)
  {
    t1 = Tcon_data[i];
    for(j=0; j<nline_data; j++)
    {
      t2 = Tline_data[j];
      Smat[ i*nall_data + (ncon_data+j)] = Smat[(j+ncon_data)*nall_data + i] =  Slc(t1, t2, theta);
      Nmat[ i*nall_data + (ncon_data+j)] = Nmat[(j+ncon_data)*nall_data + i] = 0.0;
    }
  }

/*// then line-con
  for(i=0; i<nline_data; i++)
  {
    t1 = Tline_data[i];
    for(j=0; j<ncon_data; j++)
    {
      t2 = Tcon_data[j];
      Smat[ (i+ncon_data)*n_data + j] = Slc(t2, t1, theta);
      Nmat[ (i+ncon_data)*n_data + j] = 0.0;
    }
  }  */

// then line-line
  for(i=0; i<nline_data; i++)
  {
    t1 = Tline_data[i];
    for(j=0; j<i; j++)
    {
      t2 = Tline_data[j];
      Smat[ (i+ncon_data)*nall_data + (ncon_data+j)] = Smat[ (j+ncon_data)*nall_data + (ncon_data+i)] = Sll(t1, t2, theta);
      Nmat[ (i+ncon_data)*nall_data + (ncon_data+j)] = Nmat[ (j+ncon_data)*nall_data + (ncon_data+i)] = 0.0;
    }
    nerr = Flerrs_data[i];
    Nmat[ (i+ncon_data)*nall_data + (ncon_data+i) ] = (nerr * nerr + syserr*syserr)/(sigma*sigma);
    Smat[ (i+ncon_data)*nall_data + (ncon_data+i) ] = Sll(t1, t1, theta);
  }

}

void set_covar_Umat(double *theta)
{
  double t1, t2, taud;
  int i, j;
 
  
  taud = exp(theta[1]);
  //sigma = exp(theta[0]) * sqrt(taud/2.0);

  for(i=0; i<ncon; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<ncon_data; j++)
    {
      t2 = Tcon_data[j];
      USmat[i*nall_data+j] = exp (- fabs(t1-t2) / taud );
    }

    for(j=0; j<nline_data; j++)
    {
      t2 = Tline_data[j];
      USmat[i*nall_data + (ncon_data+j)] = Slc(t1, t2, theta);
    }
  }

  for(i=0; i<nline; i++)
  {
    t1 = Tline[i];
    for(j=0; j<ncon_data; j++)
    {
      t2 = Tcon_data[j];
      USmat[ (i+ncon)*nall_data + j] = Slc(t2, t1, theta);
    }

    for(j=0; j<nline_data; j++)
    {
      t2 = Tline_data[j];
      USmat[ (i+ncon)*nall_data + (ncon_data+j)] = Sll(t1, t2, theta);
    }
  }
  return;
}

void set_covar_Amat(double *theta)
{
  double t1, t2, taud;
  int i, j;
 
  taud = exp(theta[1]);
  //sigma = exp(theta[0]) * sqrt(taud/2.0);

  for(i=0; i<ncon; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<ncon; j++)
    {
      t2 = Tcon[j];
      ASmat[i*nall+j] = exp (- fabs(t1-t2) / taud );
    }

    for(j=0; j<nline; j++)
    {
      t2 = Tline[j];
      ASmat[i*nall + (ncon+j)] = Slc(t1, t2, theta);
    }
  }

  for(i=0; i<nline; i++)
  {
    t1 = Tline[i];
    for(j=0; j<ncon; j++)
    {
      t2 = Tcon[j];
      ASmat[ (i+ncon)*nall + j] = Slc(t2, t1, theta);
    }

    for(j=0; j<nline; j++)
    {
      t2 = Tline[j];
      ASmat[ (i+ncon)*nall + (ncon+j)] = Sll(t1, t2, theta);
    }
  }
  return;
}


double Slc(double tcon, double tline, double *theta)
{
  double Dt, DT, w, tauk, fk, Sk, Stot, taud;
  int i;

  Dt = tline - tcon;
  taud = exp(theta[1]);
  w = theta[2];

  Stot = 0.0;
#ifdef TOPHAT
  for(i=0; i<nc; i++)
  {
    tauk = grid_tau[i]; //theta[3+nc+i];
    fk = exp(theta[3+i]);
    DT = Dt - tauk;
    if(DT < -w)
    {
      Sk = exp( (DT + w) / taud ) - exp( (DT - w)/taud);
    }
    else if(DT < w)
    {
      Sk = 2.0 - exp(- (DT + w) / taud) - exp( (DT - w)/taud);
    }
    else
    {
      Sk = exp( - (DT - w)/taud) - exp ( - (DT + w)/taud);
    }
    Stot += fk * Sk;
  }
  Stot *= (taud/2.0/w);
#elif defined JAVELINE 
  tauk = theta[4]; //theta[3+nc+i];
  fk = exp(theta[3]);
  DT = Dt - tauk;
  if(DT < -w)
  {
    Sk = exp( (DT + w) / taud ) - exp( (DT - w)/taud);
  }
  else if(DT < w)
  {
    Sk = 2.0 - exp(- (DT + w) / taud) - exp( (DT - w)/taud);
  }
  else
  {
    Sk = exp( - (DT - w)/taud) - exp ( - (DT + w)/taud);
  }
  Stot = fk * Sk;
  Stot *= (taud/2.0/w);
#else
  for(i=0; i<nc; i++)
  {
    tauk = grid_tau[i]; //theta[3+nc+i];
    fk = exp(theta[3+i]);
    DT = Dt - tauk;
    Sk = exp(-DT/taud) * erfc( -(DT/w - w/taud)/sqrt(2.0) ) 
        +exp( DT/taud) * erfc(  (DT/w + w/taud)/sqrt(2.0) );

    Stot += Sk * fk;
  }
  Stot *= 1.0/2.0 * exp(w*w/2.0/taud/taud);
#endif
  return Stot;
}

double Sll(double ti, double tj, double *theta)
{
  double Dt, DT, w, tauk, fk, taum, fm, Skm, Stot, taud;
  int k, m;

  Dt = ti - tj;
  taud = exp(theta[1]);
  w = theta[2];
 
  Stot = 0.0;
#ifdef TOPHAT 
  for(k=0; k<nc; k++)
  {
    tauk = grid_tau[k]; //theta[3+nc+k];
    fk = exp(theta[3+k]);
    for(m=0; m<nc; m++)
    {
      taum = grid_tau[m]; //theta[3+nc+m];
      fm = exp(theta[3+m]);

      DT = Dt - (tauk-taum);

      if(DT < -2.0*w)
      {
        Skm = exp((DT + 2.0*w)/taud) + exp((DT - 2.0*w)/taud) - 2.0*exp(DT/taud);
      }
      else if(DT < 0.0)
      {
        Skm = exp(-(DT + 2.0*w)/taud)  + exp((DT - 2.0*w)/taud)  - 2.0*exp(DT/taud) + 2.0*(DT+2.0*w)/taud;
      }
      else if(DT < 2.0*w)
      {
        Skm = exp(-(DT + 2.0*w)/taud) + exp((DT - 2.0*w)/taud)  - 2.0*exp(-DT/taud) - 2.0*(DT-2.0*w)/taud;
      }
      else
      {
        Skm = exp( -(DT + 2.0*w)/taud) + exp(-(DT - 2.0*w)/taud) - 2.0*exp(-DT/taud);
      }
      Stot += fk*fm*Skm;
    }
  } 
  Stot *= (taud*taud/4.0/w/w);

#elif defined JAVELINE
  for(k=0; k<1; k++)
  {
    tauk = theta[3+k+1]; //theta[3+nc+k];
    fk = exp(theta[3+k]);
    for(m=0; m<1; m++)
    {
      taum = theta[3+m+1]; //theta[3+nc+m];
      fm = exp(theta[3+m]);

      DT = Dt - (tauk-taum);

      if(DT < -2.0*w)
      {
        Skm = exp((DT + 2.0*w)/taud) + exp((DT - 2.0*w)/taud) - 2.0*exp(DT/taud);
      }
      else if(DT < 0.0)
      {
        Skm = exp(-(DT + 2.0*w)/taud)  + exp((DT - 2.0*w)/taud)  - 2.0*exp(DT/taud) + 2.0*(DT+2.0*w)/taud;
      }
      else if(DT < 2.0*w)
      {
        Skm = exp(-(DT + 2.0*w)/taud) + exp((DT - 2.0*w)/taud)  - 2.0*exp(-DT/taud) - 2.0*(DT-2.0*w)/taud;
      }
      else
      {
        Skm = exp( -(DT + 2.0*w)/taud) + exp(-(DT - 2.0*w)/taud) - 2.0*exp(-DT/taud);
      }
      Stot += fk*fm*Skm;
    }
  } 
  Stot *= (taud*taud/4.0/w/w);
#else  
  for(k=0; k<nc; k++)
  {
    tauk = grid_tau[k]; //theta[3+nc+k];
    fk = exp(theta[3+k]);
    for(m=0; m<nc; m++)
    {
      taum = grid_tau[m]; //theta[3+nc+m];
      fm = exp(theta[3+m]);
      DT = Dt - (tauk-taum);
      Skm = ( exp(-DT/taud) * erfc( -DT/2.0/w + w/taud )
            +  exp( DT/taud) * erfc(  DT/2.0/w + w/taud ) );
      Stot += Skm * fm*fk;
    }
  }
  Stot *= 1.0/2.0 * exp( w*w/taud/taud);
#endif
  return Stot;
} 

int reconstruct_conline()
{
  FILE *frec;
  double *Larr, *Larr_rec, *ave, *yave, *yave_rec, *ysub, *ybuf, *yrec, *yrec_err, *Cq, *ICq;
  double *Tmp1, *Tmp2, *Tmp3, *Tmp4;
  double sigma, taud;
  int i, nq, info;
  
  nq = 2*(flag_detrend + 1);

  Larr = array_malloc(nall_data * nq);
  Larr_rec = array_malloc(nall * nq);
  ybuf = array_malloc(nall_data);
  ysub = array_malloc(nall_data);
  yave = array_malloc(nall_data);
  yave_rec = array_malloc(nall);
  yrec = array_malloc(nall);
  yrec_err = array_malloc(nall);
  ave = array_malloc(nq);
  Cq = array_malloc(nq*nq);
  ICq = array_malloc(nq*nq);
  

  Tmp1 = array_malloc(max(nall,nall_data)*max(nall,nall_data));
  Tmp2 = array_malloc(max(nall,nall_data)*max(nall,nall_data));
  Tmp3 = array_malloc(max(nall,nall_data)*max(nall,nall_data));
  Tmp4 = array_malloc(max(nall,nall_data)*max(nall,nall_data));
   
  taud = exp(theta_best[1]);
  sigma = exp(theta_best[0]) * sqrt(taud/2.0);
  set_covar_mat(theta_best);
  set_covar_Umat(theta_best);
  set_covar_Amat(theta_best);

  for(i=0;i<nall_data*nall_data; i++)
  {
    Cmat[i] = Smat[i] + Nmat[i];
  }
  memcpy(ICmat, Cmat, nall_data*nall_data*sizeof(double));
  inverse_mat(ICmat, nall_data, &info);

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

  for(i=0; i<ncon; i++)
  {
    if(flag_detrend==0)
    {
      Larr_rec[i*nq+0] = 1.0;
      Larr_rec[i*nq+1] = 0.0;
    }
    else
    {
      Larr_rec[i*nq + 0] = 1.0;
      Larr_rec[i*nq + 1] = Tcon[i];
      Larr_rec[i*nq + 2] = 0.0;
      Larr_rec[i*nq + 3] = 0.0;
    }
  }
  for(i=0; i<nline; i++)
  {
    if(flag_detrend==0)
    {
      Larr_rec[(i+ncon)*nq+0] = 0.0;
      Larr_rec[(i+ncon)*nq+1] = 1.0;
    }
    else
    {
      Larr_rec[(i+ncon)*nq + 0] = 0.0;
      Larr_rec[(i+ncon)*nq + 1] = 0.0;
      Larr_rec[(i+ncon)*nq + 2] = 1.0;
      Larr_rec[(i+ncon)*nq + 3] = Tline[i];
    }
  }

  multiply_mat_MN(ICmat, Larr, Tmat1, nall_data, nq, nall_data);
  multiply_mat_MN_transposeA(Larr, Tmat1, ICq, nq, nq, nall_data);
  memcpy(Cq, ICq, nq*nq*sizeof(double));
  inverse_mat(Cq, nq, &info);

  multiply_mat_MN_transposeA(Larr, ICmat, Tmat1, nq, nall_data, nall_data);
  multiply_mat_MN(Cq, Tmat1, Tmat2, nq, nall_data, nq);
  multiply_mat_MN(Tmat2, Fall_data, ave, nq, 1, nall_data);

  multiply_mat_MN(Larr, ave, yave, nall_data, 1, nq);
  for(i=0; i<nall_data; i++)ysub[i] = Fall_data[i] - yave[i];
  multiply_matvec(ICmat, ysub, nall_data, ybuf);

  multiply_matvec_MN(USmat, nall, nall_data, ybuf, yrec);
  multiply_mat_MN(Larr_rec, ave, yave_rec, nall, 1, nq);

/* store the reconstructed, mean-substracted contionum */  
  memcpy(Fcon, yrec, ncon*sizeof(double));

  for(i=0; i<nall; i++)yrec[i] += yave_rec[i];

  // get errors
  
  multiply_mat_MN(USmat, ICmat, Tmp1, nall, nall_data, nall_data);
  multiply_mat_MN_transposeB(Tmp1, USmat, Tmp2, nall, nall, nall_data);
  
  multiply_mat_MN(Tmp1, Larr, Tmp3, nall, nq, nall_data);
  for(i=0; i<nall*nq; i++)Tmp3[i] -= Larr_rec[i];
  multiply_mat_MN(Tmp3, Cq, Tmp1, nall, nq, nq);
  multiply_mat_MN_transposeB(Tmp1, Tmp3, Tmp4, nall, nall, nq);

  for(i=0; i<nall; i++)
  {
    yrec_err[i] = sigma * sqrt(ASmat[i*nall+i] - Tmp2[i*nall+i] + Tmp4[i*nall+i]);
  }

/* store the error */
  for(i=0; i<ncon; i++)
  {
    Fcerrs[i] = sigma * sqrt(ASmat[i*nall+i] - Tmp2[i*nall+i]);
  }

  frec = fopen("data/sall_con.txt", "w");
  for(i=0; i<ncon; i++)
  {
    fprintf(frec, "%f %e %e\n", Tcon[i], yrec[i]*scale_con, yrec_err[i]*scale_con);
  }
  fclose(frec);
  frec = fopen("data/sall_line.txt", "w");
  for(i=ncon; i<nall; i++)
  {
    fprintf(frec, "%f %e %e\n", Tline[i-ncon], yrec[i]*scale_line, yrec_err[i]*scale_line);
  }
  fclose(frec);
  free(Tmp1);
  free(Tmp2);
  free(Tmp3);
  free(Tmp4);
  return 0;
}
