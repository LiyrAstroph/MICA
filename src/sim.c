#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_interp.h>

#include "allvars.h"
#include "proto.h" 

#define SNR (120.0)
int n_con_sim=301, n_line_sim=301;
double sim_error=10.0*1.0/SNR, sim_rel_error=1.0/SNR;

double *PQmat, *PSmat, *Peigens, *Peigens_mat, *Peigens_vecs;

gsl_rng_type const * gsl_T_sim;
gsl_rng * gsl_r_sim;
gsl_interp_accel *gsl_acc, *gsl_acc_error;
gsl_interp  *gsl_linear, *gsl_linear_error;

void simulate_con(double sigma, double tau, double alpha, double ave_con)
{
  FILE *fscon_out;

  double *Py, *Pybuf, *Prandvec;
  int i,j, info;

  Py = array_malloc(n_con_sim);
  Pybuf = array_malloc(n_con_sim);
  Prandvec = array_malloc(n_con_sim);

  simulate_con_init();


  set_covar_Simmat(sigma, tau, alpha);
  memcpy(PQmat, PSmat, n_con_sim*n_con_sim*sizeof(double));
  
/*  eigen_sym_mat(PQmat, n_con_sim, Peigens, &info);
  memcpy(Peigens_vecs, PQmat, n_con_sim*n_con_sim*sizeof(double));

  for(i=0; i<n_con_sim; i++)
  {
    for(j=0; j<n_con_sim; j++)
      Peigens_mat[i*n_con_sim+j] = 0.0;
    Peigens_mat[i*n_con_sim+i] = sqrt(Peigens[i]);
  }*/

  for(i=0;i<n_con_sim; i++)
  {
    Prandvec[i]  = gsl_ran_gaussian(gsl_r_sim, 1.0);
  }
  
  Chol_decomp_U(PQmat, n_con_sim, &info);
  multiply_matvec(PQmat, Prandvec, n_con_sim, Py);

//   multiply_matvec(Peigens_mat, Prandvec, n_con_sim, Pybuf);
//   multiply_matvec(Peigens_vecs, Pybuf, n_con_sim, Py);

  for(i=0;i<n_con_sim;i++)
  {
    Fcon[i] = Py[i] + ave_con;
    Fcerrs[i] = sim_error;
  }
  
  fscon_out = fopen("data/sim_con.txt", "w");

  for(i=0;i<n_con_sim;i++)
  {
    if(Tcon[i]>150.0)
      fprintf(fscon_out, "%f\t%f\t%f\n", Tcon[i], Fcon[i], Fcerrs[i]);
  }

  gsl_interp_init(gsl_linear, Tcon, Fcon, n_con_sim);
  gsl_interp_init(gsl_linear_error, Tcon, Fcerrs, n_con_sim);

  fclose(fscon_out);
}

void simulate_con_init()
{
  int i;
  double Tcon_min, Tcon_max, Tline_min, Tline_max, dTs;

  PSmat = array_malloc(n_con_sim*n_con_sim);
  PQmat = array_malloc(n_con_sim*n_con_sim);
  Peigens = array_malloc(n_con_sim);
  Peigens_vecs = array_malloc(n_con_sim*n_con_sim);
  Peigens_mat = array_malloc(n_con_sim*n_con_sim); 

  Tcon = array_malloc(n_con_sim);
  Fcon = array_malloc(n_con_sim);
  Fcerrs = array_malloc(n_con_sim);

  Tcon_min = 0.0;
  Tcon_max = 300.0;
  dTs = (Tcon_max - Tcon_min) / (n_con_sim-1.0);
  for(i=0; i<n_con_sim; i++)
  {
    Tcon[i] = dTs * i + Tcon_min;
  }

  Tline = array_malloc(n_line_sim);
  Fline = array_malloc(n_line_sim);
  Flerrs = array_malloc(n_line_sim);

  Tline_min = 0.0;
  Tline_max = 300.0;
  dTs = (Tline_max - Tline_min) / (n_line_sim-1.0);
  for(i=0; i<n_line_sim; i++)
  {
    Tline[i] = dTs * i + Tline_min;
  }

  gsl_acc = gsl_interp_accel_alloc();
  gsl_acc_error = gsl_interp_accel_alloc();

  gsl_linear = gsl_interp_alloc(gsl_interp_linear, n_con_sim);
  gsl_linear_error = gsl_interp_alloc(gsl_interp_linear, n_con_sim);

  gsl_T_sim = gsl_rng_default;
  gsl_r_sim = gsl_rng_alloc (gsl_T_sim);
  //gsl_rng_set(gsl_r_sim, 12);  // time(NULL)
  gsl_rng_set(gsl_r_sim, time(NULL)); 

  return;
}

void set_covar_Simmat(double sigma, double tau, double alpha)
{
  double t1, t2, nerr;
  int i, j;
 
  for(i=0; i<n_con_sim; i++)
  {
    t1 = Tcon[i];
    for(j=0; j<=i; j++)
    {
      t2 = Tcon[j];
      PSmat[i*n_con_sim+j] = sigma*sigma* exp (- pow (fabs(t1-t2) / tau, alpha));
      PSmat[j*n_con_sim+i] = PSmat[i*n_con_sim+j];
    }
    nerr = sim_error;
    PSmat[i*n_con_sim+i] += nerr * nerr;
  }
  return;
}


void simulate_line()
{
  int i, j;
  double flux, err, tline, tcon, taup, dt, norm, fcon, fcon_err;
  double *tf, *tau, dtau, sig;
  int type=1;

  tf = array_malloc(ntau);
  tau = array_malloc(ntau);
  dtau = 30.0/(ntau-1.0);
  for(i=0; i<ntau; i++)
  {
    tau[i] =  dtau * i;
    tf[i] = 0.0;
  }
// top-hat transfer function
  FILE *ftf;
  ftf=fopen("data/sim_tf.txt", "w");
  norm = 0.0;
  for(i=0; i<ntau; i++)
  {
    if(type==0)  // two top-hats
    {
      if(tau[i]<=5.0)
        tf[i] = 0.0;
      else if(tau[i]<10.0)
        tf[i] = 1.0;
      else if(tau[i]<15.0)
        tf[i] = 0.0;
      else if(tau[i]<20.0)
        tf[i] = 2.0;
      else
        tf[i] = 0.0;

      norm += tf[i] * dtau;
    }
    if(type==1) // two Gaussians
    {
      sig = 2.5;
      tf[i] = exp(-0.5*pow(tau[i]-10.0, 2.0)/sig/sig) + exp(-0.5*pow(tau[i]-17.0, 2.0)/sig/sig);
      norm += tf[i] * dtau;
    }
  }
  for(i=0; i<ntau; i++)
  {
    tf[i] /= norm;
    fprintf(ftf, "%f %f\n", tau[i], tf[i]);
  }
  fclose(ftf);

  FILE *fp;
  fp = fopen("data/sim_line.txt", "w");

  for(i=0; i<n_line_sim; i++)
  {
    flux = 0.0;
    err = 0.0;
    tline = Tline[i]; //4540.0 + 70.0/(nline_data-1.0) * i;
    for(j=0; j<ntau; j++)
    {
      taup = tau[j];
      tcon = tline - taup;
      if(tcon >=Tcon[0] && tcon <= Tcon[n_con_sim-1])
      {
        fcon = gsl_interp_eval(gsl_linear, Tcon, Fcon, tcon, gsl_acc);
        fcon_err = gsl_interp_eval(gsl_linear_error, Tcon, Fcerrs, tcon, gsl_acc_error);
        flux += fcon * tf[j];
        err += pow(fcon_err *tf[j] * dtau, 2);
      }
    }
    flux *= dtau;
    err = sqrt(err);
    Fline[i] = flux;
    Flerrs[i] = sim_error;
  }

  
  for(i=0; i<n_line_sim; i++)
  {
    if(Tline[i]>170.0)
      fprintf(fp, "%f %f %f\n", Tline[i], Fline[i] + Flerrs[i]*gsl_ran_gaussian(gsl_r_sim, 1.0), Flerrs[i]);
  }
  fclose(fp);
}