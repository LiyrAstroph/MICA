#include <gsl/gsl_rng.h>

#include "allvars.h"

// file names
char fname_param[100];
char fname_con[100], fname_line[100];

// transfer function parameters
double tau_lim_low, tau_lim_up;
int nc;
double *grid_tau;
const int ntau=100;
double *TF_tau, *TF;

// mcmc
const int ntheta_max = 50;
int ntheta;
int nmcmc, nbuilt;
double *workspace;
int flag_detrend;
double *Cmat, *ICmat, *Smat, *Nmat, *USmat, *ASmat;
double *ICvmat, *Tmat1, *Tmat2, *Tmat3, *Tmat4;

double *cov_matrix;
double **theta_range, *theta_best, *theta_best_var;
double *theta_input, *sigma_input;
int *theta_fixed;

const gsl_rng_type * gsl_T;
const gsl_rng * gsl_r;


// light curve data
const int ndata_max=500;
int ncon_data, nline_data, nall_data;
double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;
double *Fall_data;
double scale_con, scale_line;

int ncon=100, nline=100, nall;
double *Tcon, *Fcon, *Fcerrs, *Tline, *Fline, *Flerrs;

// error exit
char str_error_exit[100];

// mathematic functions
int *workspace_ipiv;