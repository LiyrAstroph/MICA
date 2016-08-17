#include <gsl/gsl_rng.h>

#include "allvars.h"

// file names
char fname_param[100];
char fname_con[100], fname_line[100], fname_results[100];

FILE *fp_results;

// transfer function parameters
double tau_lim_low, tau_lim_up;
int nc, nc_lim_low, nc_lim_up, nc_best;
double *grid_tau;
const int ntau=200;
double *TF_tau, *TF;

// mcmc
const int ntheta_max = 200;
int ntheta;
int nmcmc, nbuilt;
double *workspace;
int flag_detrend, flag_sim, flag_mcmc;
double *Cmat, *ICmat, *Smat, *Nmat, *USmat, *ASmat;
double *ICvmat, *Tmat1, *Tmat2, *Tmat3, *Tmat4;

double *cov_matrix;
double **theta_range, *theta_best, *theta_best_var;
double *theta_input, *sigma_input;
int *theta_fixed;

double *theta_best_con, *theta_best_var_con;

const gsl_rng_type * gsl_T;
const gsl_rng * gsl_r;


// light curve data
const int ndata_max=500;
int ncon_data, nline_data, nall_data;
double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;
double *Fall_data;
double scale_con, scale_line;
double len_con, len_line, cad_con, cad_line;

int ncon=500, nline=500, nall;
double *Tcon, *Fcon, *Fcerrs, *Tline, *Fline, *Flerrs;

// error exit
char str_error_exit[200];

// mathematic functions
int *workspace_ipiv;