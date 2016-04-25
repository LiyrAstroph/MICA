#ifndef _ALLVARS_H
#define _ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <math.h>

#define PI ((M_PI))

#define max(a, b) ((a>=b?a:b))


// file names
extern char fname_param[100];
extern char fname_con[100], fname_line[100], fname_results[100];

extern FILE *fp_results;

// transfer function parameters
extern double tau_lim_low, tau_lim_up;
extern int nc;
extern double *grid_tau;
extern const int ntau;
extern double *TF_tau, *TF;
// mcmc
extern const int ntheta_max;
extern int ntheta;
extern int nmcmc, nbuilt;
extern double *workspace;
extern int flag_detrend;
extern double *Cmat, *ICmat, *Smat, *Nmat, *USmat, *ASmat;
extern double *ICvmat, *Tmat1, *Tmat2, *Tmat3, *Tmat4;

extern double *cov_matrix;
extern double **theta_range, *theta_best, *theta_best_var;
extern double *theta_input, *sigma_input;
extern int *theta_fixed;

extern double *theta_best_con, *theta_best_var_con;

extern const gsl_rng_type * gsl_T;
extern const gsl_rng * gsl_r;

// light curve data
extern const int ndata_max;
extern int ncon_data, nline_data, nall_data;
extern double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;
extern double *Fall_data;
extern double scale_con, scale_line;
extern double len_con, len_line, cad_con, cad_line;

extern int ncon, nline, nall;
extern double *Tcon, *Fcon, *Fcerrs, *Tline, *Fline, *Flerrs;

// error exit
extern char str_error_exit[100];

// mathematic functions
extern int *workspace_ipiv;
#endif
