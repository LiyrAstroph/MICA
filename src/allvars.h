#ifndef _ALLVARS_H
#define _ALLVARS_H

// file names
extern char fname_param[100];
extern char fname_con[100], fname_line[100];

// transfer function parameters
extern double tau_lim_low, tau_lim_up;
extern int nc;

// mcmc
extern int nmcmc, nbuilt;
extern double *workspace;
extern int flag_detrend;
extern double *Cmat, *ICmat, *Smat, *Nmat, *USmat;
extern double *ICvmat, *Tmat1, *Tmat2, *Tmat3, *Tmat4, *ASmat;

// light curve data
extern const int ndata_max;
extern int ncon_data, nline_data, nall_data;
extern double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;

extern int ncon, nline, nall;
extern double *Tcon, *Fcon, *Fcerrs, *Tline, *Fline, *Flerrs;

// error exit
extern char str_error_exit[100];

// mathematic functions
extern int *workspace_ipiv;
#endif
