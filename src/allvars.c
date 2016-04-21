#include "allvars.h"

// file names
char fname_param[100];
char fname_con[100], fname_line[100];

// transfer function parameters
double tau_lim_low, tau_lim_up;
int nc;

// mcmc
int nmcmc, nbuilt;
double *workspace;
int flag_detrend;
double *Cmat, *ICmat, *Smat, *Nmat, *USmat;
double *ICvmat, *Tmat1, *Tmat2, *Tmat3, *Tmat4, *ASmat;



// light curve data
const int ndata_max=500;
int ncon_data, nline_data, nall_data;
double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;

int ncon=100, nline=100, nall;
double *Tcon, *Fcon, *Fcerrs, *Tline, *Fline, *Flerrs;

// error exit
char str_error_exit[100];

// mathematic functions
int *workspace_ipiv;