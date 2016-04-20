#include "allvars.h"

// file names
char fname_param[100];
char fname_con[100], fname_line[100];

// transfer function parameters
double tau_lim_low, tau_lim_up;
int nc;

// mcmc
int nmcmc, nbuilt;

// light curve data
const int ndata_max=500;
double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;

// error exit
char str_error_exit[100];