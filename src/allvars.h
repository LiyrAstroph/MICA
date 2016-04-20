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

// light curve data
extern const int ndata_max;
extern double *Tcon_data, *Tline_data, *Fcon_data, *Fcerrs_data, *Fline_data, *Flerrs_data;

// error exit
extern char str_error_exit[100];
#endif
