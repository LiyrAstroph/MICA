

//run.c
void run();

//read_param.c
void read_param();
void read_data();
void read_input();

//init.c
void init();
void memory_alloc();
void memory_alloc_data();
void scale_light_curves();

// mcmc_con.c
void mcmc_con_init();
void mcmc_con_run();
void set_covar_Umat_con(double sigma, double tau, double alpha);
void set_covar_mat_con(double sigma, double tau, double alpha);
double probability_con(double *);
void reconstruct_con();

//mcmc_conline.c
void mcmc_conline_init();
void mcmc_conline_run();
double Slc(double tcon, double tline, double *theta);
double Sll(double ti, double tj, double *theta);
void set_covar_mat(double *theta);
void set_covar_Umat(double *theta);
void set_covar_Amat(double *theta);
double probability_conline(double *theta);
double probability_conline_aicc(double *theta);
double cal_aicc();
int reconstruct_conline();

// mcmc.c
int mcmc_sampling(char *fname_mcmc, double (* prob_fun)(double *) );

// mcmc_conline.c
void mcmc_conline_init();
void mcmc_conline_sampling();
void mcmc_conline_run();
void simulate_line();

//mcmc_stats.c
void get_cov_matrix(double *theta, int nstep, int ntheta);
void get_cov_matrix_diag(double *theta, int nstep, int ntheta);
void mcmc_stats(char *);
int par_fit(double *theta, int n, double percent, double *p);
int fitfunc(int m, int n, double *p, double *dy, double **devc, void *vars);

// transferfunc.o
void transfer_function(double *theta);

// lineconv.o
void line_convolution();

// sim.c
void set_covar_Simmat(double sigma, double tau, double alpha);
void simulate_con_init();
void simulate_con(double sigma, double tau, double alpha, double ave_con);

//test.c
void test();

// error.c
void error_exit(int);

// mathfunc.c
void inverse_mat(double *a, int n, int *info);
double det_mat(double *a, int n, int *info);
double lndet_mat(double *a, int n, int *info);
double lndet_mat2(double *a, int n, int *info, int *sign);
double lndet_mat3(double *a, int n, int *info, int *sign);
void display_mat(double *a, int m, int n);
void multiply_mat(double * a, double *b, double *c, int n);
void multiply_mat_transposeA(double * a, double *b, double *c, int n);
void multiply_mat_transposeB(double * a, double *b, double *c, int n);
void multiply_mat_MN(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeA(double * a, double *b, double *c, int m, int n, int k);
void multiply_mat_MN_transposeB(double * a, double *b, double *c, int m, int n, int k);
int multiply_mat_MN_inverseA(double * a, double *b, int m, int n);
void multiply_matvec(double *a, double *x, int n, double *y);
void multiply_matvec_transposeA(double *a, double *x, int n, double *y);
void multiply_matvec_MN(double * a, int m, int n, double *x, double *y);
void multiply_vec2mat(double * x, double * a, int n);
void eigen_sym_mat(double *a, int n, double *val, int *info);
void Chol_decomp_U(double *a, int n, int *info);
double ** matrix_malloc(int n1, int n2);
double * array_malloc(int n);