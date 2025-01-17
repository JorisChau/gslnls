#ifndef GSLNLS_H
#define GSLNLS_H

#include <R.h>
#include <Rinternals.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>

/* nls.c */
typedef struct
{
    SEXP fn;                           // model function
    SEXP y;                            // response vector
    SEXP jac;                          // jacobian function
    SEXP fvv;                          // fvv function
    SEXP env;                          // function environment
    SEXP start;                        // start parameter values
    SEXP swts;                         // (cholesky) square root of weights
    SEXP lupars;                       // lower/upper parameter bounds
    SEXP control_int;                  // integer control paramaters
    SEXP control_dbl;                  // double control paramaters
    SEXP has_start;                    // flags if starting values fixed
    SEXP loss_config;                  // loss function configuration
    gsl_multifit_nlinear_workspace *w; // workspace
    gsl_vector *wts;                   // weights n-length vector 
    gsl_matrix *Lw;                    // unit lower triangular matrix L in LDL' 
                                       // decomposition of (n-by-n) matrix W
    gsl_qrng *q;                       // qrng workspace
    gsl_matrix *mx;                    // multistart initial value matrix
    gsl_vector *mpopt;                 // multistart optimal parameter vector
    gsl_vector *mpopt1;                // multistart back-up parameter vector
    gsl_vector *diag;                  // diagonal of scaling matrix D_k
    gsl_matrix *lu;                    // lower/upper parameter bounds
    gsl_matrix *JTJ;                   // workspace hessian matrix Jt * J
    gsl_vector *workn;                 // workspace n-length vector
    gsl_vector *workp;                 // workspace p-length vector
    gsl_matrix *worknp;                // workspace pxn-dimensional matrix
    gsl_vector *psi;                   // IRLS psi vector
    gsl_vector *psip;                  // IRLS psi' vector
} pdata;

typedef struct
{
    R_len_t p;      // number of parameters
    R_len_t n;      // number of observations
    double chisq;   // current ssr
    SEXP f;         // f language call
    SEXP df;        // df (Jacobian) language call
    SEXP fvv;       // fvv (acceleration) language call
    SEXP rho;       // environment in which to evaluate f
    SEXP y;         // observation vector y
    SEXP start;     // start parameter values
    SEXP partrace;  // parameter trace
    SEXP ssrtrace;  // ssr trace
    Rboolean warn;  // print R warnings
    int startisnum; // are start values numeric?
} fdata;

/* store info found local optima */
typedef struct
{
    double *mpall;   // all found local optima
    double *mpradii; // trust region radii of all found local optima
    R_len_t mpcount; // current number of local optima
    R_len_t mpmax;   // maximum number of local optima
} mpoptinfo;

typedef struct
{
    R_len_t n;         // number of quasirandom samples
    R_len_t p;         // number of concentration iters
    R_len_t q;         // number of points in reduced sample
    R_len_t s;         // start local search if ntix(j) >= s
    R_len_t niter;     // number of local search iters
    R_len_t max;       // maximum number of major iters
    R_len_t minsp;     // minimum number of stationary points
    int wgt_i;         // loss function identifier
    double r;          // stopping criterion nwsp >= r * nsp
    double tol;        // start local search if mssr < (1 + tol) * mssropt
    double dtol;       // discard point if det(J^T * J) < dtol
    int *ntix;         // persistent index counter across major iters
    double *qmp;       // quasirandom samples
    double *start;     // start value workspace
    double *maxlims;   // store maximum parameter limits
    int *mssr_order;   // ordered ssr's
    int mstop;         // stop status
    R_len_t mstarts;   // current number of major iteration
    R_len_t nsp;       // current number of stationary points
    R_len_t nwsp;      // current number of worse stationary points
    double mssropt[2]; // current optimal ssr + back-up
    double ssrconv[2]; // current convergence tolerance + back-up
    int all_start;     // flag if any dynamic lower/upper sampling limits
    int *has_start;    // flags fixed lower/upper sampling limits
    int *luchange;     // force change in lower/upper sampling limits
    double rejectscl;  // scaling factor parameter rejection limits
    mpoptinfo *mpopt;  // store all found local optima (currently not)
} mdata;

/* from multifit_nlinear/trust.c */
typedef struct
{
    size_t n;            /* number of observations */
    size_t p;            /* number of parameters */
    double delta;        /* trust region radius */
    double mu;           /* LM parameter */
    long nu;             /* for updating LM parameter */
    gsl_vector *diag;    /* D = diag(J^T J) */
    gsl_vector *x_trial; /* trial parameter vector */
    gsl_vector *f_trial; /* trial function vector */
    gsl_vector *workp;   /* workspace p-length vector */
    gsl_vector *workn;   /* workspace n-length vector */
    void *trs_state;     /* workspace for trust region subproblem */
    void *solver_state;  /* workspace for linear least squares solver */
    double avratio;      /* current |a| / |v| */
    /* tunable parameters */
    gsl_multifit_nlinear_parameters params;
} trust_state_t;

/* from multifit_nlinear/lm.c */
typedef struct
{
    size_t n;          /* number of observations */
    size_t p;          /* number of parameters */
    gsl_vector *fvv;   /* D_v^2 f(x), size n */
    gsl_vector *vel;   /* geodesic velocity (standard LM step), size p */
    gsl_vector *acc;   /* geodesic acceleration, size p */
    gsl_vector *workp; /* workspace, length p */
    gsl_vector *workn; /* workspace, length n */
    int accel; /* use geodesic acceleration? */
    /* tunable parameters */
    gsl_multifit_nlinear_parameters params;
} lm_state_t;

int gsl_multifit_nlinear_driver2(const R_len_t maxiter,
                                 const double xtol,
                                 const double gtol,
                                 const double ftol,
                                 void (*callback)(const R_len_t iter, void *params,
                                                  const gsl_multifit_nlinear_workspace *w),
                                 void *callback_params,
                                 int *info,
                                 double *chisq0,
                                 double *chisq1,
                                 const gsl_matrix *lu,
                                 const gsl_matrix *Lw,
                                 gsl_multifit_nlinear_workspace *w);

void gsl_multistart_driver(pdata *pars,
                           mdata *mpars,
                           gsl_multifit_nlinear_fdf *fdff,
                           SEXP mssr,
                           const double xtol,
                           const double ftol,
                           const Rboolean use_weights,
                           const Rboolean verbose);

int gsl_multifit_nlinear_rho_driver(pdata *pars,
                                    gsl_multifit_nlinear_fdf *fdff,
                                    const int wgt_i,
                                    const int maxiter,
                                    const double xtol,
                                    const double gtol,
                                    const double ftol,
                                    void *callback_params,
                                    int *info,
                                    double *chisq0,
                                    double *chisq1,
                                    double *irls_scale,
                                    R_len_t *irls_iter,
                                    int *irls_status,
                                    const Rboolean verbose);

int gsl_multifit_nlinear_winit_LD(const gsl_vector *x,
                                    const gsl_vector *DDw,
                                    const gsl_matrix *Lw,
                                    gsl_multifit_nlinear_fdf *fdf,
                                    gsl_multifit_nlinear_workspace *w);

int gsl_multifit_nlinear_eval_f_LD(gsl_multifit_nlinear_fdf *fdf,
                                     const gsl_vector *x,
                                     const gsl_vector *Dw,
                                     const gsl_matrix *Lw,
                                     gsl_vector *y);

int gsl_multifit_nlinear_eval_df_LD(const gsl_vector *x,
                                      const gsl_vector *f,
                                      const gsl_vector *Dw,
                                      const gsl_matrix *Lw,
                                      const double h,
                                      const gsl_multifit_nlinear_fdtype fdtype,
                                      gsl_multifit_nlinear_fdf *fdf,
                                      gsl_matrix *df,
                                      gsl_vector *work);

int gsl_multifit_nlinear_eval_fvv_LD(const double h,
                                       const gsl_vector *x,
                                       const gsl_vector *v,
                                       const gsl_vector *f,
                                       const gsl_matrix *J,
                                       const gsl_vector *Dw,
                                       const gsl_matrix *Lw,
                                       gsl_multifit_nlinear_fdf *fdf,
                                       gsl_vector *yvv, gsl_vector *work);

int gsl_multifit_nlinear_fdfvv_LD(const double h, const gsl_vector *x, const gsl_vector *v,
                                    const gsl_vector *f, const gsl_matrix *J, const gsl_vector *Dw,
                                    const gsl_matrix *Lw, gsl_multifit_nlinear_fdf *fdf,
                                    gsl_vector *fvv, gsl_vector *work);

int gsl_multifit_nlinear_df_LD(const double h, const gsl_multifit_nlinear_fdtype fdtype,
                                 const gsl_vector *x, const gsl_vector *Dw, 
                                 const gsl_matrix *Lw, gsl_multifit_nlinear_fdf *fdf,
                                 const gsl_vector *f, gsl_matrix *J, gsl_vector *work);

int trust_iterate_lu_LD(void *vstate, const gsl_vector *Dw, const gsl_matrix *Lw,
                          gsl_multifit_nlinear_fdf *fdf, gsl_vector *x,
                          gsl_vector *f, gsl_matrix *J, gsl_vector *g,
                          gsl_vector *dx, const gsl_matrix *lu);

int trust_init_LD(void *vstate, const gsl_vector *Dw, const gsl_matrix *Lw,
                    gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
                    gsl_vector *f, gsl_matrix *J, gsl_vector *g);

double det_cholesky_jtj(gsl_matrix *J, gsl_matrix *JTJ);

double det_eval_jtj(const gsl_multifit_nlinear_parameters params,
                    const gsl_vector *Dw, const gsl_matrix *Lw,
                    gsl_multifit_nlinear_fdf *fdf,
                    const gsl_vector *x, gsl_vector *f, gsl_matrix *J,
                    gsl_matrix *JTJ, gsl_vector *workn);

int hat_values(gsl_matrix *J, gsl_matrix *JTJ, gsl_vector *h, gsl_matrix *worknp);

int cooks_d(gsl_vector *f, gsl_matrix *J, gsl_matrix *JTJ, gsl_vector *d, gsl_matrix *worknp);

double gsl_median(double *data, const R_len_t n);

double gsl_mad(double *data, const R_len_t n, double *work);

SEXP C_nls(SEXP fn,
           SEXP y,
           SEXP jac,
           SEXP fvv,
           SEXP env,
           SEXP start,
           SEXP weigts,
           SEXP lupars,
           SEXP control_int,
           SEXP control_dbl,
           SEXP has_start,
           SEXP loss_config);

SEXP C_nls_internal(void *data);

/* nls_large.c */
typedef struct
{
    SEXP fn;                             // model function
    SEXP y;                              // response vector
    SEXP jac;                            // jacobian function
    SEXP fvv;                            // fvv function
    SEXP env;                            // function environment
    SEXP start;                          // start parameter values
    SEXP weights;                        // weights
    SEXP control_int;                    // integer control paramaters
    SEXP control_dbl;                    // double control paramaters
    gsl_multilarge_nlinear_workspace *w; // workspace
    gsl_matrix *J;                       // jacobian matrix
    gsl_spmatrix *Jsp;                   // sparse jacobian matrix
} pdata_large;

typedef struct
{
    R_len_t p;         // number of parameters
    R_len_t n;         // number of observations
    double chisq;      // current ssr
    SEXP f;            // f language call
    SEXP df;           // df (Jacobian) language call
    SEXP fvv;          // fvv (acceleration) language call
    SEXP rho;          // environment in which to evaluate f
    SEXP y;            // observation vector y
    SEXP start;        // start parameter values
    SEXP partrace;     // parameter trace
    SEXP ssrtrace;     // ssr trace
    int matclass;      // jacobian matrix class
    gsl_matrix *J;     // jacobian matrix
    gsl_spmatrix *Jsp; // sparse jacobian matrix
} fdata_large;

SEXP C_nls_large(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP weights, SEXP control_int, SEXP control_dbl);

SEXP C_nls_large_internal(void *data);

int gsl_multilarge_nlinear_driver2(const R_len_t maxiter,
                                   const double xtol,
                                   const double gtol,
                                   const double ftol,
                                   void (*callback)(const R_len_t iter, void *params,
                                                    const gsl_multilarge_nlinear_workspace *w),
                                   void *callback_params,
                                   int *info,
                                   double *chisq0,
                                   double *chisq1,
                                   gsl_multilarge_nlinear_workspace *w);

/* nls_test.c */
void F77_NAME(p00_f)(int *nprob, int *m, int *n, void *x, void *f);
void F77_NAME(p00_j)(int *nprob, int *m, int *n, void *x, void *fjac);
void F77_NAME(p00_start)(int *nprob, int *n, void *x);
void F77_NAME(p00_sol)(int *nprob, int *m, int *n, int *known, void *x);

SEXP C_nls_test_f(SEXP id, SEXP p, SEXP n, SEXP x);
SEXP C_nls_test_j(SEXP id, SEXP p, SEXP n, SEXP x);
SEXP C_nls_test_start_sol(SEXP id, SEXP p, SEXP n);

#endif
