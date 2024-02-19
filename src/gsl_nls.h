#ifndef GSLNLS_H
#define GSLNLS_H

#include <R.h>
#include <Rinternals.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
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
    SEXP swts;                         // weights
    SEXP lupars;                       // lower/upper parameter bounds
    SEXP control_int;                  // integer control paramaters
    SEXP control_dbl;                  // double control paramaters
    SEXP has_start;                    // flags if starting values fixed
    gsl_multifit_nlinear_workspace *w; // workspace
    gsl_vector *wts;                   // weights vector
    gsl_qrng *q;                       // qrng workspace
    gsl_matrix *mx;                    // multistart initial value matrix
    gsl_vector *mp;                    // multistart parameter placeholder vector
    gsl_vector *mpopt;                 // multistart optimal parameter vector
    gsl_vector *diag;                  // diagonal of scaling matrix D_k
    gsl_matrix *lu;                    // lower/upper parameter bounds
} pdata;

typedef struct
{
    R_len_t p;     // number of parameters
    R_len_t n;     // number of observations
    double chisq;  // current ssr
    SEXP f;        // f language call
    SEXP df;       // df (Jacobian) language call
    SEXP fvv;      // fvv (acceleration) language call
    SEXP rho;      // environment in which to evaluate f
    SEXP y;        // observation vector y
    SEXP start;    // start parameter values
    SEXP partrace; // parameter trace
    SEXP ssrtrace; // ssr trace
    Rboolean warn; // print R warnings
    int startisnum;  // are start values numeric?
} fdata;

typedef struct
{
    R_len_t n;        // number of quasirandom samples
    R_len_t p;        // number of concentration iters
    R_len_t q;        // number of points in reduced sample
    R_len_t s;        // start local search if ntix(j) >= s
    R_len_t niter;    // number of local search iters
    R_len_t max;      // maximum number of major iters
    R_len_t minsp;    // minimum number of stationary points
    double r;         // stopping criterion nwsp >= r * nsp
    double tol;       // start local search if mssr < (1 + tol) * mssropt
    int *ntix;        // persistent index counter across major iters
    double *qmp;      // quasirandom samples
    double *startptr; // pointer to start value matrix
    int *mssr_order;  // ordered ssr's
    int mstop;        // stop status
    R_len_t mstarts;  // current number of major iteration
    R_len_t nsp;      // current number of stationary points
    R_len_t nwsp;     // current number of worse stationary points
    double mssropt;   // current optimal ssr
    double ssrconv;   // achieved convergence tolerance
    Rboolean all_start; // flag if all fixed lower/upper sampling limits
    int *has_start;     // flags fixed lower/upper sampling limits
    int *luchange;    // force change in lower/upper sampling limits
    double rejectscl; // scaling factor parameter rejection limits
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
    gsl_vector *workp;   /* workspace, length p */
    gsl_vector *workn;   /* workspace, length n */
    void *trs_state;     /* workspace for trust region subproblem */
    void *solver_state;  /* workspace for linear least squares solver */
    double avratio;      /* current |a| / |v| */
    /* tunable parameters */
    gsl_multifit_nlinear_parameters params;
} trust_state_t;

// double wgt(double x, double *c, int i);

int trust_iterate_lu(void *vstate,
                     const gsl_vector *swts,
                     gsl_multifit_nlinear_fdf *fdf,
                     gsl_vector *x,
                     gsl_vector *f,
                     gsl_matrix *J,
                     gsl_vector *g,
                     gsl_vector *dx,
                     const gsl_matrix *lu);

SEXP C_nls(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP lupars, SEXP control_int, SEXP control_dbl, SEXP has_start);

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
    SEXP swts;                           // weights
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

SEXP C_nls_large(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl);

SEXP C_nls_large_internal(void *data);

/* nls_test.c */
void F77_NAME(p00_f)(int *nprob, int *m, int *n, void *x, void *f);
void F77_NAME(p00_j)(int *nprob, int *m, int *n, void *x, void *fjac);
void F77_NAME(p00_start)(int *nprob, int *n, void *x);
void F77_NAME(p00_sol)(int *nprob, int *m, int *n, int *known, void *x);

SEXP C_nls_test_f(SEXP id, SEXP p, SEXP n, SEXP x);
SEXP C_nls_test_j(SEXP id, SEXP p, SEXP n, SEXP x);
SEXP C_nls_test_start_sol(SEXP id, SEXP p, SEXP n);

#endif
