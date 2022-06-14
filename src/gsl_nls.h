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
    SEXP control_int;                  // integer control paramaters
    SEXP control_dbl;                  // double control paramaters
    gsl_multifit_nlinear_workspace *w; // workspace
    gsl_vector *wts;                   // weights vector
    gsl_qrng *q;                       // qrng workspace
    gsl_vector *mf;                    // multistart f vector
    gsl_vector *mpi;                   // multistart parameter vector
    gsl_vector *mpopt;                 // multistart optimal parameter vector
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
} fdata;

int gsl_f(const gsl_vector *x, void *params, gsl_vector *f);

int gsl_evalf(SEXP par, fdata *params, gsl_vector *f);

int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J);

int gsl_fvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv);

void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);

int gsl_multifit_nlinear_driver2(const size_t maxiter,
                                 const double xtol,
                                 const double gtol,
                                 const double ftol,
                                 void (*callback)(const size_t iter, void *params,
                                                  const gsl_multifit_nlinear_workspace *w),
                                 void *callback_params,
                                 int *info,
                                 double *chisq0,
                                 double *chisq1,
                                 gsl_multifit_nlinear_workspace *w);

SEXP C_nls(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl);

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

int gsl_df_large(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u, void *params, gsl_vector *v, gsl_matrix *JTJ);

void callback_large(const size_t iter, void *params, const gsl_multilarge_nlinear_workspace *w);

int gsl_multilarge_nlinear_driver2(const size_t maxiter,
                                   const double xtol,
                                   const double gtol,
                                   const double ftol,
                                   void (*callback)(const size_t iter, void *params,
                                                    const gsl_multilarge_nlinear_workspace *w),
                                   void *callback_params,
                                   int *info,
                                   double *chisq0,
                                   double *chisq1,
                                   gsl_multilarge_nlinear_workspace *w);

SEXP C_nls_large(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl);

SEXP C_nls_large_internal(void *data);

#endif
