#ifndef GSLNLS_H
#define GSLNLS_H

#include <R.h>
#include <Rinternals.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_multilarge_nlinear.h>

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
    gsl_matrix *J; // jacobian matrix
} fdata;

/* helper functions */
int gsl_f(const gsl_vector *x, void *params, gsl_vector *f);
int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J);
int gsl_fvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv);
int gsl_df_large(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u, void *params, gsl_vector *v, gsl_matrix *JTJ);
void callback(const size_t iter, void *params, const gsl_multifit_nlinear_workspace *w);
int gsl_multifit_nlinear_driver2(const size_t maxiter,
                                 const double xtol,
                                 const double gtol,
                                 const double ftol,
                                 void (*callback)(const size_t iter, void *params,
                                                  const gsl_multifit_nlinear_workspace *w),
                                 void *callback_params,
                                 int *info,
                                 double *chisq,
                                 gsl_multifit_nlinear_workspace *w);
int gsl_multilarge_nlinear_driver2(const size_t maxiter,
                                   const double xtol,
                                   const double gtol,
                                   const double ftol,
                                   void (*callback)(const size_t iter, void *params,
                                                    const gsl_multilarge_nlinear_workspace *w),
                                   void *callback_params,
                                   int *info,
                                   gsl_multilarge_nlinear_workspace *w);

/* main functions */
SEXP C_nls(SEXP f, SEXP y, SEXP df, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl);
SEXP C_nls_large(SEXP f, SEXP y, SEXP df, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl);
#endif
