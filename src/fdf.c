#include "gsl_nls.h"

/*
    these internal functions are adapted from multifit_nlinear/fdf.c
    to using a general symmetric weight matrix in order to perform
    generalized instead of weighted least squares
*/

/*
gsl_multifit_nlinear_winit_LD()
  Initialize generalized least squares solver and apply
  weighting transform using LDDL' decomposition of W if given:

  Inputs:
        x    - model parameters
        DDw  - diagonal D^2 in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDL' decomposition of W
               set to NULL for unweighted fit
        fdf  - callback function
        w    - workspace

*/
int gsl_multifit_nlinear_winit_LD(const gsl_vector *x,
                                  const gsl_vector *DDw,
                                  const gsl_matrix *Lw,
                                  gsl_multifit_nlinear_fdf *fdf,
                                  gsl_multifit_nlinear_workspace *w)
{
  const size_t n = w->f->size;

  if (n != fdf->n)
  {
    GSL_ERROR("function size does not match workspace", GSL_EBADLEN);
  }
  else if (w->x->size != x->size)
  {
    GSL_ERROR("vector length does not match workspace", GSL_EBADLEN);
  }
  else if (DDw != NULL && n != DDw->size)
  {
    GSL_ERROR("weight vector length does not match workspace", GSL_EBADLEN);
  }
  else
  {
    size_t i;

    /* initialize counters for function and Jacobian evaluations */
    fdf->nevalf = 0;
    fdf->nevaldf = 0;
    fdf->nevalfvv = 0;

    w->fdf = fdf;
    gsl_vector_memcpy(w->x, x);
    w->niter = 0;

    if (DDw)
    {
      w->sqrt_wts = w->sqrt_wts_work;

      for (i = 0; i < n; ++i)
      {
        double wi = gsl_vector_get(DDw, i);
        gsl_vector_set(w->sqrt_wts, i, sqrt(wi));
      }
    }
    else
    {
      w->sqrt_wts = NULL;
    }

    return trust_init_LD(w->state, w->sqrt_wts, Lw, w->fdf,
                         w->x, w->f, w->J, w->g);
  }

  GSL_ERROR("workspace initialization failed", GSL_EOF);

}

/*
gsl_multifit_nlinear_eval_f_LD()
  Compute residual vector y with user callback function, and apply
  weighting transform using LDDL' decomposition of W if given:

y~ = sqrt(W) y

Inputs: fdf  - callback function
        x    - model parameters
        Dw   - diagonal D in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDDL' decomposition of W
               set to NULL for unweighted fit
        y    - (output) (weighted) residual vector
*/

int gsl_multifit_nlinear_eval_f_LD(gsl_multifit_nlinear_fdf *fdf,
                                   const gsl_vector *x,
                                   const gsl_vector *Dw,
                                   const gsl_matrix *Lw,
                                   gsl_vector *y)
{
  int s = ((*((fdf)->f))(x, fdf->params, y));

  ++(fdf->nevalf);

  /* y <- sqrt(W) y */
  if (Lw)
  {
    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasUnit, Lw, y);
    if (Dw)
      gsl_vector_mul(y, Dw);
  }

  return s;
}

/*
gsl_multifit_nlinear_eval_df_LD()
  Compute Jacobian matrix J with user callback function, and apply
  weighting transform using LDDL' decomposition of W if given:

J~ = sqrt(W) J

Inputs: x      - model parameters
        f      - residual vector f(x)
        Dw     - diagonal D in LDDL' decomposition of W
        Lw     - unit lower triangular matrix L in LDDL' decomposition of W
                 set to NULL for unweighted fit
        h      - finite difference step size
        fdtype - finite difference method
        fdf    - callback function
        df     - (output) (weighted) Jacobian matrix
                 df = sqrt(W) df where df is unweighted Jacobian
        work   - workspace for finite difference, size n
*/

int gsl_multifit_nlinear_eval_df_LD(const gsl_vector *x,
                                    const gsl_vector *f,
                                    const gsl_vector *Dw,
                                    const gsl_matrix *Lw,
                                    const double h,
                                    const gsl_multifit_nlinear_fdtype fdtype,
                                    gsl_multifit_nlinear_fdf *fdf,
                                    gsl_matrix *df,
                                    gsl_vector *work)
{
  int status;

  if (fdf->df)
  {
    /* call user-supplied function */
    status = ((*((fdf)->df))(x, fdf->params, df));
    ++(fdf->nevaldf);

    /* J <- sqrt(W) J */
    if (Lw)
    {
      gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, Lw, df);
      if (Dw)
      {
        const size_t n = Dw->size;
        size_t i;
        for (i = 0; i < n; ++i)
        {
          double swi = gsl_vector_get(Dw, i);
          gsl_vector_view v = gsl_matrix_row(df, i);
          gsl_vector_scale(&v.vector, swi);
        }
      }
    }
  }
  else
  {
    /* use finite difference Jacobian approximation */
    status = gsl_multifit_nlinear_df_LD(h, fdtype, x, Dw, Lw, fdf, f, df, work);
  }

  return status;
}

/*
gsl_multifit_nlinear_eval_fvv_LD()
Compute second direction derivative vector yvv with user
callback function, and apply weighting transform using LDDL'
decomposition of W if given:

yvv~ = sqrt(W) yvv

Inputs: h    - step size for finite difference, if needed
        x    - model parameters, size p
        v    - unscaled geodesic velocity vector, size p
        f    - residual vector f(x), size n
        J    - Jacobian matrix J(x), n-by-p
        Dw   - diagonal D in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDDL' decomposition of W
               set to NULL for unweighted fit
        fdf  - callback function
        yvv  - (output) (weighted) second directional derivative vector
        work - workspace, size p
*/

int gsl_multifit_nlinear_eval_fvv_LD(const double h,
                                     const gsl_vector *x,
                                     const gsl_vector *v,
                                     const gsl_vector *f,
                                     const gsl_matrix *J,
                                     const gsl_vector *Dw,
                                     const gsl_matrix *Lw,
                                     gsl_multifit_nlinear_fdf *fdf,
                                     gsl_vector *yvv, gsl_vector *work)
{
  int status;

  if (fdf->fvv != NULL)
  {
    /* call user-supplied function */
    status = ((*((fdf)->fvv))(x, v, fdf->params, yvv));
    ++(fdf->nevalfvv);

    /* yvv <- sqrt(W) yvv */
    if (Lw)
    {
      gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasUnit, Lw, yvv);
      if (Dw)
        gsl_vector_mul(yvv, Dw);
    }
  }
  else
  {
    /* use finite difference approximation */
    status = gsl_multifit_nlinear_fdfvv_LD(h, x, v, f, J, Dw, Lw, fdf, yvv, work);
  }

  return status;
}
