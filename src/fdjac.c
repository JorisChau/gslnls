#include "gsl_nls.h"

/*
    these internal functions are adapted from multifit_nlinear/fdjac.c
    to using a general symmetric weight matrix in order to perform
    generalized instead of weighted least squares
*/

/*
forward_jac_LD()
  Compute approximate Jacobian using forward differences

Inputs: h   - finite difference step size
        x   - parameter vector
        Dw  - diagonal D in LDDL' decomposition of W
        Lw  - unit lower triangular matrix L in LDDL' decomposition of W
        fdf - fdf struct
        f   - (input) vector of function values f_i(x)
        J   - (output) Jacobian matrix

Return: success or error
*/

static int
forward_jac_LD(const double h, const gsl_vector *x, const gsl_vector *Dw, const gsl_matrix *Lw,
                 gsl_multifit_nlinear_fdf *fdf, const gsl_vector *f, gsl_matrix *J)
{
  int status = 0;
  size_t i, j;
  double delta;

  for (j = 0; j < fdf->p; ++j)
  {
    double xj = gsl_vector_get(x, j);

    /* use column j of J as temporary storage for f(x + dx) */
    gsl_vector_view v = gsl_matrix_column(J, j);

    delta = h * fabs(xj);
    if (delta == 0.0)
      delta = h;

    /* perturb x_j to compute forward difference */
    gsl_vector_set((gsl_vector *)x, j, xj + delta);

    status += gsl_multifit_nlinear_eval_f_LD(fdf, x, Dw, Lw, &v.vector);
    if (status)
      return status;

    /* restore x_j */
    gsl_vector_set((gsl_vector *)x, j, xj);

    delta = 1.0 / delta;
    for (i = 0; i < fdf->n; ++i)
    {
      double fnext = gsl_vector_get(&v.vector, i);
      double fi = gsl_vector_get(f, i);

      gsl_matrix_set(J, i, j, (fnext - fi) * delta);
    }
  }

  return status;
}

/*
center_jac_LD()
  Compute approximate Jacobian using centered differences

Inputs: h    - finite difference step size
        x    - parameter vector
        Dw   - diagonal D in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDDL' decomposition of W
        fdf  - fdf struct
        J    - (output) Jacobian matrix
        work - additional workspace, size n

Return: success or error
*/

static int
center_jac_LD(const double h, const gsl_vector *x, const gsl_vector *Dw, const gsl_matrix *Lw,
                gsl_multifit_nlinear_fdf *fdf, gsl_matrix *J, gsl_vector *work)
{
  int status = 0;
  size_t i, j;
  double delta;

  for (j = 0; j < fdf->p; ++j)
  {
    double xj = gsl_vector_get(x, j);

    /* use column j of J as temporary storage for f(x + dx) */
    gsl_vector_view v = gsl_matrix_column(J, j);

    delta = h * fabs(xj);
    if (delta == 0.0)
      delta = h;

    /* perturb x_j to compute forward difference, f(x + 1/2 delta e_j) */
    gsl_vector_set((gsl_vector *)x, j, xj + 0.5 * delta);

    status += gsl_multifit_nlinear_eval_f_LD(fdf, x, Dw, Lw, &v.vector);
    if (status)
      return status;

    /* perturb x_j to compute backward difference, f(x - 1/2 delta e_j) */
    gsl_vector_set((gsl_vector *)x, j, xj - 0.5 * delta);

    status += gsl_multifit_nlinear_eval_f_LD(fdf, x, Dw, Lw, work);
    if (status)
      return status;

    /* restore x_j */
    gsl_vector_set((gsl_vector *)x, j, xj);

    delta = 1.0 / delta;
    for (i = 0; i < fdf->n; ++i)
    {
      double fnext = gsl_vector_get(&v.vector, i);
      double fprev = gsl_vector_get(work, i);

      gsl_matrix_set(J, i, j, (fnext - fprev) * delta);
    }
  }

  return status;
}

/*
gsl_multifit_nlinear_df_LD()
  Compute approximate Jacobian using finite differences

Inputs: h      - finite difference step size
        fdtype - finite difference method
        x      - parameter vector
        Dw     - diagonal D in LDDL' decomposition of W
        Lw     - unit lower triangular matrix L in LDDL' decomposition of W
        fdf    - fdf
        f      - (input) function values f_i(x)
        J      - (output) approximate (weighted) Jacobian matrix, sqrt(W) * J
        work   - additional workspace for centered differences, size n

Return: success or error
*/

int gsl_multifit_nlinear_df_LD(const double h, const gsl_multifit_nlinear_fdtype fdtype,
                                 const gsl_vector *x, const gsl_vector *Dw,
                                 const gsl_matrix *Lw, gsl_multifit_nlinear_fdf *fdf,
                                 const gsl_vector *f, gsl_matrix *J, gsl_vector *work)
{
  int status;

  if (fdtype == GSL_MULTIFIT_NLINEAR_FWDIFF)
  {
    status = forward_jac_LD(h, x, Dw, Lw, fdf, f, J);
  }
  else if (fdtype == GSL_MULTIFIT_NLINEAR_CTRDIFF)
  {
    status = center_jac_LD(h, x, Dw, Lw, fdf, J, work);
  }
  else
  {
    GSL_ERROR("invalid specified fdtype", GSL_EINVAL);
  }

  return status;
}
