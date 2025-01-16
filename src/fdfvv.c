#include "gsl_nls.h"

/*
    these internal functions are adapted from multifit_nlinear/fdfvv.c
    to using a general symmetric weight matrix in order to perform
    generalized instead of weighted least squares
*/

/*
gsl_multifit_nlinear_fdfvv_LD()
  Compute approximate second directional derivative
using finite differences

See Eq. 19 of:

M. K. Transtrum, J. P. Sethna, Improvements to the Levenberg
Marquardt algorithm for nonlinear least-squares minimization,
arXiv:1201.5885, 2012.

Inputs: h    - step size for finite difference
        x    - parameter vector, size p
        v    - geodesic velocity, size p
        f    - function values f_i(x), size n
        J    - Jacobian matrix J(x), n-by-p
        Dw   - diagonal D in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDDL' decomposition of W
        fdf  - fdf
        fvv  - (output) approximate (weighted) second directional derivative
               vector, size n, sqrt(W) fvv
        work - workspace, size p

Return: success or error
*/

int
gsl_multifit_nlinear_fdfvv_LD(const double h, const gsl_vector *x, const gsl_vector *v,
                           const gsl_vector *f, const gsl_matrix *J, const gsl_vector *Dw,
                           const gsl_matrix *Lw, gsl_multifit_nlinear_fdf *fdf,
                           gsl_vector *fvv, gsl_vector *work)
{
  int status;
  const size_t n = fdf->n;
  const size_t p = fdf->p;
  const double hinv = 1.0 / h;
  size_t i;

  /* compute work = x + h*v */
  for (i = 0; i < p; ++i)
  {
    double xi = gsl_vector_get(x, i);
    double vi = gsl_vector_get(v, i);

    gsl_vector_set(work, i, xi + h * vi);
  }

  /* compute f(x + h*v) */
  status = gsl_multifit_nlinear_eval_f_LD(fdf, work, Dw, Lw, fvv);
  if (status)
    return status;

  for (i = 0; i < n; ++i)
  {
    double fi = gsl_vector_get(f, i);    /* f_i(x) */
    double fip = gsl_vector_get(fvv, i); /* f_i(x + h*v) */
    gsl_vector_const_view row = gsl_matrix_const_row(J, i);
    double u, fvvi;

    /* compute u = sum_{ij} J_{ij} D v_j */
    gsl_blas_ddot(&row.vector, v, &u);

    fvvi = (2.0 * hinv) * ((fip - fi) * hinv - u);

    gsl_vector_set(fvv, i, fvvi);
  }

  return status;
}
