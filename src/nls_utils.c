#include "gsl_nls.h"

/*
    several utility functions for general use
*/

/*
det_eval_jtj()
  Evaluate Jacobian J and calculate det(J^T * J) from Cholesky decomposition

  Inputs: vstate - workspace
        params - parameter workspace
        Dw   - diagonal D in LDDL' decomposition of W
        Lw   - unit lower triangular matrix L in LDDL' decomposition of W
        fdf    - user callback functions
        x      - length-p parameter vector
        f      - length-n f(x) vector
        J      - nxp jacobian matrix
        JTJ    - workspace pxp hessian matrix
        workn  - workspace length-n vector

*/
double det_eval_jtj(const gsl_multifit_nlinear_parameters params,
                    const gsl_vector *Dw, const gsl_matrix *Lw,
                    gsl_multifit_nlinear_fdf *fdf,
                    const gsl_vector *x, gsl_vector *f, gsl_matrix *J,
                    gsl_matrix *JTJ, gsl_vector *workn)
{
    int status;
    double det = 0.0;

    /* evaluate f(x) */
    if(Lw)
        status = gsl_multifit_nlinear_eval_f_LD(fdf, x, Dw, Lw, f);
    else 
        status = gsl_multifit_nlinear_eval_f(fdf, x, Dw, f);

    if (status)
        return det;

    /* evaluate J(x) */
    if(Lw) 
        status = gsl_multifit_nlinear_eval_df_LD(x, f, Dw, Lw, params.h_df, params.fdtype, fdf, J, workn);
    else
        status = gsl_multifit_nlinear_eval_df(x, f, Dw, params.h_df, params.fdtype, fdf, J, workn);

    if (status)
        return det;

    det = det_cholesky_jtj(J, JTJ);

    return det;
}

double det_cholesky_jtj(gsl_matrix *J, gsl_matrix *JTJ)
{
    int status;
    double det = 0.0;

    /* compute J^T * J */
    gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, JTJ);

    /* calculate det(J^T * J) from lower Cholesky matrix */
    status = gsl_linalg_cholesky_decomp1(JTJ);
    if (status)
        return det;

    det = 1.0;
    for (R_len_t i = 0; i < (JTJ->size1); ++i)
        det *= gsl_matrix_get(JTJ, i, i); // product of diagonal elements

    return det * det;
}

/*
hat_values()
  Evaluate hat-values as diagonal of J * (J^T * J)^-1 * J^T 

  Inputs: 
        J      - nxp Jacobian matrix
        JTJ    - pxp lower Cholesky of J^T * J
        h      - n-length vector used to write hat-values
        worknp - nxp matrix to store intermediate results

  Output: integer status

*/
int hat_values(gsl_matrix *J, gsl_matrix *JTJ, gsl_vector *h, gsl_matrix *worknp)
{
    // compute (J^T * J)^(-1)
    int status = GSL_SUCCESS;
    status = gsl_linalg_cholesky_invert(JTJ);
    if (status)
        return status;

    // compute diag(J * (J^T * J)^(-1) * J^T)
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, J, JTJ, 0.0, worknp);
    gsl_matrix_mul_elements(worknp, J);

    // diagonal of H
    for (R_len_t i = 0; i < (worknp->size1); ++i)
    {
        double hi = 0.0;
        for (R_len_t j = 0; j < (worknp->size2); ++j)
            hi += gsl_matrix_get(worknp, i, j);
        gsl_vector_set(h, i, hi);
    }

    return status;
}

/*
cooks_d()
  Evaluate Cook's distance as D = e^2 / (p * s^2) * (h / (1 - h)^2)

  Inputs:
        f      - n-length residual vector
        J      - nxp Jacobian matrix
        JTJ    - pxp lower Cholesky of J^T * J
        d      - n-length vector used to write results
        worknp - nxp matrix to store intermediate results

  Output: integer status

*/
int cooks_d(gsl_vector *f, gsl_matrix *J, gsl_matrix *JTJ, gsl_vector *d, gsl_matrix *worknp)
{
    int status = GSL_SUCCESS;

    R_len_t n = J->size1;
    R_len_t p = J->size2;

    double chisq = (double)GSL_POSINF;
    gsl_blas_ddot(f, f, &chisq);
    double s2 = chisq / (n - p);
   
    status = hat_values(J, JTJ, d, worknp);
    if (status)
        return status;

    for (R_len_t i = 0; i < n; ++i)
    {
        double ei = gsl_vector_get(f, i);             // residual i
        double hi = gsl_vector_get(d, i);             // hat-value i
        double di = (ei * ei) / (p * s2) * (hi / ((1 - hi) * (1 - hi)));
        gsl_vector_set(d, i, di);
    }

    return status;
}

/* gsl_median()
   included here for compatibility w/ GSL < 2.5. For GSL >= 2.5, 
   gsl_stats_median() is available as well

   Inputs: 
        data    - n-length numeric array
        n       - integer length 

   Output: double median
*/
double gsl_median(double *data, const R_len_t n)
{
    double median = 0.0;
    if (n == 0)
        return median;

    // order vector
    SEXP data1 = PROTECT(Rf_allocVector(REALSXP, n));
    int *order1 = (int *)R_alloc(n, sizeof(int));
    for (R_len_t i = 0; i < n; i++)
        SET_REAL_ELT(data1, i, data[i]);
    R_orderVector1(order1, n, data1, TRUE, FALSE);
    UNPROTECT(1);

    const R_len_t lhs = (n - 1) / 2;
    const R_len_t rhs = n / 2;

    if (lhs == rhs)
    {
        median = data[order1[lhs]];
    }
    else
    {
        median = (data[order1[lhs]] + data[order1[rhs]]) / 2.0;
    }

    return median;
}

/* gsl_mad()
   included here for compatibility w/ GSL < 2.5. For GSL >= 2.5,
   gsl_stats_mad() is available as well

   Inputs:
        data    - n-length numeric array
        n       - integer length
        work    - n-length workspace array

   Output: median absolute deviation
*/
double gsl_mad(double *data, const R_len_t n, double *work)
{
    double median, mad;

    /* median of input data*/
    median = gsl_median(data, n);

    /* absolute deviation from median */
    for (R_len_t i = 0; i < n; ++i)
        work[i] = fabs((double)data[i] - median);

    mad = 1.482602218505602 * gsl_median(work, n);

    return mad;

}
