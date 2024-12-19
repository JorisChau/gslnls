#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

/*
gsl_multifit_nlinear_driver2()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter  - maximum iterations to allow
        xtol     - tolerance in step x
        gtol     - tolerance in gradient
        ftol     - tolerance in ||f||
        callback - callback function to call each iteration
        callback_params - parameters to pass to callback function
        chisq0   - ssr previous iteration
        chisq1   - ssr current iteration
        info     - (output) info flag on why iteration terminated
                   1 = stopped due to small step size ||dx|
                   2 = stopped due to small gradient
                   3 = stopped due to small change in f
                   GSL_ETOLX = ||dx|| has converged to within machine
                               precision (and xtol is too small)
                   GSL_ETOLG = ||g||_inf is smaller than machine
                               precision (gtol is too small)
                   GSL_ETOLF = change in ||f|| is smaller than machine
                               precision (ftol is too small)
        lu       - (2 x p)-matrix with top row the lower
                    parameter bounds and bottom row the
                    upper parameter bounds
        w        - workspace

Return:
GSL_SUCCESS if converged
GSL_EBADFUNC if function evaluation failed
GSL_MAXITER if maxiter exceeded without converging
GSL_ENOPROG if no accepted step found on first iteration
*/
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
                                        gsl_multifit_nlinear_workspace *w)
{
    int status = GSL_CONTINUE;
    R_len_t iter = 0;
    gsl_vector *f = NULL;

    do
    {
        /* current ssr */
        chisq0[0] = chisq1[0];

        /* iterator */
        if (!lu)
            status = gsl_multifit_nlinear_iterate(w);
        else // with parameter constraints
        {
            status = trust_iterate_lu(w->state, w->sqrt_wts, w->fdf, w->x, w->f, w->J, w->g, w->dx, lu);
            w->niter++;
        }

        /* new ssr */
        f = gsl_multifit_nlinear_residual(w);
        gsl_blas_ddot(f, f, chisq1);

        if (callback)
            ((fdata *)callback_params)->chisq = chisq1[0];

        /*
         * If the solver reports no progress on the first iteration,
         * then it didn't find a single step to reduce the
         * cost function and more iterations won't help so return.
         *
         * If we get a no progress flag on subsequent iterations,
         * it means we did find a good step in a previous iteration,
         * so continue iterating since the solver has now reset
         * mu to its initial value.
         */
        if (status == GSL_EBADFUNC || (status == GSL_ENOPROG && iter == 0))
        {
            *info = status;
            return status;
        }

        ++iter;

        if (callback)
            callback(iter, callback_params, w);

        /* test for convergence */
        status = gsl_multifit_nlinear_test(xtol, gtol, ftol, info, w);
    } while (status == GSL_CONTINUE && iter < maxiter);

    /*
     * the following error codes mean that the solution has converged
     * to within machine precision, so record the error code in info
     * and return success
     */
    if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
        *info = status;
        status = GSL_SUCCESS;
    }

    /* check if max iterations reached */
    if (iter >= maxiter && status != GSL_SUCCESS)
        status = GSL_EMAXITER;

    return status;
} /* gsl_multifit_nlinear_driver() */

/*
gsl_multilarge_nlinear_driver2()
  Iterate the large-scale nonlinear least squares solver until completion

Inputs: maxiter  - maximum iterations to allow
        xtol     - tolerance in step x
        gtol     - tolerance in gradient
        ftol     - tolerance in ||f||
        callback - callback function to call each iteration
        callback_params - parameters to pass to callback function
        chisq0   - ssr previous iteration
        chisq1   - ssr current iteration
        info     - (output) info flag on why iteration terminated
                   1 = stopped due to small step size ||dx|
                   2 = stopped due to small gradient
                   3 = stopped due to small change in f
                   GSL_ETOLX = ||dx|| has converged to within machine
                               precision (and xtol is too small)
                   GSL_ETOLG = ||g||_inf is smaller than machine
                               precision (gtol is too small)
                   GSL_ETOLF = change in ||f|| is smaller than machine
                               precision (ftol is too small)
        w        - workspace

Return:
GSL_SUCCESS if converged
GSL_EBADFUNC if function evaluation failed
GSL_MAXITER if maxiter exceeded without converging
GSL_ENOPROG if no accepted step found on first iteration
*/
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
                                          gsl_multilarge_nlinear_workspace *w)
{
    int status = GSL_CONTINUE;
    R_len_t iter = 0;
    gsl_vector *f = NULL;

    do
    {
        /* current ssr */
        chisq0[0] = chisq1[0];

        status = gsl_multilarge_nlinear_iterate(w);

        /* new ssr */
        f = gsl_multilarge_nlinear_residual(w);
        gsl_blas_ddot(f, f, chisq1);

        if (callback)
            ((fdata_large *)callback_params)->chisq = chisq1[0];

        /*
         * If the solver reports no progress on the first iteration,
         * then it didn't find a single step to reduce the
         * cost function and more iterations won't help so return.
         *
         * If we get a no progress flag on subsequent iterations,
         * it means we did find a good step in a previous iteration,
         * so continue iterating since the solver has now reset
         * mu to its initial value.
         */
        if (status == GSL_EBADFUNC || (status == GSL_ENOPROG && iter == 0))
        {
            *info = status;
            return status;
        }

        ++iter;

        if (callback)
            callback(iter, callback_params, w);

        /* test for convergence */
        status = gsl_multilarge_nlinear_test(xtol, gtol, ftol, info, w);
    } while (status == GSL_CONTINUE && iter < maxiter);

    /*
     * the following error codes mean that the solution has converged
     * to within machine precision, so record the error code in info
     * and return success
     */
    if (status == GSL_ETOLF || status == GSL_ETOLX || status == GSL_ETOLG)
    {
        *info = status;
        status = GSL_SUCCESS;
    }

    /* check if max iterations reached */
    if (iter >= maxiter && status != GSL_SUCCESS)
        status = GSL_EMAXITER;

    return status;
} /* gsl_multilarge_nlinear_driver2() */
