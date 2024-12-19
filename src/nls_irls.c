/* C implementations of weight functions w(x) = psi(x)/x
   have been copied from robustbase's lmrob.c. This is easier
   and more efficient than evaluating the high-level function
    .Mwgt.psi1() in nls.c */
#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

static double wgt_huber(double x, double c)
{
    /*
     * w(x) = psi(x)/x for Huber's loss function
     */
    return (fabs(x) >= c) ? c / fabs(x) : 1.;
}

static double psip_huber(double x, double c)
{
    // psi' = rho'' : Second derivative of Huber's loss function
    return (fabs(x) >= c) ? 0. : 1.;
}

static double wgt_barron(double x, double *c)
{
    /*
     * w(x) = psi(x)/x for Barron's family of loss functions
     */
    double alpha = c[0];
    double c2 = c[1] * c[1];
    double xc2 = (x * x) / c2;

    if (fabs(alpha - 2.0) < GSL_SQRT_DBL_EPSILON)
        return 1.0 / c2;
    else if (fabs(alpha) < GSL_SQRT_DBL_EPSILON)
        return 2.0 / (x * x + 2 * c2);
    else if (alpha > -1e8)
        return 1.0 / c2 * pow((xc2 / fabs(alpha - 2.0) + 1), 0.5 * alpha - 1.0);
    else
        return 1.0 / c2 * exp(-0.5 * xc2);
}

static double psip_barron(double x, double *c)
{
    // psi' = rho'' : Second derivative for Barron's family of loss functions
    double alpha = c[0];
    double c2 = c[1] * c[1];
    double x2 = x * x;

    if (fabs(alpha - 2.0) < GSL_SQRT_DBL_EPSILON)
        return 1.0 / c2;
    else if (fabs(alpha) < GSL_SQRT_DBL_EPSILON)
        return -2. * (x2 - 2. * c2) / ((2. * c2 + x2) * (2. * c2 + x2));
    else if (alpha > -1e8)
    {
        double denom = x2 - (alpha - 2.) * c2;
        return (alpha - 2.) * ((alpha - 2.) * c2 - (alpha - 1.) * x2) * pow(1. - x2 / ((alpha - 2.) * c2), 0.5 * alpha) / (denom * denom);
    }
    else return exp(-x2 / (2. * c2)) * (c2 - x2) / (c2 * c2);
}

static double wgt_biwgt(double x, double c)
{
    /*
     * w(x) = psi(x)/x for Tukey's bisquare loss function
     */
    if (fabs(x) > c)
        return (0.);
    else
    {
        double a = x / c;
        a = (1. - a) * (1. + a);
        return (a * a);
    }
}

static double psip_biwgt(double x, double c)
{
    /*
     * Second derivative of Tukey's bisquare loss function
     */
    if (fabs(x) > c)
        return (0.);
    else
    {
        x /= c;
        double x2 = x * x;
        return ((1. - x2) * (1 - 5 * x2));
    }
}

static double wgt_gwgt(double x, double c)
{
    /*
     * w(x) = psi(x)/x for Gauss Weight (Welsch) Loss function
     */
    double a = x / c;
    return (exp(-(a * a) / 2));
}

static double psip_gwgt(double x, double c)
{
    /*
     * Gauss Weight Psi'()
     */

    // Largest x  such that  exp(-x^2/2) does not underflow :
    double MAX_Ex2 = 37.7; // ~ = sqrt(- 2. * M_LN2 * DBL_MIN_EXP);
    x /= c;
    if (fabs(x) > MAX_Ex2)
        return 0.;
    else
    {
        double ac = x * x;
        return exp(-ac / 2) * (1. - ac);
    }
}

static double wgt_opt(double x, double c)
{
    /*
     * w(x) = psi(x)/x for optimal psi Function
     */
    double ac = x / c;
    double ax = fabs(ac);
    if (ax > 3.)
        return 0.;
    else if (ax > 2.)
    {
        const double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
        ax *= ax; // = |x/c| ^ 2
        return GSL_MAX(0., R1 + ax * (R2 + ax * (R3 + ax * R4)));
    }
    else
        return 1.;
}

static double psip_opt(double x, double c)
{
    /*
     * psi'() for Optimal psi Function
     */
    double ac = x / c,
           ax = fabs(ac);
    if (ax > 3.)
        return 0;
    else if (ax > 2.)
    {
        const double R1 = -1.944, R2 = 1.728, R3 = -0.312, R4 = 0.016;
        ax *= ax; // = |x/c| ^ 2
        return R1 + ax * (3 * R2 + ax * (5 * R3 + ax * 7 * R4));
    }
    else
        return 1.;
}


static double wgt_hmpl(double x, double *k)
{
    /*
     * w(x) = psi(x)/x  for Hampel's redescending psi function
     * Hampel redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
        return (1);
    else if (u <= k[1])
        return (k[0] / u);
    else if (u <= k[2])
        return (k[0] * (k[2] - u) / (k[2] - k[1]) / u);
    else
        return (0);
}

static double psip_hmpl(double x, double *k)
{
    /*
     * psi'()  for Hampel's redescending psi function
     * constants  (a, b, r) == k[0:2]   s.t. slope of psi is 1 in the center
     */
    double u = fabs(x);

    if (u <= k[0])
        return (1);
    else if (u <= k[1])
        return (0);
    else if (u <= k[2])
        return (k[0] / (k[1] - k[2]));
    else
        return (0);
}

static double wgt_ggw(double x, double *k)
{
    /*
     * w(x) = psi(x)/x for Gauss Weight psi function with constant center
     */
    double ax = fabs(x);
    return (ax < k[2]) ? 1. : exp(-pow(ax - k[2], k[1]) / 2 / k[0]);
}

static double psip_ggw(double x, double *k)
{
    /*
     * Gauss Weight with constant center
     */
    double ax = fabs(x);
    if (ax < k[2])
        return 1.;
    else
    {
        // Largest x  such that  exp(-x) does not underflow :
        double MIN_Exp = -708.4; // ~ = M_LN2 * DBL_MIN_EXP = -log(2) * 1022 = -708.3964 */
        double a, b, c, ea;
        a = 2 * k[0];
        b = k[1];
        c = k[2];
        ea = -pow(ax - c, b) / a;
        return (ea < MIN_Exp) ? 0. : exp(ea) * (1 - b / a * ax * pow(ax - c, b - 1));
    }
}

static double wgt_lqq(double x, double *k)
{
    /*
     * w(x) = psi(x)/x for Linear-Quadratic-Quadratic ("lqq") function
     */
    double ax = fabs(x);
    if (ax <= k[1])
        return (1.);
    else
    {
        double k01 = k[0] + k[1];
        if (ax <= k01)
        {
            double s0 = ax - k[1];
            return (1. - k[2] * s0 * s0 / (2 * ax * k[0]));
        }
        else
        {
            double
                s5 = k[2] - 1.,
                s6 = -2 * k01 + k[0] * k[2];
            if (ax < k01 - s6 / s5)
            {
                double s7 = ax - k01;
                return (-(s6 / 2. + s5 * s5 / s6 * s7 * (s7 / 2. + s6 / s5)) / ax);
            }
            else
                return (0.);
        }
    }
}

static double psip_lqq(double x, double *k)
{
    double ax = fabs(x);
    if (ax <= k[1])
        return (1.);
    else
    {
        double k01 = k[0] + k[1]; // = b+c
        if (/*k[1] < ax && */ ax <= k01)
            return 1. - k[2] / k[0] * (ax - k[1]);
        else
        {
            double
                s5 = 1. - k[2], // = (1-s)
                a = (k[0] * k[2] - 2 * k01) / s5;
            if (/* k01 < ax && */ ax < k01 + a)
                return -s5 * ((ax - k01) / a - 1.);
            else
                return 0.;
        }
    }
}

double wgt(double x, double *c, int i)
{
    /*
     * return the correct wgt according to i
     * wgt: rho'(x) / x = psi(x) / x
     */
    switch (i)
    {
    default:
    case 1:
        return (wgt_huber(x, c[0])); // huber
    case 2:
        return (wgt_barron(x, c));  // barron
    case 3:
        return (wgt_biwgt(x, c[0])); // biweight
    case 4:
        return (wgt_gwgt(x, c[0])); // GaussWeight / "Welsh"
    case 5:
        return (wgt_opt(x, c[0])); // Optimal
    case 6:
        return (wgt_hmpl(x, c)); // Hampel
    case 7:
        return (wgt_ggw(x, c)); // GGW
    case 8:
        return (wgt_lqq(x, c)); // LQQ (piecewise linear psi')
    }
}

double psip(double x, double *c, int i)
{
    /*
     * return the correct rho'' according to i
     */
    switch (i)
    {
    default:
    case 1:
        return (psip_huber(x, c[0])); // huber
    case 2:
        return (psip_barron(x, c));  // barron
    case 3:
        return (psip_biwgt(x, c[0])); // biweight
    case 4:
        return (psip_gwgt(x, c[0])); // GaussWeight / "Welsh"
    case 5:
        return (psip_opt(x, c[0])); // Optimal
    case 6:
        return (psip_hmpl(x, c)); // Hampel
    case 7:
        return (psip_ggw(x, c)); // GGW
    case 8:
        return (psip_lqq(x, c)); // LQQ (piecewise linear psi')
    }
}

static int test_delta_irls(const gsl_vector *x0, const gsl_vector *x1, double xtol)
{
    int status;
    R_len_t n = x1->size;

    for (R_len_t i = 0; i < n; i++)
    {
        double xi = gsl_vector_get(x1, i);
        double dxi = fabs(gsl_vector_get(x0, i) - xi);
        if (gsl_min(dxi / fabs(xi), dxi) < xtol)
            status = GSL_SUCCESS;
        else
        {
            status = GSL_CONTINUE;
            break;
        }
    }

    return status;
}

static void callback_irls(const R_len_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    /* update traces */
    double chisq = ((fdata *)params)->chisq;
    SET_REAL_ELT(((fdata *)params)->ssrtrace, iter, chisq);
    R_len_t p = ((fdata *)params)->p;
    R_len_t n = (R_len_t)Rf_nrows(((fdata *)params)->partrace);
    double *parptr = REAL(((fdata *)params)->partrace);
    for (R_len_t k = 0; k < p; k++)
        parptr[iter + n * k] = gsl_vector_get(w->x, k);
}

/*
gsl_multifit_nlinear_rho_driver()
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
        wgt_i    - psi weight function identifier
        wgt_cc   - psi weight function parameter(s)
        irls_iter - number of IRLS iterations
        irls_delta - achieved IRLS tolerance
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
int gsl_multifit_nlinear_rho_driver(
    pdata *pars,
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
    double *irls_sigma,
    R_len_t *irls_iter,
    int *irls_status,
    Rboolean verbose)
{
    /* initialize parameters */
    SEXP wgt_cc = PROTECT(VECTOR_ELT(pars->loss_config, 1));
    double *wgt_cc_ptr = REAL(wgt_cc);
    R_len_t irls_maxiter = INTEGER_ELT(pars->control_int, 14);
    double irls_xtol = REAL_ELT(pars->control_dbl, 10);
    R_len_t n = ((fdata *) fdff->params)->n;
    R_len_t p = ((fdata *) fdff->params)->p;
    int status = GSL_CONTINUE;
    double *resid = (double *)S_alloc(n, sizeof(double));
    double *abs_resid = (double *)S_alloc(n, sizeof(double));

    do
    {
        *irls_iter += 1;

        if (*irls_iter > 1) 
        {
            /* re-initialize workspace */
            gsl_vector_memcpy(pars->workp, (pars->w)->x);
            gsl_multifit_nlinear_winit(pars->mpopt, pars->wts, fdff, pars->w);
            *chisq0 = (double)GSL_POSINF;
        }
        else 
            gsl_vector_memcpy(pars->workp, pars->mpopt); 
        
        /* WLS optimization */
        status = gsl_multifit_nlinear_driver2(
            maxiter, xtol, gtol, ftol, verbose ? callback_irls : NULL,
            verbose ? callback_params : NULL, info, chisq0, chisq1,
            pars->lu, pars->w);

        if (verbose)
        {
            /* print trace */
            Rprintf("IRLS iter %3d: ssr = %g, par = (", *irls_iter, *chisq1);
            for (R_len_t k = 0; k < p; k++)
                Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", gsl_vector_get((pars->w)->x, k));
        }

        /*
         * If the solver reports no progress on the first iteration,
         * then it didn't find a single step to reduce the
         * cost function and more iterations won't help so return.
         */
        if (status == GSL_EBADFUNC || (status == GSL_ENOPROG && *irls_iter == 1))
        {
            *info = status;
            return status;
        }

        /* update weights */
        for (R_len_t i = 0; i < n; i++)
        {
            resid[i] = gsl_vector_get((pars->w)->f, i) / gsl_vector_get((pars->w)->sqrt_wts, i);
            abs_resid[i] = fabs(resid[i]);
        }
        *irls_sigma = 1.48258 * gsl_stats_median(abs_resid, 1, n);

        for (R_len_t i = 0; i < n; i++)
        {
            double residsc = resid[i] / *irls_sigma;
            double weight = wgt(residsc, wgt_cc_ptr, wgt_i);
            gsl_vector_set(pars->wts, i, gsl_max(weight, GSL_DBL_EPSILON));
            gsl_vector_set(pars->psi, i, weight * residsc);
            gsl_vector_set(pars->psip, i, psip(residsc, wgt_cc_ptr, wgt_i));
            if (!Rf_isNull(pars->weights))
                gsl_vector_set(pars->wts, i, REAL_ELT(pars->weights, i) * gsl_vector_get(pars->wts, i));
        }

        /* check irls convergence */
        *irls_status = test_delta_irls(pars->workp, (pars->w)->x, irls_xtol);

        if(*irls_status == GSL_SUCCESS)
        {
            *info = status;
            return status;
        }

        // Rprintf("cc = %g, scale = %g, i = %d\n", wgt_cc_ptr[0], scale, wgt_i);
        // Rprintf("resid = (");
        // for (R_len_t i = 0; i < n; i++)
        //     Rprintf((i < (n - 1)) ? "%g, " : "%g)\n", resid[i] / scale);
        // Rprintf("weight = (");
        // for (R_len_t i = 0; i < n; i++)
        //     Rprintf((i < (n - 1)) ? "%g, " : "%g)\n", gsl_vector_get(pars->wts, i));

    } while (*irls_status == GSL_CONTINUE && *irls_iter < irls_maxiter);

    UNPROTECT(1);

    /* check if max iterations reached */
    if (*irls_iter >= irls_maxiter && *irls_status != GSL_SUCCESS)
    {
        *irls_status = GSL_EMAXITER;
        *info = GSL_EMAXITER;
        status = GSL_EMAXITER;
    }

    return status;
} /* gsl_multifit_nlinear_rho_driver() */
