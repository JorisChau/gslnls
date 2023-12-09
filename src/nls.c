#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

/* cleanup memory */
static void C_nls_cleanup(void *data)
{
    pdata *pars = data;

    /* free memory */
    if (pars->w)
        gsl_multifit_nlinear_free(pars->w);
    if (pars->q)
        gsl_qrng_free(pars->q);
    if (pars->wts)
        gsl_vector_free(pars->wts);
    if (pars->mx)
        gsl_matrix_free(pars->mx);
    if (pars->mp)
        gsl_vector_free(pars->mp);
    if (pars->mpopt)
        gsl_vector_free(pars->mpopt);
    if (pars->diag)
        gsl_vector_free(pars->diag);
}

/* function call w/ cleanup */
SEXP C_nls(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl)
{
    /* function arguments */
    pdata pars = {fn, y, jac, fvv, env, start, swts, control_int, control_dbl, NULL, NULL, NULL, NULL, NULL, NULL};

    /* safe function call */
    SEXP ans = R_ExecWithCleanup(C_nls_internal, &pars, C_nls_cleanup, &pars);

    return ans;
}

SEXP C_nls_internal(void *data)
{
    /* turn off error handler to avoid aborting on errors */
    gsl_set_error_handler_off();

    /* function arguments */
    pdata *pars = data;

    /* initialize parameters */
    int nprotect = 0;
    R_len_t n = Rf_length(pars->y);
    R_len_t niter = INTEGER_ELT(pars->control_int, 0);
    int verbose = INTEGER_ELT(pars->control_int, 1);
    Rboolean mstart = Rf_isMatrix(pars->start);
    SEXP startvec = NULL;
    R_len_t p;
    if (!mstart)
    {
        startvec = PROTECT(Rf_coerceVector(pars->start, REALSXP));
        p = Rf_length(startvec);
        nprotect++;
    }
    else
    {
        p = Rf_ncols(pars->start);
    }

    /* control parameters */
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    /* minimization algorithm */
    switch (INTEGER_ELT(pars->control_int, 2))
    {
    case 1:
        fdf_params.trs = gsl_multifit_nlinear_trs_lmaccel;
        break;
    case 2:
        fdf_params.trs = gsl_multifit_nlinear_trs_dogleg;
        break;
    case 3:
        fdf_params.trs = gsl_multifit_nlinear_trs_ddogleg;
        break;
    case 4:
        fdf_params.trs = gsl_multifit_nlinear_trs_subspace2D;
        break;
    default:
        fdf_params.trs = gsl_multifit_nlinear_trs_lm;
    }

    /* scaling method */
    switch (INTEGER_ELT(pars->control_int, 3))
    {
    case 1:
        fdf_params.scale = gsl_multifit_nlinear_scale_levenberg;
        break;
    case 2:
        fdf_params.scale = gsl_multifit_nlinear_scale_marquardt;
        break;
    default:
        fdf_params.scale = gsl_multifit_nlinear_scale_more;
    }

    /* solver type */
    switch (INTEGER_ELT(pars->control_int, 4))
    {
    case 1:
        fdf_params.solver = gsl_multifit_nlinear_solver_cholesky;
        break;
    case 2:
        fdf_params.solver = gsl_multifit_nlinear_solver_svd;
        break;
    default:
        fdf_params.solver = gsl_multifit_nlinear_solver_qr;
    }

    /* finite differencing type */
    fdf_params.fdtype = INTEGER_ELT(pars->control_int, 5) ? GSL_MULTIFIT_NLINEAR_CTRDIFF : GSL_MULTIFIT_NLINEAR_FWDIFF;

    /* miscellaneous parameters */
    fdf_params.factor_up = REAL_ELT(pars->control_dbl, 0);
    fdf_params.factor_down = REAL_ELT(pars->control_dbl, 1);
    fdf_params.avmax = REAL_ELT(pars->control_dbl, 2);
    fdf_params.h_df = REAL_ELT(pars->control_dbl, 3);
    fdf_params.h_fvv = REAL_ELT(pars->control_dbl, 4);
    double xtol = REAL_ELT(pars->control_dbl, 5);
    double ftol = REAL_ELT(pars->control_dbl, 6);
    double gtol = REAL_ELT(pars->control_dbl, 7);

    /* initialize data */
    SEXP xpar = Rf_install("par");
    SEXP fcall = PROTECT(Rf_lang2(pars->fn, xpar));
    SEXP parnames = NULL;
    if (!mstart)
        parnames = PROTECT(Rf_getAttrib(pars->start, R_NamesSymbol));
    else
        parnames = PROTECT(VECTOR_ELT(Rf_getAttrib(pars->start, R_DimNamesSymbol), 1));
    nprotect += 2;

    fdata params;
    params.n = n;
    params.p = p;
    params.f = fcall;
    params.df = NULL;
    params.fvv = NULL;
    params.y = pars->y;
    params.rho = pars->env;
    params.start = pars->start;
    if (!mstart)
        params.warn = TRUE;
    else
        params.warn = FALSE;

    if (verbose)
    {
        params.partrace = PROTECT(Rf_allocMatrix(REALSXP, niter + 1, p));
        params.ssrtrace = PROTECT(Rf_allocVector(REALSXP, niter + 1));
        nprotect += 2;
    }

    /* define the function to be minimized */
    gsl_multifit_nlinear_fdf fdf;
    fdf.f = gsl_f;
    fdf.df = NULL;  // set to NULL for finite-difference Jacobian
    fdf.fvv = NULL; // not using geodesic acceleration
    fdf.n = n;
    fdf.p = p;
    fdf.nevalf = 0;
    fdf.nevaldf = 0;
    fdf.nevalfvv = 0;
    fdf.params = &params;

    /* use Jacobian function */
    if (!Rf_isNull(pars->jac))
    {
        SEXP dfcall = PROTECT(Rf_lang2(pars->jac, xpar));
        params.df = dfcall;
        fdf.df = gsl_df;
        nprotect++;
    }

    /* use acceleration function */
    if (!Rf_isNull(pars->fvv))
    {
        SEXP vpar = Rf_install("v");
        SEXP fvvcall = PROTECT(Rf_lang3(pars->fvv, xpar, vpar));
        params.fvv = fvvcall;
        fdf.fvv = gsl_fvv;
        nprotect++;
    }

    /* use weights */
    pars->wts = gsl_vector_alloc(n);
    if (!Rf_isNull(pars->swts))
        for (R_len_t i = 0; i < n; i++)
            gsl_vector_set(pars->wts, i, REAL_ELT(pars->swts, i));

    /* initialize solver */
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

    /* allocate workspace with default parameters */
    pars->w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);
    pars->mpopt = gsl_vector_alloc(p);

    /* multistart algorithm */
    if (mstart)
    {
        /* allocate qrng workspace */
        if (p < 41)
            pars->q = gsl_qrng_alloc(gsl_qrng_sobol, p);
        else
            pars->q = gsl_qrng_alloc(gsl_qrng_halton, p);

        /* multistart parameters */
        mdata mpars;
        mpars.n = INTEGER_ELT(pars->control_int, 6);
        mpars.p = INTEGER_ELT(pars->control_int, 7);
        mpars.q = INTEGER_ELT(pars->control_int, 8);
        mpars.s = INTEGER_ELT(pars->control_int, 9);
        mpars.niter = INTEGER_ELT(pars->control_int, 10);
        mpars.max = INTEGER_ELT(pars->control_int, 11);
        mpars.minsp = INTEGER_ELT(pars->control_int, 12);
        mpars.r = REAL_ELT(pars->control_dbl, 8);
        mpars.tol = REAL_ELT(pars->control_dbl, 9);
        mpars.ntix = (int *)S_alloc(mpars.n, sizeof(int));
        mpars.qmp = (double *)S_alloc(p, sizeof(double));
        mpars.startptr = REAL(pars->start);
        mpars.mssr_order = (int *)R_alloc(mpars.n, sizeof(int));
        mpars.mstop = GSL_CONTINUE;
        mpars.mstarts = 0;
        mpars.nsp = 0;
        mpars.nwsp = 0;
        mpars.mssropt = GSL_POSINF;

        SEXP mssr = PROTECT(Rf_allocVector(REALSXP, mpars.n));
        SEXP mpar = PROTECT(Rf_allocVector(REALSXP, p));
        Rf_setAttrib(mpar, R_NamesSymbol, parnames);
        params.start = mpar;
        nprotect += 2;

        /* placeholder matrices/vectors */
        pars->mx = gsl_matrix_alloc(mpars.n, p);
        pars->mp = gsl_vector_alloc(p);
        pars->diag = gsl_vector_alloc(p);
        /* initially focus sampling around center points */
        for (R_len_t k = 0; k < p; k++)
            gsl_vector_set(pars->diag, k, (mpars.startptr[2 * k + 1] - mpars.startptr[2 * k]) / 20.0);

        /* multi-start global iterations */
        do
        {
            gsl_multistart_driver(pars, &mpars, &fdf, mssr, xtol, ftol, gtol, verbose);

            /* check stopping criterion */
            mpars.mstarts += 1;
            if (mpars.mstarts > mpars.max)
                mpars.mstop = GSL_EMAXITER;
            if (mpars.nsp >= mpars.minsp && mpars.nwsp >= mpars.r * mpars.nsp)
                mpars.mstop = GSL_SUCCESS;

        } while (mpars.mstop == GSL_CONTINUE);
        if (verbose)
        {
            if (mpars.mstop == GSL_SUCCESS)
                Rprintf("multi-start algorithm finished successfully (NSP = %d, NWSP = %d, # iterations = %d)\n", mpars.nsp, mpars.nwsp, mpars.mstarts);
            if (mpars.mstop == GSL_EMAXITER)
                Rprintf("multi-start algorithm reached max. number of global iterations (NSP = %d, NWSP = %d, # iterations = %d)\n", mpars.nsp, mpars.nwsp, mpars.mstarts);
            Rprintf("*******************\n");
        }
        /* add noise in case of zero residuals */
        if(mpars.mssropt == 0)
        {
            gsl_vector_set(pars->mpopt, 0, gsl_vector_get(pars->mpopt, 0) + xtol);
        }
    }
    else /* single-start */
    {
        for (R_len_t k = 0; k < p; k++)
            gsl_vector_set(pars->mpopt, k, REAL_ELT(startvec, k));
    }

    /* (re-)initialize solver w/ optimal start parameters */
    if (!Rf_isNull(pars->swts))
        gsl_multifit_nlinear_winit(pars->mpopt, pars->wts, &fdf, pars->w);
    else
        gsl_multifit_nlinear_init(pars->mpopt, &fdf, pars->w);

    /* initial cost */
    double chisq_init = GSL_POSINF;
    gsl_vector *resid = gsl_multifit_nlinear_residual(pars->w);
    gsl_blas_ddot(resid, resid, &chisq_init);
    double chisq0 = chisq_init;
    double chisq1 = chisq_init;
    params.chisq = chisq_init;

    if (verbose)
    {
        SET_REAL_ELT(params.ssrtrace, 0, chisq_init);
        double *parptr = REAL(params.partrace);
        for (R_len_t k = 0; k < p; k++)
            parptr[(niter + 1) * k] = gsl_vector_get(pars->mpopt, k);
    }

    /* solve the system  */
    int info = GSL_CONTINUE;
    int status = gsl_multifit_nlinear_driver2(niter, xtol, gtol, ftol, verbose ? callback : NULL, verbose ? &params : NULL, &info, &chisq0, &chisq1, pars->w);
    R_len_t iter = gsl_multifit_nlinear_niter(pars->w);

    /* compute covariance and cost at best fit parameters */
    gsl_matrix *J = NULL;
    gsl_matrix *cov = NULL;
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        cov = gsl_matrix_alloc(p, p);
        J = gsl_multifit_nlinear_jac(pars->w);
        gsl_multifit_nlinear_covar(J, 0.0, cov);
    }

    if (verbose)
    {
        /* print summary statistics*/
        Rprintf("*******************\nsummary from method 'multifit/%s'\n", gsl_multifit_nlinear_trs_name(pars->w));
        Rprintf("number of iterations: %d\n", iter);
        Rprintf("reason for stopping: %s\n", gsl_strerror(info));
        Rprintf("initial ssr = %g\n", chisq_init);
        Rprintf("final ssr = %g\n", chisq1);
        Rprintf("ssr/dof = %g\n", chisq1 / (n - p));
        Rprintf("ssr achieved tolerance = %g\n", chisq0 - chisq1);
        Rprintf("function evaluations: %d\n", (int)fdf.nevalf);
        Rprintf("jacobian evaluations: %d\n", (int)fdf.nevaldf);
        Rprintf("fvv evaluations: %d\n", (int)fdf.nevalfvv);
        Rprintf("status = %s\n*******************\n", gsl_strerror(status));
    }

    /* initialize result */
    SEXP ans = NULL;
    if (verbose)
    {
        const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
                                "algorithm", "neval", "partrace", "ssrtrace", ""};
        ans = PROTECT(Rf_mkNamed(VECSXP, ansnms));
    }
    else
    {
        const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
                                "algorithm", "neval", ""};
        ans = PROTECT(Rf_mkNamed(VECSXP, ansnms));
    }
    nprotect++;

    /* estimated parameters */
    SEXP anspar = PROTECT(Rf_allocVector(REALSXP, p));
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(anspar, k, gsl_vector_get((pars->w)->x, k));
    }
    else
    {
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(anspar, k, gsl_vector_get(pars->mpopt, k));
    }
    if (!Rf_isNull(parnames))
        Rf_setAttrib(anspar, R_NamesSymbol, parnames);
    SET_VECTOR_ELT(ans, 0, anspar);
    UNPROTECT(1);

    /* covariance matrix */
    SEXP anscov = PROTECT(Rf_allocMatrix(REALSXP, p, p));
    double *covptr = REAL(anscov);
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        for (R_len_t k1 = 0; k1 < p; k1++)
        {
            for (R_len_t k2 = 0; k2 < p; k2++)
                covptr[k1 + p * k2] = gsl_matrix_get(cov, k1, k2);
        }
    }
    else
    {
        for (R_len_t k1 = 0; k1 < p; k1++)
        {
            for (R_len_t k2 = 0; k2 < p; k2++)
                covptr[k1 + p * k2] = NA_REAL;
        }
    }
    if (!Rf_isNull(parnames))
    {
        SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, parnames);
        SET_VECTOR_ELT(dimnames, 1, parnames);
        Rf_setAttrib(anscov, R_DimNamesSymbol, dimnames);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(ans, 1, anscov);
    UNPROTECT(1);

    /* residuals */
    SEXP ansresid = PROTECT(Rf_allocVector(REALSXP, n));
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        for (R_len_t i = 0; i < n; i++)
            SET_REAL_ELT(ansresid, i, gsl_vector_get(resid, i));
    }
    else
    {
        for (R_len_t i = 0; i < n; i++)
            SET_REAL_ELT(ansresid, i, NA_REAL);
    }
    SET_VECTOR_ELT(ans, 2, ansresid);
    UNPROTECT(1);

    /* jacobian matrix */
    SEXP ansjac = PROTECT(Rf_allocMatrix(REALSXP, n, p));
    double *jacptr = REAL(ansjac);
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        for (R_len_t i = 0; i < n; i++)
        {
            for (R_len_t k = 0; k < p; k++)
                jacptr[i + n * k] = gsl_matrix_get(J, i, k);
        }
    }
    else
    {
        for (R_len_t i = 0; i < n; i++)
        {
            for (R_len_t k = 0; k < p; k++)
                jacptr[i + n * k] = NA_REAL;
        }
    }
    if (!Rf_isNull(parnames))
    {
        SEXP colnames = PROTECT(Rf_allocVector(VECSXP, 2));
        SET_VECTOR_ELT(colnames, 1, parnames);
        Rf_setAttrib(ansjac, R_DimNamesSymbol, colnames);
        UNPROTECT(1);
    }
    SET_VECTOR_ELT(ans, 3, ansjac);
    UNPROTECT(1);

    /* miscellaneous */
    SET_VECTOR_ELT(ans, 4, Rf_ScalarInteger(iter));
    SET_VECTOR_ELT(ans, 5, Rf_ScalarString(Rf_mkChar(gsl_strerror(status))));
    SET_VECTOR_ELT(ans, 6, Rf_ScalarInteger(status));
    SET_VECTOR_ELT(ans, 7, Rf_ScalarReal(chisq1));
    SET_VECTOR_ELT(ans, 8, Rf_ScalarReal(chisq0 - chisq1));
    SET_VECTOR_ELT(ans, 9, Rf_ScalarString(Rf_mkChar(gsl_multifit_nlinear_trs_name(pars->w))));

    const char *nms[] = {"f", "J", "fvv", ""};
    SEXP ansneval = PROTECT(Rf_mkNamed(INTSXP, nms));
    SET_INTEGER_ELT(ansneval, 0, fdf.nevalf);
    SET_INTEGER_ELT(ansneval, 1, fdf.nevaldf);
    SET_INTEGER_ELT(ansneval, 2, fdf.nevalfvv);
    SET_VECTOR_ELT(ans, 10, ansneval);
    UNPROTECT(1);

    if (verbose)
    {
        if (!Rf_isNull(parnames))
        {
            SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
            SET_VECTOR_ELT(dimnames, 0, R_NilValue);
            SET_VECTOR_ELT(dimnames, 1, parnames);
            Rf_setAttrib(params.partrace, R_DimNamesSymbol, dimnames);
            UNPROTECT(1);
        }
        SET_VECTOR_ELT(ans, 11, params.partrace);
        SET_VECTOR_ELT(ans, 12, params.ssrtrace);
    }

    /* free memory */
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
        gsl_matrix_free(cov);

    UNPROTECT(nprotect);
    return ans;
}

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
                                 gsl_multifit_nlinear_workspace *w)
{
    int status = GSL_CONTINUE;
    R_len_t iter = 0;
    gsl_vector *f = NULL;

    do
    {
        /* current ssr */
        chisq0[0] = chisq1[0];

        status = gsl_multifit_nlinear_iterate(w);

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

int gsl_f(const gsl_vector *x, void *params, gsl_vector *f)
{
    /* construct parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    if (Rf_isNumeric(start))
    {
        par = PROTECT(Rf_allocVector(REALSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(par, k, gsl_vector_get(x, k));
    }
    else
    {
        par = PROTECT(Rf_allocVector(VECSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_VECTOR_ELT(par, k, Rf_ScalarReal(gsl_vector_get(x, k)));
    }
    Rf_setAttrib(par, R_NamesSymbol, Rf_getAttrib(start, R_NamesSymbol));

    /* evaluate function f */
    SETCADR(((fdata *)params)->f, par);
    SEXP fval = PROTECT(Rf_eval(((fdata *)params)->f, ((fdata *)params)->rho));

    /* function checks */
    R_len_t n = ((fdata *)params)->n;
    if (TYPEOF(fval) != REALSXP || Rf_length(fval) != n)
    {
        if (((fdata *)params)->warn)
            Rf_warning("Evaluating fn does not return numeric vector of expected length n");
        UNPROTECT(2);
        return GSL_EBADFUNC;
    }

    /* set gsl residuals */
    double *fvalptr = REAL(fval);
    double *yptr = REAL(((fdata *)params)->y);
    for (R_len_t i = 0; i < n; i++)
    {
        if (R_IsNaN(fvalptr[i]) || !R_finite(fvalptr[i]))
            gsl_vector_set(f, i, GSL_POSINF);
        else
            gsl_vector_set(f, i, fvalptr[i] - yptr[i]);
    }

    UNPROTECT(2);
    return GSL_SUCCESS;
}

int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    /* construct parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    if (Rf_isNumeric(start))
    {
        par = PROTECT(Rf_allocVector(REALSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(par, k, gsl_vector_get(x, k));
    }
    else
    {
        par = PROTECT(Rf_allocVector(VECSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_VECTOR_ELT(par, k, Rf_ScalarReal(gsl_vector_get(x, k)));
    }
    Rf_setAttrib(par, R_NamesSymbol, Rf_getAttrib(start, R_NamesSymbol));

    /* evaluate Jacobian function */
    SETCADR(((fdata *)params)->df, par);
    SEXP dfval = PROTECT(Rf_eval(((fdata *)params)->df, ((fdata *)params)->rho));

    /* Jacobian checks */
    R_len_t n = ((fdata *)params)->n;
    if (TYPEOF(dfval) != REALSXP || !Rf_isMatrix(dfval) || Rf_ncols(dfval) != p || Rf_nrows(dfval) != n)
    {
        if (((fdata *)params)->warn)
            Rf_warning("Evaluating jac does not return numeric matrix of dimensions n x p");
        UNPROTECT(2);
        return GSL_EBADFUNC;
    }

    double *jacptr = REAL(dfval);
    for (R_len_t i = 0; i < n; i++)
        for (R_len_t k = 0; k < p; k++)
            if (R_IsNaN(jacptr[i + n * k]) || !R_finite(jacptr[i + n * k]))
            {
                if (((fdata *)params)->warn)
                    Rf_warning("Missing/infinite values not allowed when evaluating jac");
                UNPROTECT(2);
                return GSL_EBADFUNC;
            }

    /* set gsl jacobian matrix */
    for (R_len_t i = 0; i < n; i++)
        for (R_len_t k = 0; k < p; k++)
            gsl_matrix_set(J, i, k, jacptr[i + n * k]);

    UNPROTECT(2);
    return GSL_SUCCESS;
}

int gsl_fvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv)
{
    /* populate parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    SEXP parnames = PROTECT(Rf_getAttrib(start, R_NamesSymbol));
    if (Rf_isNumeric(start))
    {
        par = PROTECT(Rf_allocVector(REALSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(par, k, gsl_vector_get(x, k));
    }
    else
    {
        par = PROTECT(Rf_allocVector(VECSXP, p));
        for (R_len_t k = 0; k < p; k++)
            SET_VECTOR_ELT(par, k, Rf_ScalarReal(gsl_vector_get(x, k)));
    }
    Rf_setAttrib(par, R_NamesSymbol, parnames);

    /* populate v vector */
    SEXP vpar = PROTECT(Rf_allocVector(REALSXP, p));
    for (R_len_t k = 0; k < p; k++)
        SET_REAL_ELT(vpar, k, gsl_vector_get(v, k));
    Rf_setAttrib(vpar, R_NamesSymbol, parnames);

    /* evaluate fvv function */
    SETCADR(((fdata *)params)->fvv, par);
    SETCADDR(((fdata *)params)->fvv, vpar);
    SEXP fvvval = PROTECT(Rf_eval(((fdata *)params)->fvv, ((fdata *)params)->rho));

    /* function checks */
    R_len_t n = ((fdata *)params)->n;
    if (TYPEOF(fvvval) != REALSXP || Rf_length(fvvval) != n)
    {
        if (((fdata *)params)->warn)
            Rf_warning("Evaluating fvv does not return numeric vector of expected length n");
        UNPROTECT(4);
        return GSL_EBADFUNC;
    }

    double *fvvvalptr = REAL(fvvval);
    for (R_len_t i = 0; i < n; i++)
    {
        if (R_IsNaN(fvvvalptr[i]) || !R_finite(fvvvalptr[i]))
        {
            if (((fdata *)params)->warn)
                Rf_warning("Missing/infinite values not allowed when evaluating fvv");
            UNPROTECT(4);
            return GSL_EBADFUNC;
        }
    }

    /* set fvv vector */
    for (R_len_t i = 0; i < n; i++)
        gsl_vector_set(fvv, i, fvvvalptr[i]);

    UNPROTECT(4);
    return GSL_SUCCESS;
}

void callback(const R_len_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
{
    /* update traces */
    double chisq = ((fdata *)params)->chisq;
    SET_REAL_ELT(((fdata *)params)->ssrtrace, iter, chisq);
    R_len_t p = ((fdata *)params)->p;
    R_len_t n = (R_len_t)Rf_nrows(((fdata *)params)->partrace);
    double *parptr = REAL(((fdata *)params)->partrace);
    for (R_len_t k = 0; k < p; k++)
        parptr[iter + n * k] = gsl_vector_get(w->x, k);

    /* print trace */
    Rprintf("iter %3d: ssr = %g, par = (", iter, chisq);
    for (R_len_t k = 0; k < p; k++)
        Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", parptr[iter + n * k]);
}

void gsl_multistart_driver(pdata *pars,
                           mdata *mpars,
                           gsl_multifit_nlinear_fdf *fdff,
                           SEXP mssr,
                           const double xtol,
                           const double ftol,
                           const double gtol,
                           Rboolean verbose)
{
    // initialize variables
    int minfo;
    double kd, l0, l1, kd0, kd1;
    double rcond, mchisq0 = GSL_POSINF, mchisq1 = GSL_POSINF;
    trust_state_t *trust_state = (trust_state_t *)(pars->w->state);

    /* sample initial points */
    for (R_len_t nn = 0; nn < mpars->n; nn++)
    {
        if ((mpars->ntix)[nn] == 0)
        {
            gsl_qrng_get(pars->q, mpars->qmp);
            for (R_len_t k = 0; k < (pars->mp)->size; k++)
            {
                l0 = (mpars->startptr)[2 * k];
                l1 = (mpars->startptr)[2 * k + 1];

                if (l1 > l0)
                {
                    kd = gsl_vector_get(pars->diag, k);
                    kd0 = 1.0 / (1.0 + exp(kd / 2.0));
                    kd1 = 1.0 / (1.0 + exp(-kd / 2.0));
                    // gsl_matrix_set(pars->mx, nn, k, l0 + (l1 - l0) * (mpars->qmp)[k]);
                    (mpars->qmp)[k] = kd0 + (kd1 - kd0) * (mpars->qmp)[k];
                    gsl_matrix_set(pars->mx, nn, k, (l1 + l0) / 2.0 + (l1 - l0) / kd * log((mpars->qmp)[k] / (1 - (mpars->qmp)[k])));
                }
                else
                {
                    gsl_matrix_set(pars->mx, nn, k, l0);
                }
            }
        }

        /* concentrate points */
        gsl_matrix_get_row(pars->mp, pars->mx, nn);
        if (!Rf_isNull(pars->swts))
            gsl_multifit_nlinear_winit(pars->mp, pars->wts, fdff, pars->w);
        else
            gsl_multifit_nlinear_init(pars->mp, fdff, pars->w);
        gsl_multifit_nlinear_driver2(mpars->p, xtol, 1e-3, ftol, NULL, NULL, &minfo, &mchisq0, &mchisq1, pars->w);
        gsl_matrix_set_row(pars->mx, nn, (pars->w)->x);
        if (mchisq1 != GSL_POSINF)
            SET_REAL_ELT(mssr, nn, mchisq1);
        else
            SET_REAL_ELT(mssr, nn, NA_REAL);
    }
    /* reduce sample points */
    R_orderVector1(mpars->mssr_order, mpars->n, mssr, TRUE, FALSE);
    for (R_len_t nn = 0; nn < mpars->n; nn++)
        (mpars->ntix)[(mpars->mssr_order)[nn]] = (nn < (mpars->q)) ? (mpars->ntix)[(mpars->mssr_order)[nn]] + 1 : 0;

    /* local optimization stage */
    for (R_len_t nn = 0; nn < (mpars->n); nn++)
    {
        if ((mpars->ntix)[nn] >= (mpars->s) && !R_IsNA(REAL_ELT(mssr, nn)))
        {
            (mpars->ntix)[nn] = 0;
            if ((mpars->nsp) == 0 || REAL_ELT(mssr, nn) < (1 + (mpars->tol)) * (mpars->mssropt))
            {
                gsl_matrix_get_row(pars->mp, pars->mx, nn);
                if (!Rf_isNull(pars->swts))
                    gsl_multifit_nlinear_winit(pars->mp, pars->wts, fdff, pars->w);
                else
                    gsl_multifit_nlinear_init(pars->mp, fdff, pars->w);
                mchisq1 = REAL_ELT(mssr, nn);
                gsl_multifit_nlinear_driver2(mpars->niter, xtol, 1e-3, ftol, NULL, NULL, &minfo, &mchisq0, &mchisq1, pars->w);
                gsl_multifit_nlinear_rcond(&rcond, pars->w);
                if (mchisq1 < (mpars->mssropt) && ((!gsl_isnan(rcond) && rcond > gtol) || mchisq1 < 2 * ftol))
                {
                    mpars->mssropt = mchisq1;
                    mpars->nsp += 1;
                    mpars->nwsp = 0;
                    gsl_vector_memcpy(pars->mpopt, (pars->w)->x);
                    gsl_vector_memcpy(pars->diag, trust_state->diag);
                    if (verbose)
                    {
                        size_t p = pars->mp->size;
                        Rprintf("mstart ssr* = %g, cond(J) = %g, NSP = %d, NWSP = %d, par = (", mpars->mssropt, 1.0 / rcond, mpars->nsp, mpars->nwsp);
                        for (R_len_t k = 0; k < p; k++)
                            Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", gsl_vector_get((pars->w)->x, k));
                    }
                }
                else
                    mpars->nwsp += 1;
            }
            else
                mpars->nwsp += 1;
        }
    }
}

/* deprecated function */
// int gsl_evalf(SEXP par, fdata *params, gsl_vector *f)
// {
//     /* evaluate function f */
//     SETCADR(params->f, par);
//     SEXP fval = PROTECT(Rf_eval(params->f, params->rho));

//     /* function checks */
//     R_len_t n = params->n;
//     if (TYPEOF(fval) != REALSXP || Rf_length(fval) != n)
//     {
//         UNPROTECT(1);
//         return GSL_EBADFUNC;
//     }

//     /* set gsl residuals */
//     double *fvalptr = REAL(fval);
//     double *yptr = REAL(params->y);
//     for (R_len_t i = 0; i < n; i++)
//     {
//         if (R_IsNaN(fvalptr[i]) || !R_finite(fvalptr[i]))
//             gsl_vector_set(f, i, GSL_POSINF);
//         else
//             gsl_vector_set(f, i, fvalptr[i] - yptr[i]);
//     }

//     UNPROTECT(1);
//     return GSL_SUCCESS;
// }