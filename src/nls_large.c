#define R_NO_REMAP

#include "gsl_nls.h"

SEXP C_nls_large(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl)
{
    /* turn off error handler to avoid aborting on errors */
    gsl_set_error_handler_off();

    /* initialize parameters */
    int nprotect = 1;
    SEXP startvec = PROTECT(Rf_coerceVector(start, REALSXP));
    R_len_t p = Rf_length(startvec);
    R_len_t n = Rf_length(y);
    R_len_t niter = INTEGER_ELT(control_int, 0);
    int verbose = INTEGER_ELT(control_int, 1);

    /* control parameters */
    gsl_multilarge_nlinear_parameters fdf_params = gsl_multilarge_nlinear_default_parameters();

    /* minimization algorithm */
    switch (INTEGER_ELT(control_int, 2))
    {
    case 1:
        fdf_params.trs = gsl_multilarge_nlinear_trs_lm;
    case 2:
        fdf_params.trs = gsl_multilarge_nlinear_trs_lmaccel;
        break;
    case 3:
        fdf_params.trs = gsl_multilarge_nlinear_trs_dogleg;
        break;
    case 4:
        fdf_params.trs = gsl_multilarge_nlinear_trs_ddogleg;
        break;
    case 5:
        fdf_params.trs = gsl_multilarge_nlinear_trs_subspace2D;
        break;
    default:
        fdf_params.trs = gsl_multilarge_nlinear_trs_cgst;
    }

    /* scaling method */
    switch (INTEGER_ELT(control_int, 3))
    {
    case 1:
        fdf_params.scale = gsl_multilarge_nlinear_scale_levenberg;
        break;
    case 2:
        fdf_params.scale = gsl_multilarge_nlinear_scale_marquardt;
        break;
    default:
        fdf_params.scale = gsl_multilarge_nlinear_scale_more;
    }

    /* finite differencing type */
    fdf_params.fdtype = INTEGER_ELT(control_int, 5) ? GSL_MULTILARGE_NLINEAR_CTRDIFF : GSL_MULTILARGE_NLINEAR_FWDIFF;

    /* miscellaneous parameters */
    fdf_params.factor_up = REAL_ELT(control_dbl, 0);
    fdf_params.factor_down = REAL_ELT(control_dbl, 1);
    fdf_params.avmax = REAL_ELT(control_dbl, 2);
    fdf_params.h_df = REAL_ELT(control_dbl, 3);
    fdf_params.h_fvv = REAL_ELT(control_dbl, 4);
    double xtol = REAL_ELT(control_dbl, 5);
    double ftol = REAL_ELT(control_dbl, 6);
    double gtol = REAL_ELT(control_dbl, 7);

    /* initialize data */
    SEXP xpar = Rf_install("par");
    SEXP fcall = PROTECT(Rf_lang2(fn, xpar));
    SEXP parnames = PROTECT(Rf_getAttrib(start, R_NamesSymbol));
    nprotect += 2;

    fdata params;
    params.n = n;
    params.p = p;
    params.f = fcall;
    params.df = NULL;
    params.fvv = NULL;
    params.y = y;
    params.rho = env;
    params.start = start;

    if (verbose)
    {
        params.partrace = PROTECT(Rf_allocMatrix(REALSXP, niter + 1, p));
        params.ssrtrace = PROTECT(Rf_allocVector(REALSXP, niter + 1));
        nprotect += 2;
    }

    /* define the function to be minimized */
    gsl_multilarge_nlinear_fdf fdf;
    fdf.f = gsl_f;
    fdf.df = NULL;  // set to NULL for finite-difference Jacobian
    fdf.fvv = NULL; // not using geodesic acceleration
    fdf.n = n;
    fdf.p = p;
    fdf.nevalf = 0;
    fdf.nevaldfu = 0;
    fdf.nevaldf2 = 0;
    fdf.nevalfvv = 0;
    fdf.params = &params;

    /* use Jacobian function */
    // if (!Rf_isNull(jac))
    // {
    //     SEXP dfcall = PROTECT(Rf_lang2(jac, xpar));
    //     params.df = dfcall;
    //     fdf.df = gsl_df;
    //     nprotect++;
    // }

    /* use acceleration function */
    // if (!Rf_isNull(fvv))
    // {
    //     SEXP vpar = Rf_install("v");
    //     SEXP fvvcall = PROTECT(Rf_lang3(fvv, xpar, vpar));
    //     params.fvv = fvvcall;
    //     fdf.fvv = gsl_fvv;
    //     nprotect++;
    // }

    /* set nls optimization parameters */
    double *start1 = (double *)S_alloc(p, sizeof(double));
    for (R_len_t k = 0; k < p; k++)
        start1[k] = REAL_ELT(startvec, k);
    gsl_vector_view par = gsl_vector_view_array(start1, p);

    /* initialize solver */
    const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;
    gsl_multilarge_nlinear_workspace *w;

    /* allocate workspace with default parameters */
    w = gsl_multilarge_nlinear_alloc(T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    gsl_multilarge_nlinear_init(&par.vector, &fdf, w);

    // if (!Rf_isNull(swts))
    // {
    //     double *swts1 = (double *)S_alloc(n, sizeof(double));
    //     for (R_len_t i = 0; i < n; i++)
    //         swts1[i] = REAL_ELT(swts, i);
    //     gsl_vector_view wts = gsl_vector_view_array(swts1, n);
    //     gsl_multifit_nlinear_winit(&par.vector, &wts.vector, &fdf, w);
    // }
    // else
    // {
    //     gsl_multifit_nlinear_init(&par.vector, &fdf, w);
    // }

    /* compute initial cost function */
    double chisq_init = 0.0;
    gsl_vector *resid = gsl_multilarge_nlinear_residual(w);
    gsl_blas_ddot(resid, resid, &chisq_init);
    double chisq0 = chisq_init;
    double chisq1 = chisq_init;

    /* solve the system  */
    int info = GSL_CONTINUE;
    int status = gsl_multilarge_nlinear_driver2(niter, xtol, gtol, ftol, NULL, NULL, &info, w);
    R_len_t iter = gsl_multilarge_nlinear_niter(w);

    /* compute covariance and cost at best fit parameters */
    gsl_matrix *J = NULL;
    gsl_matrix *cov = NULL;
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        cov = gsl_matrix_alloc(p, p);
        gsl_multilarge_nlinear_covar(cov, w);
    }

    // if (verbose)
    // {
    //     /* print summary statistics*/
    //     Rprintf("*******************\nsummary from method '%s/%s'\n", gsl_multifit_nlinear_name(w), gsl_multifit_nlinear_trs_name(w));
    //     Rprintf("number of iterations: %d\n", iter);
    //     Rprintf("reason for stopping: %s\n", gsl_strerror(info));
    //     Rprintf("initial ssr = %g\n", chisq_init);
    //     Rprintf("final ssr = %g\n", chisq1);
    //     Rprintf("ssr/dof = %g\n", chisq1 / (n - p));
    //     Rprintf("ssr achieved tolerance = %g\n", chisq0 - chisq1);
    //     Rprintf("function evaluations: %d\n", fdf.nevalf);
    //     Rprintf("Jacobian evaluations: %d\n", fdf.nevaldf);
    //     Rprintf("fvv evaluations: %d\n", fdf.nevalfvv);
    //     Rprintf("status = %s\n*******************\n", gsl_strerror(status));
    // }

    /* initialize result */
    SEXP ans = NULL;
    // if (verbose)
    // {
    //     const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
    //                             "algorithm", "neval", "partrace", "ssrtrace", ""};
    //     ans = PROTECT(Rf_mkNamed(VECSXP, ansnms));
    // }
    // else
    // {
        const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
                                "algorithm", "neval", ""};
        ans = PROTECT(Rf_mkNamed(VECSXP, ansnms));
    // }
    nprotect++;

    /* estimated parameters */
    SEXP anspar = PROTECT(Rf_allocVector(REALSXP, p));
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(anspar, k, gsl_vector_get(w->x, k));
    }
    else
    {
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(anspar, k, REAL_ELT(startvec, k));
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
    SET_VECTOR_ELT(ans, 9, Rf_ScalarString(Rf_mkChar(gsl_multilarge_nlinear_trs_name(w))));

    const char *nms[] = {"f", "dfu", "df2", "fvv", ""};
    SEXP ansneval = PROTECT(Rf_mkNamed(INTSXP, nms));
    SET_INTEGER_ELT(ansneval, 0, fdf.nevalf);
    SET_INTEGER_ELT(ansneval, 1, fdf.nevaldfu);
    SET_INTEGER_ELT(ansneval, 2, fdf.nevaldf2);
    SET_INTEGER_ELT(ansneval, 3, fdf.nevalfvv);
    SET_VECTOR_ELT(ans, 10, ansneval);
    UNPROTECT(1);

    // if (verbose)
    // {
    //     if (!Rf_isNull(parnames))
    //     {
    //         SEXP dimnames = PROTECT(Rf_allocVector(VECSXP, 2));
    //         SET_VECTOR_ELT(dimnames, 0, R_NilValue);
    //         SET_VECTOR_ELT(dimnames, 1, parnames);
    //         Rf_setAttrib(params.partrace, R_DimNamesSymbol, dimnames);
    //         UNPROTECT(1);
    //     }
    //     SET_VECTOR_ELT(ans, 11, params.partrace);
    //     SET_VECTOR_ELT(ans, 12, params.ssrtrace);
    // }

    /* free memory */
    UNPROTECT(nprotect);
    gsl_multilarge_nlinear_free(w);
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
        gsl_matrix_free(cov);

    return ans;
}

/*
gsl_multilarge_nlinear_driver()
  Iterate the nonlinear least squares solver until completion

Inputs: maxiter  - maximum iterations to allow
        xtol     - tolerance in step x
        gtol     - tolerance in gradient
        ftol     - tolerance in ||f||
        callback - callback function to call each iteration
        callback_params - parameters to pass to callback function
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
GSL_MAXITER if maxiter exceeded without converging
GSL_ENOPROG if no accepted step found on first iteration
*/
int gsl_multilarge_nlinear_driver2(const size_t maxiter,
                                  const double xtol,
                                  const double gtol,
                                  const double ftol,
                                  void (*callback)(const size_t iter, void *params,
                                                   const gsl_multilarge_nlinear_workspace *w),
                                  void *callback_params,
                                  int *info,
                                  gsl_multilarge_nlinear_workspace *w)
{
    int status = GSL_CONTINUE;
    size_t iter = 0;

    /* call user callback function prior to any iterations
   * with initial system state */
    if (callback)
        callback(iter, callback_params, w);

    do
    {
        status = gsl_multilarge_nlinear_iterate(w);

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
} /* gsl_multilarge_nlinear_driver() */

// int gsl_df_large(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u, void *params, gsl_vector *v, gsl_matrix *JTJ);
// {
//     /* construct parameter vector */
//     SEXP par = NULL;
//     R_len_t p = ((fdata *)params)->p;
//     SEXP start = ((fdata *)params)->start;
//     if (Rf_isNumeric(start))
//     {
//         par = PROTECT(Rf_allocVector(REALSXP, p));
//         for (R_len_t k = 0; k < p; k++)
//             SET_REAL_ELT(par, k, gsl_vector_get(x, k));
//     }
//     else
//     {
//         par = PROTECT(Rf_allocVector(VECSXP, p));
//         for (R_len_t k = 0; k < p; k++)
//             SET_VECTOR_ELT(par, k, Rf_ScalarReal(gsl_vector_get(x, k)));
//     }
//     Rf_setAttrib(par, R_NamesSymbol, Rf_getAttrib(start, R_NamesSymbol));

//     /* evaluate Jacobian function */
//     SETCADR(((fdata *)params)->df, par);
//     SEXP dfval = PROTECT(Rf_eval(((fdata *)params)->df, ((fdata *)params)->rho));

//     /* Jacobian checks */
//     R_len_t n = ((fdata *)params)->n;
//     if (TYPEOF(dfval) != REALSXP || !Rf_isMatrix(dfval) || Rf_ncols(dfval) != p || Rf_nrows(dfval) != n)
//     {
//         Rf_warning("Evaluating jac does not return numeric matrix of dimensions n x p");
//         UNPROTECT(2);
//         return GSL_EBADFUNC;
//     }

//     double *jacptr = REAL(dfval);
//     for (R_len_t i = 0; i < n; i++)
//         for (R_len_t k = 0; k < p; k++)
//             if (R_IsNaN(jacptr[i + n * k]) || !R_finite(jacptr[i + n * k]))
//             {
//                 Rf_warning("Missing/infinite values not allowed when evaluating jac");
//                 UNPROTECT(2);
//                 return GSL_EBADFUNC;
//             }

//     /* set gsl jacobian matrix */
//     for (R_len_t i = 0; i < n; i++)
//         for (R_len_t k = 0; k < p; k++)
//             gsl_matrix_set(J, i, k, jacptr[i + n * k]);

//     UNPROTECT(2);
//     return GSL_SUCCESS;
// }

// void callback_large(const size_t iter, void *params, const gsl_multilarge_nlinear_workspace *w)
// {
//     gsl_vector *f = gsl_multilarge_nlinear_residual(w);
//     double chisq = 0.0;
//     gsl_blas_ddot(f, f, &chisq); // current ssr

//     /* update traces */
//     SET_REAL_ELT(((fdata *)params)->ssrtrace, (R_len_t)iter, chisq);
//     R_len_t p = ((fdata *)params)->p;
//     R_len_t n = (R_len_t)Rf_nrows(((fdata *)params)->partrace);
//     double *parptr = REAL(((fdata *)params)->partrace);
//     for (R_len_t k = 0; k < p; k++)
//         parptr[iter + n * k] = gsl_vector_get(w->x, k);

//     /* print trace */
//     Rprintf("iter %3d: ssr = %g, par = (", iter, chisq);
//     for (R_len_t k = 0; k < p; k++)
//         Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", parptr[iter + n * k]);
// }