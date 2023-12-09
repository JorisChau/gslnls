#define R_NO_REMAP

#include <string.h>
#include "gsl_nls.h"

/* helper to match sparse matrix class */
static int match_dg_class(SEXP obj)
{
    int matclass = -2;
    SEXP klass = PROTECT(Rf_getAttrib(obj, R_ClassSymbol));
    R_len_t n = Rf_length(klass);
    if (n > 0)
    {
        for (R_len_t i = 0; i < n; i++)
        {
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgTMatrix") == 0)
            {
                matclass = 0;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgCMatrix") == 0)
            {
                matclass = 1;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgRMatrix") == 0)
            {
                matclass = 2;
                break;
            }
            if (strcmp(CHAR(STRING_ELT(klass, i)), "dgeMatrix") == 0)
            {
                matclass = -1;
                break;
            }
        }
    }
    UNPROTECT(1);
    return matclass;
}

/* cleanup memory */
static void C_nls_large_cleanup(void *data)
{
    pdata_large *pars = data;

    /* free memory */
    if (pars->w)
        gsl_multilarge_nlinear_free(pars->w);
    if (pars->J)
        gsl_matrix_free(pars->J);
    if (pars->Jsp)
        gsl_spmatrix_free(pars->Jsp);
}

/* function call w/ cleanup */
SEXP C_nls_large(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl)
{
    /* function arguments */
    pdata_large pars = {fn, y, jac, fvv, env, start, swts, control_int, control_dbl, NULL, NULL, NULL};

    /* safe function call */
    SEXP ans = R_ExecWithCleanup(C_nls_large_internal, &pars, C_nls_large_cleanup, &pars);

    return ans;
}

SEXP C_nls_large_internal(void *data)
{
    /* turn off error handler to avoid aborting on errors */
    gsl_set_error_handler_off();

    /* function arguments */
    pdata_large *pars = data;

    /* initialize parameters */
    int nprotect = 1;
    SEXP startvec = PROTECT(Rf_coerceVector(pars->start, REALSXP));
    R_len_t p = Rf_length(startvec);
    R_len_t n = Rf_length(pars->y);
    R_len_t niter = INTEGER_ELT(pars->control_int, 0);
    int verbose = INTEGER_ELT(pars->control_int, 1);

    /* control parameters */
    gsl_multilarge_nlinear_parameters fdf_params = gsl_multilarge_nlinear_default_parameters();

    /* minimization algorithm */
    switch (INTEGER_ELT(pars->control_int, 2))
    {
    case 1:
        fdf_params.trs = gsl_multilarge_nlinear_trs_lmaccel;
        break;
    case 2:
        fdf_params.trs = gsl_multilarge_nlinear_trs_dogleg;
        break;
    case 3:
        fdf_params.trs = gsl_multilarge_nlinear_trs_ddogleg;
        break;
    case 4:
        fdf_params.trs = gsl_multilarge_nlinear_trs_subspace2D;
        break;
    case 5:
        fdf_params.trs = gsl_multilarge_nlinear_trs_cgst;
        break;
    default:
        fdf_params.trs = gsl_multilarge_nlinear_trs_lm;
    }

    /* scaling method */
    switch (INTEGER_ELT(pars->control_int, 3))
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
    fdf_params.fdtype = INTEGER_ELT(pars->control_int, 4) ? GSL_MULTILARGE_NLINEAR_CTRDIFF : GSL_MULTILARGE_NLINEAR_FWDIFF;

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
    SEXP dfcall = PROTECT(Rf_lang2(pars->jac, xpar));
    SEXP parnames = PROTECT(Rf_getAttrib(pars->start, R_NamesSymbol));
    nprotect += 3;

    fdata_large params;
    params.n = n;
    params.p = p;
    params.f = fcall;
    params.df = dfcall;
    params.fvv = NULL;
    params.y = pars->y;
    params.rho = pars->env;
    params.start = pars->start;
    params.matclass = INTEGER_ELT(pars->control_int, 5);
    params.J = NULL;
    params.Jsp = NULL;

    /* allocate matrix to store Jacobian */
    if (params.matclass < 0)
    {
        pars->J = gsl_matrix_alloc(n, p);
        params.J = pars->J;
    }
    else
    {
        R_len_t nzmax = INTEGER_ELT(pars->control_int, 6);
        pars->Jsp = gsl_spmatrix_alloc_nzmax(n, p, nzmax, GSL_SPMATRIX_TRIPLET);
        params.Jsp = pars->Jsp;
    }

    if (verbose)
    {
        params.partrace = PROTECT(Rf_allocMatrix(REALSXP, niter + 1, p));
        params.ssrtrace = PROTECT(Rf_allocVector(REALSXP, niter + 1));
        nprotect += 2;
    }

    /* define the function to be minimized */
    gsl_multilarge_nlinear_fdf fdf;
    fdf.f = gsl_f;
    fdf.df = gsl_df_large;
    fdf.fvv = NULL; // not using geodesic acceleration
    fdf.n = n;
    fdf.p = p;
    fdf.params = &params;

    /* use acceleration function */
    if (!Rf_isNull(pars->fvv))
    {
        SEXP vpar = Rf_install("v");
        SEXP fvvcall = PROTECT(Rf_lang3(pars->fvv, xpar, vpar));
        params.fvv = fvvcall;
        fdf.fvv = gsl_fvv;
        nprotect++;
    }

    /* set nls optimization parameters */
    double *start1 = (double *)S_alloc(p, sizeof(double));
    for (R_len_t k = 0; k < p; k++)
        start1[k] = REAL_ELT(startvec, k);
    gsl_vector_view par = gsl_vector_view_array(start1, p);

    /* initialize solver */
    const gsl_multilarge_nlinear_type *T = gsl_multilarge_nlinear_trust;

    /* allocate workspace with default parameters */
    pars->w = gsl_multilarge_nlinear_alloc(T, &fdf_params, n, p);

    /* initialize solver with starting point and weights */
    if (!Rf_isNull(pars->swts))
    {
        double *swts1 = (double *)S_alloc(n, sizeof(double));
        for (R_len_t i = 0; i < n; i++)
            swts1[i] = REAL_ELT(pars->swts, i);
        gsl_vector_view wts = gsl_vector_view_array(swts1, n);
        gsl_multilarge_nlinear_winit(&par.vector, &wts.vector, &fdf, pars->w);
    }
    else
    {
        gsl_multilarge_nlinear_init(&par.vector, &fdf, pars->w);
    }

    /* compute initial cost function */
    double chisq_init = GSL_POSINF;
    gsl_vector *resid = gsl_multilarge_nlinear_residual(pars->w);
    gsl_blas_ddot(resid, resid, &chisq_init);
    double chisq0 = chisq_init;
    double chisq1 = chisq_init;
    params.chisq = chisq_init;

    if (verbose)
    {
        SET_REAL_ELT(params.ssrtrace, 0, chisq_init);
        double *parptr = REAL(params.partrace);
        for (R_len_t k = 0; k < p; k++)
            parptr[(niter + 1) * k] = start1[k];
    }

    /* solve the system  */
    int info = GSL_CONTINUE;
    int status = gsl_multilarge_nlinear_driver2(niter, xtol, gtol, ftol, verbose ? callback_large : NULL, verbose ? &params : NULL, &info, &chisq0, &chisq1, pars->w);
    R_len_t iter = gsl_multilarge_nlinear_niter(pars->w);

    /* compute covariance and cost at best fit parameters */
    gsl_matrix *cov = NULL;
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        cov = gsl_matrix_alloc(p, p);
        gsl_multilarge_nlinear_covar(cov, pars->w);
    }

    if (verbose)
    {
        /* print summary statistics*/
        Rprintf("*******************\nsummary from method 'multilarge/%s'\n", gsl_multilarge_nlinear_trs_name(pars->w));
        Rprintf("number of iterations: %d\n", iter);
        Rprintf("reason for stopping: %s\n", gsl_strerror(info));
        Rprintf("initial ssr = %g\n", chisq_init);
        Rprintf("final ssr = %g\n", chisq1);
        Rprintf("ssr/dof = %g\n", chisq1 / (n - p));
        Rprintf("ssr achieved tolerance = %g\n", chisq0 - chisq1);
        Rprintf("function evaluations: %d\n", (int)fdf.nevalf);
        Rprintf("jacobian-vector product evaluations: %d\n", (int)fdf.nevaldfu);
        Rprintf("jacobian-jacobian product evaluations: %d\n", (int)fdf.nevaldf2);
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
            SET_REAL_ELT(anspar, k, gsl_vector_get(pars->w->x, k));
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
        if (params.matclass < 0) // dense Jacobian
        {
            for (R_len_t i = 0; i < n; i++)
                for (R_len_t k = 0; k < p; k++)
                    jacptr[i + n * k] = gsl_matrix_get(params.J, i, k);
        }
        else // sparse Jacobian
        {
            for (R_len_t i = 0; i < n; i++)
                for (R_len_t k = 0; k < p; k++)
                    jacptr[i + n * k] = gsl_spmatrix_get(params.Jsp, i, k);
        }
    }
    else
    {
        for (R_len_t i = 0; i < n; i++)
            for (R_len_t k = 0; k < p; k++)
                jacptr[i + n * k] = NA_REAL;
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
    SET_VECTOR_ELT(ans, 9, Rf_ScalarString(Rf_mkChar(gsl_multilarge_nlinear_trs_name(pars->w))));

    const char *nms[] = {"f", "dfu", "df2", "fvv", ""};
    SEXP ansneval = PROTECT(Rf_mkNamed(INTSXP, nms));
    SET_INTEGER_ELT(ansneval, 0, fdf.nevalf);
    SET_INTEGER_ELT(ansneval, 1, fdf.nevaldfu);
    SET_INTEGER_ELT(ansneval, 2, fdf.nevaldf2);
    SET_INTEGER_ELT(ansneval, 3, fdf.nevalfvv);
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
} /* gsl_multilarge_nlinear_driver() */

int gsl_df_large(CBLAS_TRANSPOSE_t TransJ, const gsl_vector *x, const gsl_vector *u, void *params, gsl_vector *v, gsl_matrix *JTJ)
{
    /* construct parameter vector */
    int nprotect = 2;
    SEXP par = NULL;
    R_len_t p = ((fdata_large *)params)->p;
    R_len_t n = ((fdata_large *)params)->n;
    SEXP start = ((fdata_large *)params)->start;
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
    SETCADR(((fdata_large *)params)->df, par);
    SEXP dfval = PROTECT(Rf_eval(((fdata_large *)params)->df, ((fdata_large *)params)->rho));

    /* populate Jacobian */
    int jacclass = ((fdata_large *)params)->matclass;
    int dgclass = -2;

    if (jacclass == -2) // dense numeric matrix
    {
        /* Jacobian checks */
        if (TYPEOF(dfval) != REALSXP || !Rf_isMatrix(dfval) || Rf_ncols(dfval) != p || Rf_nrows(dfval) != n)
        {
            Rf_warning("Evaluating jac does not return numeric matrix of size n x p");
            UNPROTECT(nprotect);
            return GSL_EBADFUNC;
        }

        double *jacptr = REAL(dfval);
        for (R_len_t i = 0; i < n; i++)
            for (R_len_t k = 0; k < p; k++)
                if (R_IsNaN(jacptr[i + n * k]) || !R_finite(jacptr[i + n * k]))
                {
                    Rf_warning("Missing/infinite values not allowed when evaluating jac");
                    UNPROTECT(nprotect);
                    return GSL_EBADFUNC;
                }

        for (R_len_t i = 0; i < n; i++)
            for (R_len_t k = 0; k < p; k++)
                gsl_matrix_set(((fdata_large *)params)->J, i, k, jacptr[i + n * k]);
    }
    else
    {
        if (!Rf_isS4(dfval) || !R_has_slot(dfval, Rf_install("x")) || !R_has_slot(dfval, Rf_install("Dim")))
        {
            Rf_warning("Evaluating jac does not return a \"dgCMatrix\", \"dgRMatrix\", \"dgTMatrix\", or \"dgeMatrix\"");
            UNPROTECT(nprotect);
            return GSL_EBADFUNC;
        }

        SEXP dfdim = PROTECT(R_do_slot(dfval, Rf_install("Dim")));
        nprotect++;

        if (INTEGER_ELT(dfdim, 0) != n || INTEGER_ELT(dfdim, 1) != p)
        {
            Rf_warning("Evaluating jac does not return matrix of size n x p");
            UNPROTECT(nprotect);
            return GSL_EBADFUNC;
        }

        SEXP dfx = PROTECT(R_do_slot(dfval, Rf_install("x")));
        R_len_t nz = Rf_length(dfx);
        double *xptr = REAL(dfx);
        nprotect++;

        for (R_len_t m = 0; m < nz; m++)
            if (R_IsNaN(xptr[m]) || !R_finite(xptr[m]))
            {
                Rf_warning("Missing/infinite values not allowed when evaluating jac");
                UNPROTECT(nprotect);
                return GSL_EBADFUNC;
            }

        /* re-detect sparse matrix class */
        dgclass = match_dg_class(dfval);

        if (dgclass == -2)
        {
            Rf_warning("Evaluating jac does not return a \"dgCMatrix\", \"dgRMatrix\", \"dgTMatrix\", or \"dgeMatrix\"");
            UNPROTECT(nprotect);
            return GSL_EBADFUNC;
        }
        else if (dgclass == -1) // dgeMatrix
        {
            for (R_len_t i = 0; i < n; i++)
                for (R_len_t k = 0; k < p; k++)
                    gsl_matrix_set(((fdata_large *)params)->J, i, k, xptr[i + n * k]);
        }
        else if (dgclass == 0) // dgTMatrix
        {
            SEXP dfi = PROTECT(R_do_slot(dfval, Rf_install("i")));
            SEXP dfj = PROTECT(R_do_slot(dfval, Rf_install("j")));
            R_len_t *iptr = INTEGER(dfi);
            R_len_t *jptr = INTEGER(dfj);
            gsl_spmatrix_set_zero(((fdata_large *)params)->Jsp);

            for (R_len_t m = 0; m < nz; m++)
            {
                gsl_spmatrix_set(((fdata_large *)params)->Jsp, iptr[m], jptr[m], xptr[m]);
            }
            nprotect += 2;
        }
        else if (dgclass == 1) // dgCMatrix
        {
            SEXP dfi = PROTECT(R_do_slot(dfval, Rf_install("i")));
            SEXP dfp = PROTECT(R_do_slot(dfval, Rf_install("p")));
            R_len_t *iptr = INTEGER(dfi);
            R_len_t *pptr = INTEGER(dfp);
            gsl_spmatrix_set_zero(((fdata_large *)params)->Jsp);

            R_len_t j = 0;
            for (R_len_t m = 0; m < nz; m++)
            {
                while (pptr[j] <= m)
                    j += 1;
                gsl_spmatrix_set(((fdata_large *)params)->Jsp, iptr[m], j - 1, xptr[m]);
            }
            nprotect += 2;
        }
        else if (dgclass == 2) // dgRMatrix
        {
            SEXP dfp = PROTECT(R_do_slot(dfval, Rf_install("p")));
            SEXP dfj = PROTECT(R_do_slot(dfval, Rf_install("j")));
            R_len_t *pptr = INTEGER(dfp);
            R_len_t *jptr = INTEGER(dfj);
            gsl_spmatrix_set_zero(((fdata_large *)params)->Jsp);

            R_len_t i = 0;
            for (R_len_t m = 0; m < nz; m++)
            {
                while (pptr[i] <= m)
                    i += 1;
                gsl_spmatrix_set(((fdata_large *)params)->Jsp, i - 1, jptr[m], xptr[m]);
            }
            nprotect += 2;
        }
    }

    if (jacclass == -2 || dgclass == -1) // non-sparse jacobian
    {
        /* calculate J * u or J' * u and return in v */
        if (v)
            gsl_blas_dgemv(TransJ, 1.0, ((fdata_large *)params)->J, u, 0.0, v);

        /* calculate J'J and return in JTJ */
        if (JTJ)
            gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, ((fdata_large *)params)->J, 0.0, JTJ);
    }
    else // sparse jacobian
    {
        /* calculate J * u or J' * u and return in v */
        if (v)
            gsl_spblas_dgemv(TransJ, 1.0, ((fdata_large *)params)->Jsp, u, 0.0, v);

        /* calculate J'J and return in JTJ */
        if (JTJ)
        {
            gsl_matrix *J = gsl_matrix_alloc(n, p);
            gsl_spmatrix_sp2d(J, ((fdata_large *)params)->Jsp);
            gsl_blas_dsyrk(CblasLower, CblasTrans, 1.0, J, 0.0, JTJ);
            gsl_matrix_free(J);
        }
    }

    UNPROTECT(nprotect);
    return GSL_SUCCESS;
}

void callback_large(const R_len_t iter, void *params, const gsl_multilarge_nlinear_workspace *w)
{
    /* ssr trace */
    double chisq = ((fdata_large *)params)->chisq;
    SET_REAL_ELT(((fdata_large *)params)->ssrtrace, iter, chisq);

    /* parameter trace */
    R_len_t p = ((fdata_large *)params)->p;
    R_len_t n = (R_len_t)Rf_nrows(((fdata_large *)params)->partrace);
    double *parptr = REAL(((fdata_large *)params)->partrace);
    gsl_vector *x = gsl_multilarge_nlinear_position(w);
    for (R_len_t k = 0; k < p; k++)
        parptr[iter + n * k] = gsl_vector_get(x, k);

    /* parameter norm */
    double xsq;
    gsl_blas_ddot(x, x, &xsq);

    /* cond(J) */
    double rcond;
    gsl_multilarge_nlinear_rcond(&rcond, w);

    /* print trace */
    Rprintf("iter %3d: ssr = %g, |x|^2 = %g, cond(J) = %g\n", iter, chisq, xsq, 1.0 / rcond);
}
