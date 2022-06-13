#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

/* helper to evaluate f at par */
static int gsl_evalf(SEXP par, fdata *params, gsl_vector *f)
{
    /* evaluate function f */
    SETCADR(params->f, par);
    SEXP fval = PROTECT(Rf_eval(params->f, params->rho));

    /* function checks */
    R_len_t n = params->n;
    if (TYPEOF(fval) != REALSXP || Rf_length(fval) != n)
    {
        UNPROTECT(1);
        return GSL_EBADFUNC;
    }

    /* set gsl residuals */
    double *fvalptr = REAL(fval);
    double *yptr = REAL(params->y);
    for (R_len_t i = 0; i < n; i++)
    {
        if (R_IsNaN(fvalptr[i]) || !R_finite(fvalptr[i]))
            gsl_vector_set(f, i, GSL_POSINF);
        else
            gsl_vector_set(f, i, fvalptr[i] - yptr[i]);
    }

    UNPROTECT(1);
    return GSL_SUCCESS;
}

SEXP C_nls_multistart(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP control_int, SEXP control_dbl)
{
    /* turn off error handler to avoid aborting on errors */
    gsl_set_error_handler_off();

    /* initialize parameters */
    int nprotect = 1;
    R_len_t p = Rf_nrows(start);
    R_len_t n = Rf_length(y);
    R_len_t niter = INTEGER_ELT(control_int, 0);
    int verbose = INTEGER_ELT(control_int, 1);

    /* control parameters */
    gsl_multifit_nlinear_parameters fdf_params = gsl_multifit_nlinear_default_parameters();

    /* minimization algorithm */
    switch (INTEGER_ELT(control_int, 2))
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
    switch (INTEGER_ELT(control_int, 3))
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
    switch (INTEGER_ELT(control_int, 4))
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
    fdf_params.fdtype = INTEGER_ELT(control_int, 5) ? GSL_MULTIFIT_NLINEAR_CTRDIFF : GSL_MULTIFIT_NLINEAR_FWDIFF;

    /* miscellaneous parameters */
    fdf_params.factor_up = REAL_ELT(control_dbl, 0);
    fdf_params.factor_down = REAL_ELT(control_dbl, 1);
    fdf_params.avmax = REAL_ELT(control_dbl, 2);
    fdf_params.h_df = REAL_ELT(control_dbl, 3);
    fdf_params.h_fvv = REAL_ELT(control_dbl, 4);
    double xtol = REAL_ELT(control_dbl, 5);
    double ftol = REAL_ELT(control_dbl, 6);
    double gtol = REAL_ELT(control_dbl, 7);
    R_len_t ntest = INTEGER_ELT(control_int, 6);
    R_len_t nseed = INTEGER_ELT(control_int, 7);

    /* initialize data */
    SEXP xpar = Rf_install("par");
    SEXP fcall = PROTECT(Rf_lang2(fn, xpar));
    SEXP parnames = PROTECT(Rf_GetRowNames(Rf_getAttrib(start, R_DimNamesSymbol)));
    nprotect += 2;

    fdata params;
    params.n = n;
    params.p = p;
    params.f = fcall;
    params.df = NULL;
    params.fvv = NULL;
    params.y = y;
    params.rho = env;
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
    if (!Rf_isNull(jac))
    {
        SEXP dfcall = PROTECT(Rf_lang2(jac, xpar));
        params.df = dfcall;
        fdf.df = gsl_df;
        nprotect++;
    }

    /* use acceleration function */
    if (!Rf_isNull(fvv))
    {
        SEXP vpar = Rf_install("v");
        SEXP fvvcall = PROTECT(Rf_lang3(fvv, xpar, vpar));
        params.fvv = fvvcall;
        fdf.fvv = gsl_fvv;
        nprotect++;
    }

    /* use weights */
    gsl_vector *wts = gsl_vector_alloc(n);
    if (!Rf_isNull(swts))
        for (R_len_t i = 0; i < n; i++)
            gsl_vector_set(wts, i, REAL_ELT(swts, i));
    else 
        gsl_vector_set_all(wts, 1.0);

    /* initialize solver */
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

    /* allocate workspace with default parameters */
    gsl_multifit_nlinear_workspace *w;
    w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);

    /* allocate qrng workspace */
    gsl_qrng *q;
    if (p < 41)
        q = gsl_qrng_alloc(gsl_qrng_sobol, p);
    else if (p < 1230)
        q = gsl_qrng_alloc(gsl_qrng_halton, p);
    else
        Rf_error("GSL quasi-random Halton sequences are only available up to 1229 parameters"); /* stop with error if p too large */

    /* initialize multistart parameters */
    int statustest;
    double chisqtest, pkval;
    gsl_matrix *partestmat = gsl_matrix_alloc(ntest, p);
    double *partestarray = (double *)R_alloc(p, sizeof(double));
    gsl_vector *ftest = gsl_vector_alloc(n);

    double *startptr = REAL(start);
    SEXP partest = PROTECT(Rf_allocVector(REALSXP, p));
    Rf_setAttrib(partest, R_NamesSymbol, parnames);
    SEXP chisqtestvec = PROTECT(Rf_allocVector(REALSXP, ntest));
    nprotect += 2;

    /* global stage: evaluate f at each test parameter vector */
    for (R_len_t i = 0; i < ntest; i++)
    {
        /* get parameters */
        gsl_qrng_get(q, partestarray);
        /* evaluate f at parameters*/
        for (R_len_t k = 0; k < p; k++)
        {
            pkval = startptr[k] + (startptr[p + k] - startptr[k]) * partestarray[k];
            gsl_matrix_set(partestmat, i, k, pkval);
            SET_REAL_ELT(partest, k, pkval);
        }
        statustest = gsl_evalf(partest, &params, ftest);
        /* evaluate ssr */
        if (statustest == GSL_SUCCESS)
        {
            gsl_blas_ddot(ftest, ftest, &chisqtest);
            SET_REAL_ELT(chisqtestvec, i, chisqtest);
        }
        else
            SET_REAL_ELT(chisqtestvec, i, NA_REAL);

        /* for testing */
        Rprintf("global iter: %d, ssr: %g\n", i, chisqtest);
    }

    /* order ssr vector */
    int *chisqidx = (int *)R_alloc(ntest, sizeof(int));
    R_orderVector1(chisqidx, ntest, chisqtestvec, TRUE, FALSE);

    /* local stage: optimize ssr at each seed parameter vector */
    int seed, seedinfo, seedstatus;
    double theta, chisqseed0, chisqseed = GSL_POSINF, chisqopt = GSL_POSINF;
    gsl_vector *seedpar = gsl_vector_alloc(p);
    gsl_vector *seedparopt = gsl_vector_alloc(p);
    gsl_matrix_get_row(seedparopt, partestmat, chisqidx[0]);

    /* update fdf params */
    params.start = partest;

    for (R_len_t i = 0; i < nseed; i++)
    {
        /* get seed parameters */
        seed = chisqidx[i];
        gsl_matrix_get_row(seedpar, partestmat, seed);

        /* tiktak convex combination */
        theta = sqrt(((double)i) / ((double)nseed));
        theta = theta < 0.995 ? (theta > 0.1 ? theta : 0.1) : 0.995;
        for (R_len_t k = 0; k < p; k++)
            gsl_vector_set(seedpar, k, (1.0 - theta) * gsl_vector_get(seedpar, k) + theta * gsl_vector_get(seedparopt, k));

        /* re-initialize solver */
        if (!Rf_isNull(swts))
            gsl_multifit_nlinear_winit(seedpar, wts, &fdf, w);
        else
            gsl_multifit_nlinear_init(seedpar, &fdf, w);

        /* optimize at seed parameters*/
        chisqseed = REAL_ELT(chisqtestvec, seed);
        if (R_IsNA(chisqseed))
            chisqseed = GSL_POSINF;

        seedstatus = gsl_multifit_nlinear_driver2(niter / 2, xtol, 1e-3, ftol, NULL, NULL, &seedinfo, &chisqseed0, &chisqseed, w);

        /* update optimal parameters */
        if ((seedstatus == GSL_SUCCESS || seedstatus == GSL_EMAXITER) && chisqseed < chisqopt)
        {
            chisqopt = chisqseed;
            gsl_vector_memcpy(seedparopt, w->x);
        }

        /* for testing */
        Rprintf("local iter: %d, status: %d, ssr: %g, ssr opt: %g, p0: %g, p1: %g\n", i, seedstatus, chisqseed, chisqopt, gsl_vector_get(w->x, 0), gsl_vector_get(w->x, 1));
    }

    /* re-execute local search at optimal parameters */
    /* re-initialize solver */
    if (!Rf_isNull(swts))
        gsl_multifit_nlinear_winit(seedparopt, wts, &fdf, w);
    else
        gsl_multifit_nlinear_init(seedparopt, &fdf, w);

    /* compute initial cost function */
    double chisq_init = GSL_POSINF;
    gsl_vector *resid = gsl_multifit_nlinear_residual(w);
    gsl_blas_ddot(resid, resid, &chisq_init);
    double chisq0 = chisq_init;
    double chisq1 = chisq_init;
    params.chisq = chisq_init;

    if (verbose)
    {
        SET_REAL_ELT(params.ssrtrace, 0, chisq_init);
        double *parptr = REAL(params.partrace);
        for (R_len_t k = 0; k < p; k++)
            parptr[(niter + 1) * k] = gsl_vector_get(seedparopt, k);
    }

    /* solve the system  */
    int info = GSL_CONTINUE;
    int status = gsl_multifit_nlinear_driver2(niter, xtol, gtol, ftol, verbose ? callback : NULL, verbose ? &params : NULL, &info, &chisq0, &chisq1, w);
    R_len_t iter = gsl_multifit_nlinear_niter(w);

    /* compute covariance and cost at best fit parameters */
    gsl_matrix *J = NULL;
    gsl_matrix *cov = NULL;
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
    {
        cov = gsl_matrix_alloc(p, p);
        J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, cov);
    }

    if (verbose)
    {
        /* print summary statistics*/
        Rprintf("*******************\nsummary from method 'multifit/%s'\n", gsl_multifit_nlinear_trs_name(w));
        Rprintf("number of iterations: %d\n", iter);
        Rprintf("reason for stopping: %s\n", gsl_strerror(info));
        Rprintf("initial ssr = %g\n", chisq_init);
        Rprintf("final ssr = %g\n", chisq1);
        Rprintf("ssr/dof = %g\n", chisq1 / (n - p));
        Rprintf("ssr achieved tolerance = %g\n", chisq0 - chisq1);
        Rprintf("function evaluations: %d\n", fdf.nevalf);
        Rprintf("jacobian evaluations: %d\n", fdf.nevaldf);
        Rprintf("fvv evaluations: %d\n", fdf.nevalfvv);
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
            SET_REAL_ELT(anspar, k, gsl_vector_get(w->x, k));
    }
    else
    {
        for (R_len_t k = 0; k < p; k++)
            SET_REAL_ELT(anspar, k, gsl_vector_get(seedparopt, k));
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
    SET_VECTOR_ELT(ans, 9, Rf_ScalarString(Rf_mkChar(gsl_multifit_nlinear_trs_name(w))));

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
    gsl_qrng_free(q);
    gsl_vector_free(ftest);
    gsl_vector_free(seedpar);
    gsl_vector_free(wts);
    gsl_vector_free(seedparopt);
    gsl_matrix_free(partestmat);

    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
        gsl_matrix_free(cov);

    UNPROTECT(nprotect);
    return ans;
}
