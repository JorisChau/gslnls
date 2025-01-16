#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

/* static function declarations */
static int gsl_f(const gsl_vector *x, void *params, gsl_vector *f);

static int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J);

static int gsl_fvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv);

static void callback(const R_len_t iter, void *params, const gsl_multifit_nlinear_workspace *w);

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
    if (pars->Lw)
        gsl_matrix_free(pars->Lw);
    if (pars->mx)
        gsl_matrix_free(pars->mx);
    if (pars->mpopt)
        gsl_vector_free(pars->mpopt);
    if (pars->diag)
        gsl_vector_free(pars->diag);
    if (pars->lu)
        gsl_matrix_free(pars->lu);
    if (pars->JTJ)
        gsl_matrix_free(pars->JTJ);
    if (pars->workn)
        gsl_vector_free(pars->workn);
    if (pars->workp)
        gsl_vector_free(pars->workp);
    if (pars->worknp)
        gsl_matrix_free(pars->worknp);
    if (pars->mpopt1)
        gsl_vector_free(pars->mpopt1);
    if (pars->psi)
        gsl_vector_free(pars->psi);
    if (pars->psip)
        gsl_vector_free(pars->psip);
}

/* function call w/ cleanup */
SEXP C_nls(SEXP fn, SEXP y, SEXP jac, SEXP fvv, SEXP env, SEXP start, SEXP swts, SEXP lupars, SEXP control_int, SEXP control_dbl, SEXP has_start, SEXP loss_config)
{
    /* function arguments */
    pdata pars = {fn, y, jac, fvv, env, start, swts, lupars, control_int, control_dbl, has_start, loss_config,
                  NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

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
    params.startisnum = INTEGER_ELT(pars->control_int, 13);

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
    int wgt_i = INTEGER_ELT(VECTOR_ELT(pars->loss_config, 0), 0);
    pars->wts = gsl_vector_alloc(n);
    gsl_vector_set_all(pars->wts, 1.0);
    if (!Rf_isNull(pars->swts))
    {
        if (Rf_isMatrix(pars->swts))
        {
            /* move diagonal to pars->wts */
            pars->Lw = gsl_matrix_alloc(n, n);
            double *swts_ptr = REAL(pars->swts);
            for (R_len_t n1 = 0; n1 < n; n1++)
            {
                gsl_vector_set(pars->wts, n1, swts_ptr[n1 + n * n1] * swts_ptr[n1 + n * n1]);
                for (R_len_t n2 = 0; n2 < n; n2++)
                    gsl_matrix_set(pars->Lw, n1, n2, swts_ptr[n1 + n * n2] / swts_ptr[n1 + n * n1]);
            }
        }
        else
        {
            double *swts_ptr = REAL(pars->swts);
            for (R_len_t i = 0; i < n; i++)
                gsl_vector_set(pars->wts, i, swts_ptr[i] * swts_ptr[i]);
        }
    }

    /* initialize solver */
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;

    /* parameter constraints */
    if (Rf_isMatrix(pars->lupars))
    {
        double *luptr = REAL(pars->lupars);
        pars->lu = gsl_matrix_alloc(2, p);
        for (R_len_t k = 0; k < p; k++)
        {
            if (R_finite(luptr[2 * k]))
                gsl_matrix_set(pars->lu, 0, k, luptr[2 * k]);
            else
                gsl_matrix_set(pars->lu, 0, k, (double)GSL_NEGINF);
            if (R_finite(luptr[2 * k + 1]))
                gsl_matrix_set(pars->lu, 1, k, luptr[2 * k + 1]);
            else
                gsl_matrix_set(pars->lu, 1, k, (double)GSL_POSINF);
        }
    }

    /* allocate workspace with default parameters */
    pars->w = gsl_multifit_nlinear_alloc(T, &fdf_params, n, p);
    pars->mpopt = gsl_vector_alloc(p);
    pars->workp = gsl_vector_alloc(p);
    pars->workn = gsl_vector_alloc(n);
    pars->worknp = gsl_matrix_alloc(n, p);
    gsl_vector_set_zero(pars->mpopt);

    /* multistart algorithm */
    if (mstart)
    {
        /* allocate qrng workspace */
        if (p < 41)
            pars->q = gsl_qrng_alloc(gsl_qrng_sobol, p);
        else
            pars->q = gsl_qrng_alloc(gsl_qrng_halton, p);

        /* manually initialize workspace */
        (pars->w)->fdf = &fdf;
        (pars->w)->sqrt_wts = (pars->w)->sqrt_wts_work;
        if (!Rf_isNull(pars->swts))
        {
            for (R_len_t i = 0; i < n; i++)
            {
                double wi = gsl_vector_get(pars->wts, i);
                gsl_vector_set((pars->w)->sqrt_wts, i, sqrt(wi));
            }
        }
        else
            gsl_vector_set_all((pars->w)->sqrt_wts, 1.0);

        /* multistart parameters */
        mdata mpars;
        mpars.n = INTEGER_ELT(pars->control_int, 6);
        mpars.p = INTEGER_ELT(pars->control_int, 7);
        mpars.q = INTEGER_ELT(pars->control_int, 8);
        mpars.s = INTEGER_ELT(pars->control_int, 9);
        mpars.niter = INTEGER_ELT(pars->control_int, 10);
        mpars.max = INTEGER_ELT(pars->control_int, 11);
        mpars.minsp = INTEGER_ELT(pars->control_int, 12);
        mpars.wgt_i = wgt_i;
        mpars.all_start = TRUE;
        mpars.has_start = LOGICAL(pars->has_start);
        mpars.r = REAL_ELT(pars->control_dbl, 8);
        mpars.tol = REAL_ELT(pars->control_dbl, 9);
        mpars.dtol = 1.0e-6;
        mpars.ntix = (int *)S_alloc(mpars.n, sizeof(int));
        mpars.qmp = (double *)S_alloc(p, sizeof(double));
        mpars.mssr_order = (int *)R_alloc(mpars.n, sizeof(int));
        mpars.mstop = GSL_CONTINUE;
        mpars.mstarts = 0;
        mpars.nsp = 0;
        mpars.nwsp = 0;
        mpars.luchange = (int *)S_alloc(p, sizeof(int));
        mpars.rejectscl = 1.25;
        mpars.mssropt[0] = (double)GSL_POSINF;
        mpars.mssropt[1] = (double)GSL_POSINF;
        mpars.ssrconv[0] = 1.0;
        mpars.ssrconv[1] = 1.0;
        mpars.start = (double *)R_alloc(2 * p, sizeof(double));
        mpars.maxlims = (double *)R_alloc(2 * p, sizeof(double));

        // mpoptinfo mpoptall;
        // mpoptall.mpcount = 0;
        // mpoptall.mpmax = 100;
        // mpoptall.mpall = (double *)S_alloc(mpars.p * mpoptall.mpmax, sizeof(double));
        // mpoptall.mpradii = (double *)S_alloc(mpars.p * mpoptall.mpmax, sizeof(double));
        // mpars.mpopt = &mpoptall;

        SEXP mssr = PROTECT(Rf_allocVector(REALSXP, mpars.n));
        SEXP mpar = PROTECT(Rf_allocVector(REALSXP, p));
        Rf_setAttrib(mpar, R_NamesSymbol, parnames);
        params.start = mpar;
        nprotect += 2;

        /* workspace matrices/vectors */
        pars->mx = gsl_matrix_alloc(mpars.n, p);
        pars->diag = gsl_vector_alloc(p);
        pars->JTJ = gsl_matrix_alloc(p, p);
        pars->mpopt1 = gsl_vector_alloc(p);

        double *startptr = REAL(pars->start);
        for (R_len_t k = 0; k < p; k++)
        {
            mpars.start[2 * k] = startptr[2 * k];
            mpars.start[2 * k + 1] = startptr[2 * k + 1];
            mpars.maxlims[2 * k] = startptr[2 * k];
            mpars.maxlims[2 * k + 1] = startptr[2 * k + 1];
        }

        /* initially focus sampling around center points */
        for (R_len_t k = 0; k < p; k++)
        {
            if (!(mpars.has_start)[2 * k] || !(mpars.has_start)[2 * k + 1])
            {
                gsl_vector_set(pars->diag, k, 1.0);
                mpars.all_start = FALSE;
            }
            else
            {
                gsl_vector_set(pars->diag, k, 0.75);
                if ((mpars.start)[2 * k] + xtol > (mpars.start)[2 * k + 1])
                    mpars.rejectscl = -1.0; // starting range is not an interval
            }
        }

        /* multi-start global iterations */
        do
        {
            gsl_multistart_driver(pars, &mpars, &fdf, mssr, xtol, ftol, FALSE, verbose);

            /* check stopping criterion */
            mpars.mstarts += 1;

            if (mpars.mstarts > mpars.max)
                mpars.mstop = GSL_EMAXITER;

            if (mpars.nsp >= mpars.minsp && mpars.nwsp > (mpars.r + sqrt(mpars.r) * mpars.nsp))
                mpars.mstop = GSL_SUCCESS;

            if (!(mpars.mstarts % 10) && !(mpars.mssropt[0] < (double)GSL_POSINF))
            {
                /* reduce determinant tolerance */
                mpars.dtol = gsl_max(0.5 * mpars.dtol, GSL_DBL_EPSILON);
                if (!(mpars.mstarts % 100)) // reset start ranges
                {
                    for (R_len_t k = 0; k < p; k++)
                    {
                        mpars.start[2 * k] = startptr[2 * k];
                        mpars.start[2 * k + 1] = startptr[2 * k + 1];
                    }
                }
            }

        } while (mpars.mstop == GSL_CONTINUE);
        /* second pass multi-start global iterations */
        if (wgt_i)
        {
            /*
                outliers are identified as observations with Di > min(4 / n, 5 * mad(Di)),
                associated weights are set to zero in the second pass of the multi-start procedure.
            */
            if (mpars.mssropt[1] < mpars.mssropt[0])
                gsl_vector_memcpy(pars->mpopt, pars->mpopt1);
            // compute J and J^T * J using pars->mpopt
            if (pars->Lw)
                gsl_multifit_nlinear_winit_LD(pars->mpopt, pars->wts, pars->Lw, &fdf, pars->w);
            else
                gsl_multifit_nlinear_winit(pars->mpopt, pars->wts, &fdf, pars->w);

            det_eval_jtj((pars->w)->params, (pars->w)->sqrt_wts, pars->Lw, (pars->w)->fdf,
                         (pars->w)->x, (pars->w)->f, (pars->w)->J, pars->JTJ, pars->workn);

            // compute cook's d
            mpars.mstop = cooks_d((pars->w)->f, (pars->w)->J, pars->JTJ, pars->workn, pars->worknp);

            if (!mpars.mstop)
            {
                // update weights
                R_len_t noutlier = 0;
                double *workn = (double *)R_alloc(n, sizeof(double));
                double mad = gsl_stats_mad((pars->workn)->data, 1, n, workn);
                double thresh = gsl_min(4.0 / n, 5 * mad);
                (pars->w)->sqrt_wts = (pars->w)->sqrt_wts_work;
                for (R_len_t i = 0; i < n; i++)
                {
                    double di = gsl_vector_get(pars->workn, i);
                    if (di > thresh)
                    {
                        gsl_vector_set(pars->wts, i, 0.0);
                        gsl_vector_set((pars->w)->sqrt_wts, i, 0.0);
                        noutlier += 1;
                    }
                    else
                    {
                        double wi = gsl_vector_get(pars->wts, i);
                        gsl_vector_set((pars->w)->sqrt_wts, i, sqrt(wi));
                    }
                }
                if ((noutlier > 0) && (noutlier < (n - p)))
                {
                    /* reset qrng workspace */
                    gsl_qrng_init(pars->q);
                    // reset counters
                    mpars.mstop = GSL_CONTINUE;
                    mpars.mstarts = 0;
                    mpars.nsp = 0;
                    mpars.nwsp = 0;
                    mpars.dtol = 1.0e-6;
                    mpars.rejectscl = 1.25;
                    mpars.mssropt[0] = (double)GSL_POSINF;
                    mpars.mssropt[1] = (double)GSL_POSINF;
                    mpars.ssrconv[0] = 1.0;
                    mpars.ssrconv[1] = 1.0;
                    memset(mpars.ntix, 0, mpars.n * sizeof(int));
                    memset(mpars.luchange, 0, p * sizeof(int));
                    // second pass multi-start
                    do
                    {
                        gsl_multistart_driver(pars, &mpars, &fdf, mssr, xtol, ftol, TRUE, verbose);

                        /* check stopping criterion */
                        mpars.mstarts += 1;

                        if (mpars.mstarts > mpars.max)
                            mpars.mstop = GSL_EMAXITER;

                        if (mpars.nsp >= mpars.minsp && mpars.nwsp > (mpars.r + sqrt(mpars.r) * mpars.nsp))
                            mpars.mstop = GSL_SUCCESS;

                        if (!(mpars.mstarts % 10) && !(mpars.mssropt[0] < (double)GSL_POSINF))
                        {
                            /* reduce determinant tolerance */
                            mpars.dtol = gsl_max(0.5 * mpars.dtol, GSL_DBL_EPSILON);
                            if (!(mpars.mstarts % 100)) // reset start ranges
                            {
                                for (R_len_t k = 0; k < p; k++)
                                {
                                    mpars.start[2 * k] = startptr[2 * k];
                                    mpars.start[2 * k + 1] = startptr[2 * k + 1];
                                }
                            }
                        }
                    } while (mpars.mstop == GSL_CONTINUE);
                }
                // reset original weights
                // sqrt_wts is initialized by gsl_multifit_nlinear_winit
                gsl_vector_set_all(pars->wts, 1.0);
                if (!Rf_isNull(pars->swts))
                {
                    if (Rf_isMatrix(pars->swts))
                    {
                        double *swts_ptr = REAL(pars->swts);
                        for (R_len_t i = 0; i < n; i++)
                            gsl_vector_set(pars->wts, i, swts_ptr[i + n * i] * swts_ptr[i + n * i]);
                    }
                    else
                    {
                        double *swts_ptr = REAL(pars->swts);
                        for (R_len_t i = 0; i < n; i++)
                            gsl_vector_set(pars->wts, i, swts_ptr[i] * swts_ptr[i]);
                    }
                }
            }
        }
        if (verbose)
        {
            if (mpars.mstop == GSL_SUCCESS)
                Rprintf("multi-start algorithm finished successfully (NSP = %d, NWSP = %d, # iterations = %d)\n", mpars.nsp, mpars.nwsp, mpars.mstarts);
            if (mpars.mstop == GSL_EMAXITER)
                Rprintf("multi-start algorithm reached max. number of global iterations (NSP = %d, NWSP = %d, # iterations = %d)\n", mpars.nsp, mpars.nwsp, mpars.mstarts);
            Rprintf("*******************\n");
        }
        if (mpars.mssropt[1] < mpars.mssropt[0])
        {
            mpars.mssropt[0] = mpars.mssropt[1];
            mpars.ssrconv[0] = mpars.ssrconv[1];
            gsl_vector_memcpy(pars->mpopt, pars->mpopt1);
        }
        /* add small jitter to parameters */
        if (mpars.mssropt[0] < ftol || mpars.ssrconv[0] < ftol)
        {
            if (pars->lu)
                gsl_vector_set(pars->mpopt, 0, gsl_min(gsl_vector_get(pars->mpopt, 0) + 1.0e-4, gsl_matrix_get(pars->lu, 1, 0)));
            else
                gsl_vector_set(pars->mpopt, 0, gsl_vector_get(pars->mpopt, 0) + 1.0e-4);
        }
    }
    else /* single-start */
    {
        for (R_len_t k = 0; k < p; k++)
            gsl_vector_set(pars->mpopt, k, REAL_ELT(startvec, k));
    }

    /* (re-)initialize solver w/ optimal start parameters */
    if (pars->Lw)
        gsl_multifit_nlinear_winit_LD(pars->mpopt, pars->wts, pars->Lw, &fdf, pars->w);
    else if (!Rf_isNull(pars->swts) || wgt_i)
        gsl_multifit_nlinear_winit(pars->mpopt, pars->wts, &fdf, pars->w);
    else
        gsl_multifit_nlinear_init(pars->mpopt, &fdf, pars->w);

    /* initial cost */
    double chisq_init = (double)GSL_POSINF;
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
    int status = GSL_FAILURE;
    int irls_status = GSL_FAILURE;
    R_len_t irls_iter = 0;
    double irls_delta = 0.0;
    double irls_sigma = 1.0;

    if (!wgt_i) /* default loss function */
    {
        status = gsl_multifit_nlinear_driver2(niter, xtol, gtol, ftol, verbose ? callback : NULL,
                                              verbose ? &params : NULL, &info, &chisq0, &chisq1,
                                              pars->lu, pars->Lw, pars->w);
    }
    else /* non-default loss function */
    {
        /* irls psi and psi' */
        pars->psi = gsl_vector_alloc(n);
        pars->psip = gsl_vector_alloc(n);

        status = gsl_multifit_nlinear_rho_driver(pars, &fdf, wgt_i, niter, xtol, gtol, ftol, &params,
                                                 &info, &chisq0, &chisq1, &irls_sigma, &irls_iter, &irls_status, verbose);

        /* copy irls weights */
        gsl_vector_memcpy(pars->wts, pars->workn);

        /* irls x tolerance */
        for (R_len_t k = 0; k < p; k++)
        {
            double x0 = gsl_vector_get(pars->workp, k);
            double x1 = gsl_vector_get((pars->w)->x, k);
            irls_delta = gsl_max(irls_delta, fabs(x0 - x1));
        }
    }

    R_len_t iter = (R_len_t)gsl_multifit_nlinear_niter(pars->w);

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
        if (wgt_i)
        {
            Rprintf("IRLS number of iterations: %d\n", irls_iter);
            Rprintf("IRLS achieved tolerance: %g\n", irls_delta);
            Rprintf("IRLS convergence status: %s\n", gsl_strerror(irls_status));
        }
        Rprintf("number of iterations: %d\n", iter);
        // Rprintf("reason for stopping: %s\n", gsl_strerror(info));
        Rprintf("initial ssr: %g\n", chisq_init);
        Rprintf("final ssr: %g\n", chisq1);
        Rprintf("ssr/dof: %g\n", chisq1 / (n - p));
        Rprintf("ssr achieved tolerance: %g\n", chisq0 - chisq1);
        Rprintf("function evaluations: %d\n", (int)fdf.nevalf);
        Rprintf("jacobian evaluations: %d\n", (int)fdf.nevaldf);
        Rprintf("fvv evaluations: %d\n", (int)fdf.nevalfvv);
        Rprintf("status: %s\n*******************\n", gsl_strerror(status));
    }

    /* initialize result */
    SEXP ans = NULL;
    if (verbose)
    {
        const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
                                "algorithm", "neval", "irls", "partrace", "ssrtrace", ""};
        ans = PROTECT(Rf_mkNamed(VECSXP, ansnms));
    }
    else
    {
        const char *ansnms[] = {"par", "covar", "resid", "grad", "niter", "status", "conv", "ssr", "ssrtol",
                                "algorithm", "neval", "irls", ""};
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
    SET_INTEGER_ELT(ansneval, 0, (R_len_t)fdf.nevalf);
    SET_INTEGER_ELT(ansneval, 1, (R_len_t)fdf.nevaldf);
    SET_INTEGER_ELT(ansneval, 2, (R_len_t)fdf.nevalfvv);
    SET_VECTOR_ELT(ans, 10, ansneval);
    UNPROTECT(1);

    /* irls weights */
    if (wgt_i)
    {
        const char *irlsnms[] = {"irls_weights", "irls_psi", "irls_dpsi", "irls_sigma", "irls_status", "irls_niter", "irls_tol", "irls_conv", ""};
        SEXP ansirls = PROTECT(Rf_mkNamed(VECSXP, irlsnms));
        SEXP irlswts = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP irlspsi = PROTECT(Rf_allocVector(REALSXP, n));
        SEXP irlspsip = PROTECT(Rf_allocVector(REALSXP, n));
        if (status == GSL_SUCCESS || status == GSL_EMAXITER)
        {
            for (R_len_t i = 0; i < n; i++)
            {
                SET_REAL_ELT(irlswts, i, gsl_vector_get(pars->wts, i));
                SET_REAL_ELT(irlspsi, i, gsl_vector_get(pars->psi, i));
                SET_REAL_ELT(irlspsip, i, gsl_vector_get(pars->psip, i));
            }
        }
        else
        {
            for (R_len_t i = 0; i < n; i++)
            {
                SET_REAL_ELT(irlswts, i, NA_REAL);
                SET_REAL_ELT(irlspsi, i, NA_REAL);
                SET_REAL_ELT(irlspsip, i, NA_REAL);
            }
        }
        SET_VECTOR_ELT(ansirls, 0, irlswts);
        SET_VECTOR_ELT(ansirls, 1, irlspsi);
        SET_VECTOR_ELT(ansirls, 2, irlspsip);
        SET_VECTOR_ELT(ansirls, 3, Rf_ScalarReal(irls_sigma));
        SET_VECTOR_ELT(ansirls, 4, Rf_ScalarString(Rf_mkChar(gsl_strerror(irls_status))));
        SET_VECTOR_ELT(ansirls, 5, Rf_ScalarInteger(irls_iter));
        SET_VECTOR_ELT(ansirls, 6, Rf_ScalarReal(irls_delta));
        SET_VECTOR_ELT(ansirls, 7, Rf_ScalarInteger(irls_status));
        SET_VECTOR_ELT(ans, 11, ansirls);
        UNPROTECT(4);
    }

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
        SET_VECTOR_ELT(ans, 12, params.partrace);
        SET_VECTOR_ELT(ans, 13, params.ssrtrace);
    }

    /* free memory */
    if (status == GSL_SUCCESS || status == GSL_EMAXITER)
        gsl_matrix_free(cov);

    UNPROTECT(nprotect);
    return ans;
}

static int gsl_f(const gsl_vector *x, void *params, gsl_vector *f)
{
    /* construct parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    if (((fdata *)params)->startisnum)
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
            gsl_vector_set(f, i, (double)GSL_POSINF);
        else
            gsl_vector_set(f, i, fvalptr[i] - yptr[i]);
    }

    UNPROTECT(2);
    return GSL_SUCCESS;
}

static int gsl_df(const gsl_vector *x, void *params, gsl_matrix *J)
{
    /* construct parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    if (((fdata *)params)->startisnum)
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

static int gsl_fvv(const gsl_vector *x, const gsl_vector *v, void *params, gsl_vector *fvv)
{
    /* populate parameter vector */
    SEXP par = NULL;
    R_len_t p = ((fdata *)params)->p;
    SEXP start = ((fdata *)params)->start;
    SEXP parnames = PROTECT(Rf_getAttrib(start, R_NamesSymbol));
    if (((fdata *)params)->startisnum)
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

static void callback(const R_len_t iter, void *params, const gsl_multifit_nlinear_workspace *w)
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
