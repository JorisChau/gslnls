#define R_NO_REMAP

#include <math.h>
#include "gsl_nls.h"

void gsl_multistart_driver(pdata *pars,
                                  mdata *mpars,
                                  gsl_multifit_nlinear_fdf *fdff,
                                  SEXP mssr,
                                  const double xtol,
                                  const double ftol,
                                  Rboolean verbose)
{
    // initialize variables
    int minfo;
    R_len_t p = (R_len_t)((pars->workp)->size);
    double kd, l0, l1, diagmin;
    // double rcond;
    double det_jtj, mchisq0 = (double)GSL_POSINF, mchisq1 = (double)GSL_POSINF;
    trust_state_t *trust_state = (trust_state_t *)(pars->w->state);

    /* sample initial points */
    for (R_len_t nn = 0; nn < mpars->n; nn++)
    {
        SET_REAL_ELT(mssr, nn, NA_REAL);

        if ((mpars->ntix)[nn] == 0)
        {
            gsl_qrng_get(pars->q, mpars->qmp);
            for (R_len_t k = 0; k < p; k++)
            {
                l0 = (mpars->start)[2 * k];
                l1 = (mpars->start)[2 * k + 1];
                if (l1 > l0)
                {
                    kd = gsl_vector_get(pars->diag, k);
                    (mpars->qmp)[k] = l0 + (l1 - l0) * (mpars->qmp)[k];
                    if (l0 > 0.0)
                        gsl_matrix_set(pars->mx, nn, k, (pow((mpars->qmp)[k] - l0 + 1.0, kd) - 1.0) / kd + l0);
                    else if (l1 < 0.0)
                        gsl_matrix_set(pars->mx, nn, k, -(pow(-(mpars->qmp)[k] + l1 + 1.0, kd) - 1.0) / kd + l1);
                    else if ((mpars->qmp)[k] > 0.0)
                        gsl_matrix_set(pars->mx, nn, k, (pow((mpars->qmp)[k] + 1.0, kd) - 1.0) / kd);
                    else
                        gsl_matrix_set(pars->mx, nn, k, -(pow(-(mpars->qmp)[k] + 1.0, kd) - 1.0) / kd);
                }
                else
                {
                    gsl_matrix_set(pars->mx, nn, k, l0);
                }
            }
        }
        /* calculate det(J^T * J) */
        gsl_vector_view nnx = gsl_matrix_row(pars->mx, nn);
        gsl_vector_memcpy((pars->w)->x, &nnx.vector);
        det_jtj = det_eval_jtj((pars->w)->params, (pars->w)->sqrt_wts, (pars->w)->fdf, (pars->w)->x, (pars->w)->f, (pars->w)->J, pars->JTJ, pars->workn);

        if (det_jtj > mpars->dtol)
        {
            gsl_matrix_get_row(pars->workp, pars->mx, nn);

            /* concentrate point */
            if (!Rf_isNull(pars->weights))
                gsl_multifit_nlinear_winit(pars->workp, pars->wts, fdff, pars->w);
            else
                gsl_multifit_nlinear_init(pars->workp, fdff, pars->w);
            gsl_multifit_nlinear_driver2(mpars->p, xtol, 1e-3, ftol, NULL, NULL, &minfo, &mchisq0, &mchisq1, pars->lu, pars->w);

            det_jtj = det_cholesky_jtj((pars->w)->J, pars->JTJ);

            if (mchisq1 < (double)GSL_POSINF)
            {
                if (det_jtj > mpars->dtol)
                {
                    gsl_matrix_set_row(pars->mx, nn, (pars->w)->x);
                    SET_REAL_ELT(mssr, nn, mchisq1);
                    if (mchisq1 < 0.99 * gsl_min((mpars->mssropt)[0], (mpars->mssropt)[1]))
                    {
                        (mpars->mssropt)[0] = mchisq1;
                        (mpars->ssrconv)[0] = mchisq0 - mchisq1;
                        gsl_vector_memcpy(pars->mpopt, (pars->w)->x);
                    }
                }
                else if (mchisq1 < 0.99 * gsl_min((mpars->mssropt)[0], (mpars->mssropt)[1]))
                {
                    (mpars->mssropt)[1] = mchisq1;
                    (mpars->ssrconv)[1] = mchisq0 - mchisq1;
                    gsl_vector_memcpy(pars->mpopt1, (pars->w)->x);
                }
            }
        }
        else if (!((mpars->mssropt)[0] < (double)GSL_POSINF) && det_jtj > GSL_DBL_EPSILON)
        {
            // back-up in case no stationary points found
            gsl_blas_ddot((pars->w)->f, (pars->w)->f, &mchisq1);
            if (mchisq1 < 0.99 * (mpars->mssropt)[1])
            {
                (mpars->mssropt)[1] = mchisq1;
                (mpars->ssrconv)[1] = mchisq0 - mchisq1;
                gsl_vector_memcpy(pars->mpopt1, (pars->w)->x);
            }
        }
    }

    /* reduce sample points */
    R_orderVector1(mpars->mssr_order, mpars->n, mssr, TRUE, FALSE);
    for (R_len_t nn = 0; nn < mpars->n; nn++)
    {
        if (nn < (mpars->q) && !R_IsNA(REAL_ELT(mssr, (mpars->mssr_order)[nn])))
            (mpars->ntix)[(mpars->mssr_order)[nn]] += 1;
        else
            (mpars->ntix)[(mpars->mssr_order)[nn]] = 0;
    }

    /* dynamic lower/upper limits */
    if (!(mpars->all_start))
    {
        double pk, pmin = 0.0, pmax = 1.0;
        double mssr_diff = REAL_ELT(mssr, (mpars->mssr_order)[0]);

        if (!R_IsNA(mssr_diff))
        {
            for (R_len_t nn = mpars->n - 1; nn > 0; nn--)
            {
                if (!R_IsNA(REAL_ELT(mssr, (mpars->mssr_order)[nn])))
                {
                    mssr_diff -= REAL_ELT(mssr, (mpars->mssr_order)[nn]);
                    break;
                }
            }
        }
        if (R_IsNA(mssr_diff) || fabs(mssr_diff) < 1e-5)
        {
            for (R_len_t k = 0; k < p; k++)
                (mpars->luchange)[k] += 1;
        }

        for (R_len_t k = 0; k < p; k++)
        {
            /* evaluate min/max parameter values for reduced sample */
            if ((mpars->mssropt)[0] < (double)GSL_POSINF)
            {
                pmin = gsl_vector_get(((mpars->mssropt)[1] < (mpars->mssropt)[0]) ? pars->mpopt1 : pars->mpopt, k);
                pmax = gsl_vector_get(((mpars->mssropt)[1] < (mpars->mssropt)[0]) ? pars->mpopt1 : pars->mpopt, k);
            }
            for (R_len_t nn = 0; nn < mpars->q; nn++)
            {
                if ((mpars->ntix)[(mpars->mssr_order)[nn]] > 0 && REAL_ELT(mssr, (mpars->mssr_order)[nn]) < 1.25 * (mpars->mssropt)[0])
                {
                    pk = gsl_matrix_get(pars->mx, (mpars->mssr_order)[nn], k);
                    pmin = (pk < pmin) ? pk : pmin;
                    pmax = (pk > pmax) ? pk : pmax;
                }
            }

            /* rescale current limits */
            l0 = (mpars->start)[2 * k];
            l1 = (mpars->start)[2 * k + 1];
            int luchange_add = 0;

            // lower limit
            if (!(mpars->has_start)[2 * k])
            {
                if (pmin < 0.9 * l0 || (mpars->luchange)[k] > 4) // enlarge
                {
                    (mpars->start)[2 * k] = l0 < 0 ? gsl_max(l0 / pow(-1e-5 * (l0 - 1.0), 0.1) - 1.0, -1.0E5) : -0.1;
                    if (pars->lu)
                        (mpars->start)[2 * k] = gsl_max((mpars->start)[2 * k], gsl_matrix_get(pars->lu, 0, k));
                    (mpars->maxlims)[2 * k] = gsl_min((mpars->start)[2 * k], (mpars->maxlims)[2 * k]);
                    luchange_add = -1;
                }
                else if (pmin > 0.2 * l0) // shrink
                {
                    (mpars->start)[2 * k] = gsl_min(l0 / pow(-0.05 * (l0 - 1.0), 0.05), -0.01);
                    if (pars->lu)
                        (mpars->start)[2 * k] = gsl_max((mpars->start)[2 * k], gsl_matrix_get(pars->lu, 0, k));
                    luchange_add = ((mpars->mssropt)[0] < (double)GSL_POSINF) ? -1 : 1;
                }
                else
                    luchange_add = 1;
            }
            // upper limit
            if (!(mpars->has_start)[2 * k + 1])
            {
                if (pmax > 0.9 * l1 || (mpars->luchange)[k] > 4) // enlarge
                {
                    (mpars->start)[2 * k + 1] = gsl_min(l1 / pow(1e-5 * (l1 + 1.0), 0.1) + 1.0, 1.0E5);
                    if (pars->lu)
                        (mpars->start)[2 * k + 1] = gsl_min((mpars->start)[2 * k + 1], gsl_matrix_get(pars->lu, 1, k));
                    (mpars->maxlims)[2 * k + 1] = gsl_max((mpars->start)[2 * k + 1], (mpars->maxlims)[2 * k + 1]);
                    luchange_add = -1;
                }
                else if (pmax < 0.2 * l1) // shrink
                {
                    (mpars->start)[2 * k + 1] = gsl_max(l1 / pow(0.05 * (l1 + 1.0), 0.05), 0.1);
                    if (pars->lu)
                        (mpars->start)[2 * k + 1] = gsl_min((mpars->start)[2 * k + 1], gsl_matrix_get(pars->lu, 1, k));
                    luchange_add = ((mpars->mssropt)[0] < (double)GSL_POSINF) ? -1 : 1;
                }
                else
                    luchange_add = 1;
            }
            if (luchange_add)
            {
                (mpars->luchange)[k] = (luchange_add > 0) ? (mpars->luchange)[k] + 1 : 0;
            }
        }
    }

    /* local optimization stage */
    for (R_len_t nn = 0; nn < (mpars->n); nn++)
    {
        if ((mpars->ntix)[nn] >= (mpars->s))
        {
            (mpars->ntix)[nn] = 0;
            mpars->nwsp += 1;

            if ((mpars->nsp) == 0 || REAL_ELT(mssr, nn) < (1 + (mpars->tol)) * (mpars->mssropt)[0])
            {
                gsl_matrix_get_row(pars->workp, pars->mx, nn);
                if (!Rf_isNull(pars->weights))
                    gsl_multifit_nlinear_winit(pars->workp, pars->wts, fdff, pars->w);
                else
                    gsl_multifit_nlinear_init(pars->workp, fdff, pars->w);
                mchisq1 = REAL_ELT(mssr, nn);
                gsl_multifit_nlinear_driver2(mpars->niter, xtol, 1e-3, ftol, NULL, NULL, &minfo, &mchisq0, &mchisq1, pars->lu, pars->w);
                det_jtj = det_cholesky_jtj((pars->w)->J, pars->JTJ);

                // gsl_multifit_nlinear_rcond(&rcond, pars->w);

                /* save result local optimizer */
                // if ((mpars->mpopt)->mpcount > (mpars->mpopt)->mpmax)
                // {
                //     (mpars->mpopt)->mpall = (double *)S_realloc((char *)(mpars->mpopt)->mpall, 2 * (mpars->mpopt)->mpmax * p, (mpars->mpopt)->mpmax * p, sizeof(double));
                //     (mpars->mpopt)->mpradii = (double *)S_realloc((char *)(mpars->mpopt)->mpradii, 2 * (mpars->mpopt)->mpmax * p, (mpars->mpopt)->mpmax * p, sizeof(double));
                //     (mpars->mpopt)->mpmax *= 2;
                // }
                // for (R_len_t k = 0; k < p; k++)
                // {
                //     ((mpars->mpopt)->mpall)[(mpars->mpopt)->mpcount * p + k] = gsl_vector_get((pars->w)->x, k);
                //     ((mpars->mpopt)->mpradii)[(mpars->mpopt)->mpcount * p + k] = gsl_vector_get(pars->workp, k);
                // }

                // if (1 && verbose)
                // {
                //     Rprintf("%d opt: (", (mpars->mpopt)->mpcount);
                //     for (R_len_t k = 0; k < p; k++)
                //         Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", ((mpars->mpopt)->mpall)[(mpars->mpopt)->mpcount * p + k]);

                //     Rprintf("%d radius: (", (mpars->mpopt)->mpcount);
                //     for (R_len_t k = 0; k < p; k++)
                //         Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", ((mpars->mpopt)->mpradii)[(mpars->mpopt)->mpcount * p + k]);
                // }

                // (mpars->mpopt)->mpcount += 1;

                // if (0 && verbose)
                // {
                //     Rprintf("mssr*0:%g, mssr*1: %g, mchisq1: %g, det(JTJ): %g, rejectscale: %g\n", (mpars->mssropt)[0], (mpars->mssropt)[1], mchisq1, det_jtj, mpars->rejectscl);
                //     Rprintf("{opt, lwr, upr} = {");
                //     for (R_len_t k = 0; k < p; k++)
                //         Rprintf("(%g, %g, %g)%s", gsl_vector_get((pars->w)->x, k), (mpars->start)[2 * k], (mpars->start)[2 * k + 1], k < (p - 1) ? "" : "}\n");
                // }

                if (mchisq1 < (double)GSL_POSINF && ((mpars->nsp) == 0 || mchisq1 < 0.99 * (mpars->mssropt)[0]) && (det_jtj > (mpars->dtol) || mchisq1 < (2 * ftol)))
                {
                    int reject = 0;
                    if (mpars->rejectscl > 0)
                    {
                        for (R_len_t k = 0; k < p; k++)
                        {
                            double xk = gsl_vector_get((pars->w)->x, k);
                            if (mpars->all_start)
                                reject += (xk > gsl_max((mpars->maxlims)[2 * k + 1], 1.0) || xk < gsl_min((mpars->maxlims)[2 * k], -1.0));
                            else
                                // reject += (xk > pow((mpars->start)[2 * k + 1], mpars->rejectscl) || xk < -pow(-(mpars->start)[2 * k], mpars->rejectscl));
                                reject += (xk > gsl_max(pow((mpars->maxlims)[2 * k + 1], mpars->rejectscl), 1.0) || xk < gsl_min(-pow(-(mpars->maxlims)[2 * k], mpars->rejectscl), -1.0));
                            if (reject > 0)
                                break;
                        }
                        if (!(mpars->all_start))
                            mpars->rejectscl += 0.05;
                    }

                    if (!reject)
                    {
                        (mpars->mssropt)[0] = mchisq1;
                        (mpars->ssrconv)[0] = mchisq0 - mchisq1;
                        gsl_vector_memcpy(pars->mpopt, (pars->w)->x);
                        mpars->nsp += 1;
                        mpars->nwsp = 0;
                        if (mpars->rejectscl > 0)
                            mpars->rejectscl = 1.25;
                        if (mpars->all_start)
                        {
                            diagmin = gsl_vector_min(trust_state->diag);
                            for (R_len_t k = 0; k < p; k++)
                                gsl_vector_set(pars->diag, k, pow(diagmin / gsl_vector_get(trust_state->diag, k), 0.25));
                            // Rprintf("diag = (");
                            // for (R_len_t k = 0; k < p; k++)
                            //     Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", gsl_vector_get(pars->diag, k));
                        }
                        if (verbose)
                        {
                            Rprintf("mstart ssr* = %g, det(JTJ) = %g, NSP = %d, NWSP = %d, par = (", (mpars->mssropt)[0], det_jtj, mpars->nsp, mpars->nwsp);
                            for (R_len_t k = 0; k < p; k++)
                                Rprintf((k < (p - 1)) ? "%g, " : "%g)\n", gsl_vector_get((pars->w)->x, k));
                        }
                    }
                }
                else if (mchisq1 < 0.99 * gsl_min((mpars->mssropt)[0], (mpars->mssropt)[1]))
                {
                    // back-up in case no stationary points found
                    (mpars->mssropt)[1] = mchisq1;
                    (mpars->ssrconv)[1] = mchisq0 - mchisq1;
                    gsl_vector_memcpy(pars->mpopt1, (pars->w)->x);
                }
            }
        }
    }
}
