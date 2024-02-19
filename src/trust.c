#include "gsl_nls.h"

/* 
    these internal functions are copied from multifit_nlinear/trust.c 
    in order to modify the internal function trust_iterate() to include 
    parameter bounds constraints
*/

/* compute x_trial = x + dx */
static void
trust_trial_step_lu(const gsl_vector *x, const gsl_vector *dx,
                    gsl_vector *x_trial, const gsl_matrix *lu, double delta)
{
    size_t i, N = x->size;

    double dxi, xi, xi_trial, xi_lwr, xi_upr;
    for (i = 0; i < N; i++)
    {
        dxi = gsl_vector_get(dx, i);
        xi = gsl_vector_get(x, i);
        xi_trial = xi + dxi;

        xi_lwr = gsl_matrix_get(lu, 0, i);
        xi_upr = gsl_matrix_get(lu, 1, i);

        if (xi_trial < xi_lwr)
            xi_trial = xi + (dxi / gsl_max(fabs(dxi), delta) * fabs(xi - xi_lwr));
        else if (xi_trial > xi_upr)
            xi_trial = xi + (dxi / gsl_max(fabs(dxi), delta) * fabs(xi - xi_upr));

        gsl_vector_set(x_trial, i, xi_trial);
    }
}

static double
trust_calc_rho(const gsl_vector *f, const gsl_vector *f_trial,
               const gsl_vector *g, const gsl_matrix *J,
               const gsl_vector *dx, trust_state_t *state)
{
    int status;
    const gsl_multifit_nlinear_parameters *params = &(state->params);
    const gsl_multifit_nlinear_trs *trs = params->trs;
    const double normf = gsl_blas_dnrm2(f);
    const double normf_trial = gsl_blas_dnrm2(f_trial);
    gsl_multifit_nlinear_trust_state trust_state;
    double rho;
    double actual_reduction;
    double pred_reduction;
    double u;

    /* if ||f(x+dx)|| > ||f(x)|| reject step immediately */
    if (normf_trial >= normf)
        return -1.0;

    trust_state.x = NULL;
    trust_state.f = f;
    trust_state.g = g;
    trust_state.J = J;
    trust_state.diag = state->diag;
    trust_state.sqrt_wts = NULL;
    trust_state.mu = &(state->mu);
    trust_state.params = params;
    trust_state.solver_state = state->solver_state;
    trust_state.fdf = NULL;
    trust_state.avratio = &(state->avratio);

    /* compute numerator of rho (actual reduction) */
    u = normf_trial / normf;
    actual_reduction = 1.0 - u * u;

    /*
     * compute denominator of rho (predicted reduction); this is calculated
     * inside each trust region subproblem, since it depends on the local
     * model used, which can vary according to each TRS
     */
    status = (trs->preduction)(&trust_state, dx, &pred_reduction, state->trs_state);
    if (status)
        return -1.0;

    if (pred_reduction > 0.0)
        rho = actual_reduction / pred_reduction;
    else
        rho = -1.0;

    return rho;
}

/*
trust_eval_step()
  Evaluate proposed step to determine if it should be
accepted or rejected
*/

static int
trust_eval_step(const gsl_vector *f, const gsl_vector *f_trial,
                const gsl_vector *g, const gsl_matrix *J,
                const gsl_vector *dx, double *rho, trust_state_t *state)
{
    int status = GSL_SUCCESS;
    const gsl_multifit_nlinear_parameters *params = &(state->params);

    if (params->trs == gsl_multifit_nlinear_trs_lmaccel)
    {
        /* reject step if acceleration is too large compared to velocity */
        if (state->avratio > params->avmax)
            status = GSL_FAILURE;
    }

    /* compute rho */
    *rho = trust_calc_rho(f, f_trial, g, J, dx, state);
    if (*rho <= 0.0)
        status = GSL_FAILURE;

    return status;
}

static int
nielsen_accept(const double rho, double *mu, long *nu)
{
    double b;

    /* reset nu */
    *nu = 2;

    b = 2.0 * rho - 1.0;
    b = 1.0 - b * b * b;
    *mu *= gsl_max(0.333333333333333, b);

    return GSL_SUCCESS;
}

static int
nielsen_reject(double *mu, long *nu)
{
    *mu *= (double)*nu;

    /* nu := 2*nu */
    *nu <<= 1;

    return GSL_SUCCESS;
}

/*
trust_iterate_lu()
  This function performs 1 iteration of the trust region algorithm
  including lower/upper parameter constraints.
  It calls a user-specified method for computing the next step
  (LM or dogleg), then tests if the computed step is acceptable.

Args: vstate - trust workspace
      swts   - data weights (NULL if unweighted)
      fdf    - function and Jacobian pointers
      x      - on input, current parameter vector
               on output, new parameter vector x + dx
      f      - on input, f(x)
               on output, f(x + dx)
      J      - on input, J(x)
               on output, J(x + dx)
      g      - on input, g(x) = J(x)' f(x)
               on output, g(x + dx) = J(x + dx)' f(x + dx)
      dx     - (output only) parameter step vector
      lu     - (2 x p)-matrix with top row the lower
                parameter bounds and bottom row the
                upper parameter bounds

Return:
1) GSL_SUCCESS if we found a step which reduces the cost
function

2) GSL_ENOPROG if 15 successive attempts were to made to
find a good step without success

3) If a scaling matrix D is used, inputs and outputs are
set to the unscaled quantities (ie: J and g)
*/
int trust_iterate_lu(void *vstate, const gsl_vector *swts,
               gsl_multifit_nlinear_fdf *fdf, gsl_vector *x,
               gsl_vector *f, gsl_matrix *J, gsl_vector *g,
               gsl_vector *dx, const gsl_matrix *lu)
{
    int status;
    trust_state_t *state = (trust_state_t *)vstate;
    const gsl_multifit_nlinear_parameters *params = &(state->params);
    const gsl_multifit_nlinear_trs *trs = params->trs;
    gsl_multifit_nlinear_trust_state trust_state;
    gsl_vector *x_trial = state->x_trial; /* trial x + dx */
    gsl_vector *f_trial = state->f_trial; /* trial f(x + dx) */
    gsl_vector *diag = state->diag;       /* diag(D) */
    double rho;                           /* ratio actual_reduction/predicted_reduction */
    int foundstep = 0;                    /* found step dx */
    int bad_steps = 0;                    /* consecutive rejected steps */

    /* store all state parameters needed by low level methods */
    trust_state.x = x;
    trust_state.f = f;
    trust_state.g = g;
    trust_state.J = J;
    trust_state.diag = state->diag;
    trust_state.sqrt_wts = swts;
    trust_state.mu = &(state->mu);
    trust_state.params = params;
    trust_state.solver_state = state->solver_state;
    trust_state.fdf = fdf;
    trust_state.avratio = &(state->avratio);

    /* initialize trust region subproblem with this Jacobian */
    status = (trs->preloop)(&trust_state, state->trs_state);
    if (status)
        return status;

    /* loop until we find an acceptable step dx */
    while (!foundstep)
    {
        /* calculate new step */
        status = (trs->step)(&trust_state, state->delta, dx, state->trs_state);

        /* occasionally the iterative methods (ie: CG Steihaug) can fail to find a step,
         * so in this case skip rho calculation and count it as a rejected step */

        if (status == GSL_SUCCESS)
        {
            /* compute x_trial = x + dx */
            trust_trial_step_lu(x, dx, x_trial, lu, state->delta);

            /* compute f_trial = f(x + dx) */
            status = gsl_multifit_nlinear_eval_f(fdf, x_trial, swts, f_trial);
            if (status)
                return status;

            /* check if step should be accepted or rejected */
            status = trust_eval_step(f, f_trial, g, J, dx, &rho, state);
            if (status == GSL_SUCCESS)
                foundstep = 1;
        }
        else
        {
            /* an iterative TRS method failed to find a step vector */
            rho = -1.0;
        }

        /*
         * update trust region radius: if rho is large,
         * then the quadratic model is a good approximation
         * to the objective function, enlarge trust region.
         * If rho is small (or negative), the model function
         * is a poor approximation so decrease trust region. This
         * can happen even if the step is accepted.
         */
        if (rho > 0.75)
            state->delta *= params->factor_up;
        else if (rho < 0.25)
            state->delta /= params->factor_down;

        if (foundstep)
        {
            /* step was accepted */

            /* compute J <- J(x + dx) */
            status = gsl_multifit_nlinear_eval_df(x_trial, f_trial, swts,
                                                  params->h_df, params->fdtype,
                                                  fdf, J, state->workn);
            if (status)
                return status;

            /* update x <- x + dx */
            gsl_vector_memcpy(x, x_trial);

            /* update f <- f(x + dx) */
            gsl_vector_memcpy(f, f_trial);

            /* compute new g = J^T f */
            gsl_blas_dgemv(CblasTrans, 1.0, J, f, 0.0, g);

            /* update scaling matrix D */
            (params->scale->update)(J, diag);

            /* step accepted, decrease LM parameter */
            status = nielsen_accept(rho, &(state->mu), &(state->nu));
            if (status)
                return status;

            bad_steps = 0;
        }
        else
        {
            /* step rejected, increase LM parameter */
            status = nielsen_reject(&(state->mu), &(state->nu));
            if (status)
                return status;

            if (++bad_steps > 15)
            {
                /* if more than 15 consecutive rejected steps, report no progress */
                return GSL_ENOPROG;
            }
        }
    }

    return GSL_SUCCESS;
} /* trust_iterate_lu() */
