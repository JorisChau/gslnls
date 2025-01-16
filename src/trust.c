#include "gsl_nls.h"

/*
    these internal functions are adapted from multifit_nlinear/trust.c
    modifying trust_iterate() to allow for parameter bounds constraints
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

static void
trust_trial_step_default(const gsl_vector *x, const gsl_vector *dx, gsl_vector *x_trial)
{
    size_t i, N = x->size;

    for (i = 0; i < N; i++)
    {
        double dxi = gsl_vector_get(dx, i);
        double xi = gsl_vector_get(x, i);
        gsl_vector_set(x_trial, i, xi + dxi);
    }
}

/* compute || diag(D) a || */
static double
trust_scaled_norm(const gsl_vector *D, const gsl_vector *a)
{
    const size_t n = a->size;
    double e2 = 0.0;
    size_t i;

    for (i = 0; i < n; ++i)
    {
        double Di = gsl_vector_get(D, i);
        double ai = gsl_vector_get(a, i);
        double u = Di * ai;

        e2 += u * u;
    }

    return sqrt(e2);
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
nielsen_init(const gsl_matrix *J, const gsl_vector *diag,
             double *mu, long *nu)
{
    const double mu0 = 1.0e-3;
    const size_t p = J->size2;
    size_t j;
    double max = -1.0;

    *nu = 2;

    /* set mu = mu0 * max(diag(J~^T J~)), with J~ = J D^{-1} */

    for (j = 0; j < p; ++j)
    {
        gsl_vector_const_view v = gsl_matrix_const_column(J, j);
        double dj = gsl_vector_get(diag, j);
        double norm = gsl_blas_dnrm2(&v.vector) / dj;
        max = GSL_MAX(max, norm);
    }

    *mu = mu0 * max * max;

    return GSL_SUCCESS;
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

/* compute z = alpha*x + beta*y */
static void
scaled_addition(const double alpha, const gsl_vector *x,
                const double beta, const gsl_vector *y, gsl_vector *z)
{
    const size_t N = z->size;
    size_t i;

    for (i = 0; i < N; i++)
    {
        double xi = gsl_vector_get(x, i);
        double yi = gsl_vector_get(y, i);
        gsl_vector_set(z, i, alpha * xi + beta * yi);
    }
}

/*
lm_step_LD()
  Calculate a new step vector by solving the linear
  least squares system and apply weighting transform using 
  LDDL' decomposition of W if given:
*/
static int
lm_step_LD(const void *vtrust_state, const double delta,
        gsl_vector *dx, void *vstate, const gsl_vector *Dw, const gsl_matrix *Lw)
{
    int status;
    const gsl_multifit_nlinear_trust_state *trust_state =
        (const gsl_multifit_nlinear_trust_state *)vtrust_state;
    lm_state_t *state = (lm_state_t *)vstate;
    const gsl_multifit_nlinear_parameters *params = trust_state->params;
    const double mu = *(trust_state->mu);

    (void)delta;

    /* prepare the linear solver with current LM parameter mu */
    status = (params->solver->presolve)(mu, trust_state, trust_state->solver_state);
    if (status)
        return status;

    /*
     * solve: [     J      ] v = - [ f ]
     *        [ sqrt(mu)*D ]       [ 0 ]
     */
    status = (params->solver->solve)(trust_state->f,
                                     state->vel,
                                     trust_state,
                                     trust_state->solver_state);
    if (status)
        return status;

    if (state->accel)
    {
        double anorm, vnorm;

        /* compute geodesic acceleration */
        status = gsl_multifit_nlinear_eval_fvv_LD(params->h_fvv,
                                               trust_state->x,
                                               state->vel,
                                               trust_state->f,
                                               trust_state->J,
                                               Dw,
                                               Lw, 
                                               trust_state->fdf,
                                               state->fvv,
                                               state->workp);
        if (status)
            return status;

        /*
         * solve: [     J      ] a = - [ fvv ]
         *        [ sqrt(mu)*D ]       [  0  ]
         */
        status = (params->solver->solve)(state->fvv,
                                         state->acc,
                                         trust_state,
                                         trust_state->solver_state);
        if (status)
            return status;

        anorm = gsl_blas_dnrm2(state->acc);
        vnorm = gsl_blas_dnrm2(state->vel);

        /* store |a| / |v| */
        *(trust_state->avratio) = anorm / vnorm;
    }

    /* compute step dx = v + 1/2 a */
    scaled_addition(1.0, state->vel, 0.5, state->acc, dx);

    return GSL_SUCCESS;
}

/*
trust_init_LD()
  Initialize trust region solver with general weight matrix

Inputs: vstate - workspace
        Dw     - diagonal D in LDDL' decomposition of W
        Lw     - unit lower triangular matrix L in LDDL' decomposition of W
                 set to NULL for unweighted fit
        fdf    - user callback functions
        x      - initial parameter values
        f      - (output) f(x) vector
        J      - (output) J(x) matrix
        g      - (output) J(x)' f(x) vector

Return: success/error
*/

int trust_init_LD(void *vstate, const gsl_vector *Dw, const gsl_matrix *Lw,
           gsl_multifit_nlinear_fdf *fdf, const gsl_vector *x,
           gsl_vector *f, gsl_matrix *J, gsl_vector *g)
{
    int status;
    trust_state_t *state = (trust_state_t *)vstate;
    const gsl_multifit_nlinear_parameters *params = &(state->params);
    double Dx;

    /* evaluate function and Jacobian at x and apply weight transform */
    status = gsl_multifit_nlinear_eval_f_LD(fdf, x, Dw, Lw, f);
    if (status)
        return status;

    status = gsl_multifit_nlinear_eval_df_LD(x, f, Dw, Lw, params->h_df,
                                          params->fdtype, fdf, J, state->workn);
    if (status)
        return status;

    /* compute g = J^T f */
    gsl_blas_dgemv(CblasTrans, 1.0, J, f, 0.0, g);

    /* initialize diagonal scaling matrix D */
    (params->scale->init)(J, state->diag);

    /* compute initial trust region radius */
    Dx = trust_scaled_norm(state->diag, x);
    state->delta = 0.3 * GSL_MAX(1.0, Dx);

    /* initialize LM parameter */
    status = nielsen_init(J, state->diag, &(state->mu), &(state->nu));
    if (status)
        return status;

    /* initialize trust region method solver */
    {
        gsl_multifit_nlinear_trust_state trust_state;

        trust_state.x = x;
        trust_state.f = f;
        trust_state.g = g;
        trust_state.J = J;
        trust_state.diag = state->diag;
        trust_state.sqrt_wts = Dw;
        trust_state.mu = &(state->mu);
        trust_state.params = params;
        trust_state.solver_state = state->solver_state;
        trust_state.fdf = fdf;
        trust_state.avratio = &(state->avratio);

        status = (params->trs->init)(&trust_state, state->trs_state);

        if (status)
            return status;
    }

    /* set default parameters */

    state->avratio = 0.0;

    return GSL_SUCCESS;
}

/*
trust_iterate_lu_LD()
  This function performs 1 iteration of the trust region algorithm
  including parameter constraints and/or a general weight matrix.
  It calls a user-specified method for computing the next step
  (LM or dogleg), then tests if the computed step is acceptable.

Args: vstate - trust workspace
      Dw     - diagonal D in LDDL' decomposition of W
      Lw     - optional unit lower triangular matrix L in LDDL' decomposition of W
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
      lu     - optional (2 x p)-matrix with top row the lower
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
int trust_iterate_lu_LD(void *vstate, const gsl_vector *Dw,
                     const gsl_matrix *Lw,
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
    trust_state.sqrt_wts = Dw;
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
        if (Lw && trs == gsl_multifit_nlinear_trs_lmaccel)
            status = lm_step_LD(&trust_state, state->delta, dx, state->trs_state, Dw, Lw);
        else
            status = (trs->step)(&trust_state, state->delta, dx, state->trs_state);

        /* occasionally the iterative methods (ie: CG Steihaug) can fail to find a step,
         * so in this case skip rho calculation and count it as a rejected step */

        if (status == GSL_SUCCESS)
        {
            /* compute x_trial = x + dx */
            if(lu)
                trust_trial_step_lu(x, dx, x_trial, lu, state->delta);
            else
                trust_trial_step_default(x, dx, x_trial);

            /* compute f_trial = f(x + dx) */
            if(Lw)
                status = gsl_multifit_nlinear_eval_f_LD(fdf, x_trial, Dw, Lw, f_trial);
            else
                status = gsl_multifit_nlinear_eval_f(fdf, x_trial, Dw, f_trial);

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
            if(Lw)
                status = gsl_multifit_nlinear_eval_df_LD(x_trial, f_trial, Dw, Lw,
                                                  params->h_df, params->fdtype,
                                                  fdf, J, state->workn);
            else
                status = gsl_multifit_nlinear_eval_df(x_trial, f_trial, Dw,
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
