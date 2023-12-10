#define R_NO_REMAP

#include "gsl_nls.h"

SEXP C_nls_test_f(SEXP id, SEXP p, SEXP n, SEXP x)
{
    // inputs/outputs
    int nprob = INTEGER_ELT(id, 0) - 27; // problem number
    int m_f = INTEGER_ELT(n, 0);         // # equations
    int n_f = INTEGER_ELT(p, 0);         // # parameters
    double *x_f = (double *)R_alloc(n_f, sizeof(double));
    double *fval_f = (double *)R_alloc(m_f, sizeof(double));
    for (int k = 0; k < n_f; k++)
        x_f[k] = REAL_ELT(x, k);

    // resolve function values
    F77_CALL(p00_f)
    (&nprob, &m_f, &n_f, x_f, fval_f);

    // return function values
    // if undefined returns NA (can this actually happen?)
    SEXP fval = PROTECT(Rf_allocVector(REALSXP, m_f));
    for (int i = 0; i < m_f; i++)
        SET_REAL_ELT(fval, i, fval_f[i]);

    UNPROTECT(1);
    return fval;
}

SEXP C_nls_test_j(SEXP id, SEXP p, SEXP n, SEXP x)
{
    // inputs/outputs
    int nprob = INTEGER_ELT(id, 0) - 27; // problem number
    int m_f = INTEGER_ELT(n, 0);         // # equations
    int n_f = INTEGER_ELT(p, 0);         // # parameters
    double *x_f = (double *)R_alloc(n_f, sizeof(double));
    double *jac_f = (double *)R_alloc(m_f * n_f, sizeof(double));
    for (int k = 0; k < n_f; k++)
        x_f[k] = REAL_ELT(x, k);

    // resolve function values
    F77_CALL(p00_j)
    (&nprob, &m_f, &n_f, x_f, jac_f);

    // return function values
    // if undefined returns NA (can this actually happen?)
    SEXP jac = PROTECT(Rf_allocMatrix(REALSXP, m_f, n_f));
    double *jacptr = REAL(jac);
    for (int k = 0; k < n_f; k++)
        for (int i = 0; i < m_f; i++)
                jacptr[i + k * m_f] = jac_f[i + k * m_f];

    UNPROTECT(1);
    return jac;
}

SEXP C_nls_test_start_sol(SEXP id, SEXP p, SEXP n)
{
    // inputs/outputs
    int nprob = INTEGER_ELT(id, 0) - 27; // problem number
    int m_f = INTEGER_ELT(n, 0);         // # equations
    int n_f = INTEGER_ELT(p, 0);         // # parameters
    int known_f;
    double *start_f = (double *)R_alloc(n_f, sizeof(double));
    double *target_f = (double *)S_alloc(n_f, sizeof(double));

    // resolve start values
    F77_CALL(p00_start)
    (&nprob, &n_f, start_f);

    // resolve target values (if any)
    F77_CALL(p00_sol)
    (&nprob, &m_f, &n_f, &known_f, target_f);

    // return named list
    // if target is not known returns vector of NAs
    SEXP start = PROTECT(Rf_allocVector(REALSXP, n_f));
    SEXP target = PROTECT(Rf_allocVector(REALSXP, n_f));

    for (int k = 0; k < n_f; k++)
    {
        SET_REAL_ELT(start, k, start_f[k]);
        if (known_f)
            SET_REAL_ELT(target, k, target_f[k]);
        else
            SET_REAL_ELT(target, k, NA_REAL);
    }

    SEXP start_sol = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(start_sol, 0, start);
    SET_VECTOR_ELT(start_sol, 1, target);

    SEXP nms = PROTECT(Rf_allocVector(STRSXP, 2));
    SET_STRING_ELT(nms, 0, Rf_mkChar("start"));
    SET_STRING_ELT(nms, 1, Rf_mkChar("target"));
    Rf_setAttrib(start_sol, R_NamesSymbol, nms);

    UNPROTECT(4);
    return start_sol;
}
