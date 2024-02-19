#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
void R_init_gslnls(DllInfo *dll);
extern SEXP C_nls(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_large(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_test_f(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_test_j(SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_test_start_sol(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_nls", (DL_FUNC)&C_nls, 11},
    {"C_nls_large", (DL_FUNC)&C_nls_large, 9},
    {"C_nls_test_f", (DL_FUNC)&C_nls_test_f, 4},
    {"C_nls_test_j", (DL_FUNC)&C_nls_test_j, 4},
    {"C_nls_test_start_sol", (DL_FUNC)&C_nls_test_start_sol, 3},
    {NULL, NULL, 0}};

void R_init_gslnls(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
