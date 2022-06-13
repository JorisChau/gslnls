#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
void R_init_gslnls(DllInfo *dll);
extern SEXP C_nls(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_large(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_nls_multistart(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"C_nls", (DL_FUNC)&C_nls, 9},
    {"C_nls_large", (DL_FUNC)&C_nls_large, 9},
    {"C_nls_multistart", (DL_FUNC)&C_nls_multistart, 9},
    {NULL, NULL, 0}};

void R_init_gslnls(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
