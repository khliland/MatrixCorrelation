#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP MatrixCorrelation_RV2cpp(SEXP, SEXP);
extern SEXP MatrixCorrelation_significantRcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP MatrixCorrelation_significantRepRcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"MatrixCorrelation_RV2cpp",             (DL_FUNC) &MatrixCorrelation_RV2cpp,             2},
    {"MatrixCorrelation_significantRcpp",    (DL_FUNC) &MatrixCorrelation_significantRcpp,    4},
    {"MatrixCorrelation_significantRepRcpp", (DL_FUNC) &MatrixCorrelation_significantRepRcpp, 8},
    {NULL, NULL, 0}
};

void R_init_MatrixCorrelation(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
