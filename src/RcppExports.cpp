// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// significantRcpp
NumericVector significantRcpp(SEXP smi1, SEXP T1, SEXP U1, SEXP B1);
RcppExport SEXP MatrixCorrelation_significantRcpp(SEXP smi1SEXP, SEXP T1SEXP, SEXP U1SEXP, SEXP B1SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type smi1(smi1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type U1(U1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type B1(B1SEXP);
    __result = Rcpp::wrap(significantRcpp(smi1, T1, U1, B1));
    return __result;
END_RCPP
}
// significantRepRcpp
NumericVector significantRepRcpp(SEXP smi1, SEXP T1, SEXP U1, SEXP B1, SEXP perm, SEXP nrep, SEXP nseg, SEXP inperm);
RcppExport SEXP MatrixCorrelation_significantRepRcpp(SEXP smi1SEXP, SEXP T1SEXP, SEXP U1SEXP, SEXP B1SEXP, SEXP permSEXP, SEXP nrepSEXP, SEXP nsegSEXP, SEXP inpermSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type smi1(smi1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type T1(T1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type U1(U1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type perm(permSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nrep(nrepSEXP);
    Rcpp::traits::input_parameter< SEXP >::type nseg(nsegSEXP);
    Rcpp::traits::input_parameter< SEXP >::type inperm(inpermSEXP);
    __result = Rcpp::wrap(significantRepRcpp(smi1, T1, U1, B1, perm, nrep, nseg, inperm));
    return __result;
END_RCPP
}
