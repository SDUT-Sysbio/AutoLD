// Generated interface for src/em_core.cpp
#include <Rcpp.h>
using namespace Rcpp;

Rcpp::List EM_full_cpp(Rcpp::NumericVector counts, int c_ploidy, double tol, int max_iter);
Rcpp::List EM_partial_cpp(Rcpp::NumericVector counts, int c_ploidy, double tol, int max_iter);

RcppExport SEXP _AutoLD_EM_full_cpp(SEXP countsSEXP, SEXP c_ploidySEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< int >::type c_ploidy(c_ploidySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    return Rcpp::wrap(EM_full_cpp(counts, c_ploidy, tol, max_iter));
END_RCPP
}

RcppExport SEXP _AutoLD_EM_partial_cpp(SEXP countsSEXP, SEXP c_ploidySEXP, SEXP tolSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type counts(countsSEXP);
    Rcpp::traits::input_parameter< int >::type c_ploidy(c_ploidySEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    return Rcpp::wrap(EM_partial_cpp(counts, c_ploidy, tol, max_iter));
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_AutoLD_EM_full_cpp", (DL_FUNC) &_AutoLD_EM_full_cpp, 4},
    {"_AutoLD_EM_partial_cpp", (DL_FUNC) &_AutoLD_EM_partial_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_AutoLD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
