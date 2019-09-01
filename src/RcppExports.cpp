// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bias2
NumericVector bias2(NumericVector lambda, NumericVector mu, double alpha, std::string filt);
RcppExport SEXP _rlystop_bias2(SEXP lambdaSEXP, SEXP muSEXP, SEXP alphaSEXP, SEXP filtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< std::string >::type filt(filtSEXP);
    rcpp_result_gen = Rcpp::wrap(bias2(lambda, mu, alpha, filt));
    return rcpp_result_gen;
END_RCPP
}
// variance
NumericVector variance(NumericVector lambda, double delta, double alpha, std::string filt);
RcppExport SEXP _rlystop_variance(SEXP lambdaSEXP, SEXP deltaSEXP, SEXP alphaSEXP, SEXP filtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< std::string >::type filt(filtSEXP);
    rcpp_result_gen = Rcpp::wrap(variance(lambda, delta, alpha, filt));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_rlystop_bias2", (DL_FUNC) &_rlystop_bias2, 4},
    {"_rlystop_variance", (DL_FUNC) &_rlystop_variance, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_rlystop(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
