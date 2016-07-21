// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cmr
arma::mat cmr(arma::mat psi, arma::mat X, int n, arma::mat K, int reps, arma::vec margin);
RcppExport SEXP mumocomp_cmr(SEXP psiSEXP, SEXP XSEXP, SEXP nSEXP, SEXP KSEXP, SEXP repsSEXP, SEXP marginSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type reps(repsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type margin(marginSEXP);
    __result = Rcpp::wrap(cmr(psi, X, n, K, reps, margin));
    return __result;
END_RCPP
}
