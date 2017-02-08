// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpptnm
SEXP rcpptnm(SEXP sigma, SEXP a, SEXP mu, SEXP sig, SEXP ssig, SEXP low, SEXP high);
RcppExport SEXP CNTTN_rcpptnm(SEXP sigmaSEXP, SEXP aSEXP, SEXP muSEXP, SEXP sigSEXP, SEXP ssigSEXP, SEXP lowSEXP, SEXP highSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type a(aSEXP);
    Rcpp::traits::input_parameter< SEXP >::type mu(muSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sig(sigSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ssig(ssigSEXP);
    Rcpp::traits::input_parameter< SEXP >::type low(lowSEXP);
    Rcpp::traits::input_parameter< SEXP >::type high(highSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpptnm(sigma, a, mu, sig, ssig, low, high));
    return rcpp_result_gen;
END_RCPP
}
// rcppcomdel
SEXP rcppcomdel(SEXP W, SEXP U, SEXP H);
RcppExport SEXP CNTTN_rcppcomdel(SEXP WSEXP, SEXP USEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type W(WSEXP);
    Rcpp::traits::input_parameter< SEXP >::type U(USEXP);
    Rcpp::traits::input_parameter< SEXP >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppcomdel(W, U, H));
    return rcpp_result_gen;
END_RCPP
}
// rcppng
SEXP rcppng(SEXP sigma, SEXP sigw, SEXP isigt, SEXP invL, SEXP muW, SEXP rngN);
RcppExport SEXP CNTTN_rcppng(SEXP sigmaSEXP, SEXP sigwSEXP, SEXP isigtSEXP, SEXP invLSEXP, SEXP muWSEXP, SEXP rngNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type sigw(sigwSEXP);
    Rcpp::traits::input_parameter< SEXP >::type isigt(isigtSEXP);
    Rcpp::traits::input_parameter< SEXP >::type invL(invLSEXP);
    Rcpp::traits::input_parameter< SEXP >::type muW(muWSEXP);
    Rcpp::traits::input_parameter< SEXP >::type rngN(rngNSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppng(sigma, sigw, isigt, invL, muW, rngN));
    return rcpp_result_gen;
END_RCPP
}
