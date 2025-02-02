// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// brentMin
Rcpp::List brentMin(Function f, NumericVector interval, double lower, double upper, double tol, int maxiter, int trace);
RcppExport SEXP _smoothemplik_brentMin(SEXP fSEXP, SEXP intervalSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(brentMin(f, interval, lower, upper, tol, maxiter, trace));
    return rcpp_result_gen;
END_RCPP
}
// brentZero
List brentZero(Function f, NumericVector interval, double lower, double upper, Nullable<NumericVector> f_lower, Nullable<NumericVector> f_upper, std::string extendInt, bool check_conv, double tol, int maxiter, int trace);
RcppExport SEXP _smoothemplik_brentZero(SEXP fSEXP, SEXP intervalSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP f_lowerSEXP, SEXP f_upperSEXP, SEXP extendIntSEXP, SEXP check_convSEXP, SEXP tolSEXP, SEXP maxiterSEXP, SEXP traceSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type interval(intervalSEXP);
    Rcpp::traits::input_parameter< double >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< double >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type f_lower(f_lowerSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type f_upper(f_upperSEXP);
    Rcpp::traits::input_parameter< std::string >::type extendInt(extendIntSEXP);
    Rcpp::traits::input_parameter< bool >::type check_conv(check_convSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type trace(traceSEXP);
    rcpp_result_gen = Rcpp::wrap(brentZero(f, interval, lower, upper, f_lower, f_upper, extendInt, check_conv, tol, maxiter, trace));
    return rcpp_result_gen;
END_RCPP
}
// kernelFunCPP
arma::vec kernelFunCPP(arma::vec x, std::string kernel, int order, bool convolution);
RcppExport SEXP _smoothemplik_kernelFunCPP(SEXP xSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelFunCPP(x, kernel, order, convolution));
    return rcpp_result_gen;
END_RCPP
}
// kernelWeightsOneCPP
arma::mat kernelWeightsOneCPP(arma::vec x, arma::vec xout, double bw, std::string kernel, int order, bool convolution);
RcppExport SEXP _smoothemplik_kernelWeightsOneCPP(SEXP xSEXP, SEXP xoutSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelWeightsOneCPP(x, xout, bw, kernel, order, convolution));
    return rcpp_result_gen;
END_RCPP
}
// sparseKernelWeightsOneCPP
arma::sp_mat sparseKernelWeightsOneCPP(arma::vec x, arma::vec xout, double bw, std::string kernel, int order, bool convolution);
RcppExport SEXP _smoothemplik_sparseKernelWeightsOneCPP(SEXP xSEXP, SEXP xoutSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseKernelWeightsOneCPP(x, xout, bw, kernel, order, convolution));
    return rcpp_result_gen;
END_RCPP
}
// kernelWeightsCPP
arma::mat kernelWeightsCPP(arma::mat x, arma::mat xout, arma::vec bw, std::string kernel, int order, bool convolution);
RcppExport SEXP _smoothemplik_kernelWeightsCPP(SEXP xSEXP, SEXP xoutSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelWeightsCPP(x, xout, bw, kernel, order, convolution));
    return rcpp_result_gen;
END_RCPP
}
// sparseKernelWeightsCPP
arma::sp_mat sparseKernelWeightsCPP(arma::mat x, arma::mat xout, arma::vec bw, std::string kernel, int order, bool convolution);
RcppExport SEXP _smoothemplik_sparseKernelWeightsCPP(SEXP xSEXP, SEXP xoutSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseKernelWeightsCPP(x, xout, bw, kernel, order, convolution));
    return rcpp_result_gen;
END_RCPP
}
// kernelDensityCPP
NumericVector kernelDensityCPP(arma::mat x, arma::mat xout, arma::vec weights, arma::vec bw, std::string kernel, int order, bool convolution, int chunks);
RcppExport SEXP _smoothemplik_kernelDensityCPP(SEXP xSEXP, SEXP xoutSEXP, SEXP weightsSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP convolutionSEXP, SEXP chunksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    Rcpp::traits::input_parameter< int >::type chunks(chunksSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelDensityCPP(x, xout, weights, bw, kernel, order, convolution, chunks));
    return rcpp_result_gen;
END_RCPP
}
// kernelSmoothCPP
NumericVector kernelSmoothCPP(arma::mat x, arma::vec y, arma::mat xout, arma::vec weights, arma::vec bw, std::string kernel, int order, bool LOO, bool convolution, int chunks);
RcppExport SEXP _smoothemplik_kernelSmoothCPP(SEXP xSEXP, SEXP ySEXP, SEXP xoutSEXP, SEXP weightsSEXP, SEXP bwSEXP, SEXP kernelSEXP, SEXP orderSEXP, SEXP LOOSEXP, SEXP convolutionSEXP, SEXP chunksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xout(xoutSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< std::string >::type kernel(kernelSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type LOO(LOOSEXP);
    Rcpp::traits::input_parameter< bool >::type convolution(convolutionSEXP);
    Rcpp::traits::input_parameter< int >::type chunks(chunksSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelSmoothCPP(x, y, xout, weights, bw, kernel, order, LOO, convolution, chunks));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_smoothemplik_brentMin", (DL_FUNC) &_smoothemplik_brentMin, 7},
    {"_smoothemplik_brentZero", (DL_FUNC) &_smoothemplik_brentZero, 11},
    {"_smoothemplik_kernelFunCPP", (DL_FUNC) &_smoothemplik_kernelFunCPP, 4},
    {"_smoothemplik_kernelWeightsOneCPP", (DL_FUNC) &_smoothemplik_kernelWeightsOneCPP, 6},
    {"_smoothemplik_sparseKernelWeightsOneCPP", (DL_FUNC) &_smoothemplik_sparseKernelWeightsOneCPP, 6},
    {"_smoothemplik_kernelWeightsCPP", (DL_FUNC) &_smoothemplik_kernelWeightsCPP, 6},
    {"_smoothemplik_sparseKernelWeightsCPP", (DL_FUNC) &_smoothemplik_sparseKernelWeightsCPP, 6},
    {"_smoothemplik_kernelDensityCPP", (DL_FUNC) &_smoothemplik_kernelDensityCPP, 8},
    {"_smoothemplik_kernelSmoothCPP", (DL_FUNC) &_smoothemplik_kernelSmoothCPP, 10},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_smoothemplik(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
