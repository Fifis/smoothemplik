// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// kernelWeightsOneCPP
NumericMatrix kernelWeightsOneCPP(NumericVector x, NumericVector xgrid, double bw, bool gaussian);
RcppExport SEXP _smoothemplik_kernelWeightsOneCPP(SEXP xSEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelWeightsOneCPP(x, xgrid, bw, gaussian));
    return rcpp_result_gen;
END_RCPP
}
// kernelDensityOneCPP
NumericVector kernelDensityOneCPP(NumericVector x, NumericVector xgrid, double bw, bool gaussian);
RcppExport SEXP _smoothemplik_kernelDensityOneCPP(SEXP xSEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelDensityOneCPP(x, xgrid, bw, gaussian));
    return rcpp_result_gen;
END_RCPP
}
// kernelSmoothOneCPP
NumericVector kernelSmoothOneCPP(NumericVector x, NumericVector y, NumericVector xgrid, double bw, bool gaussian, bool LOO);
RcppExport SEXP _smoothemplik_kernelSmoothOneCPP(SEXP xSEXP, SEXP ySEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP, SEXP LOOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< double >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    Rcpp::traits::input_parameter< bool >::type LOO(LOOSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelSmoothOneCPP(x, y, xgrid, bw, gaussian, LOO));
    return rcpp_result_gen;
END_RCPP
}
// kernelWeightsMultiCPP
NumericMatrix kernelWeightsMultiCPP(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian);
RcppExport SEXP _smoothemplik_kernelWeightsMultiCPP(SEXP xSEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelWeightsMultiCPP(x, xgrid, bw, gaussian));
    return rcpp_result_gen;
END_RCPP
}
// kernelDensityMultiCPPold
NumericVector kernelDensityMultiCPPold(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian);
RcppExport SEXP _smoothemplik_kernelDensityMultiCPPold(SEXP xSEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelDensityMultiCPPold(x, xgrid, bw, gaussian));
    return rcpp_result_gen;
END_RCPP
}
// kernelDensityMultiCPP
NumericVector kernelDensityMultiCPP(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian);
RcppExport SEXP _smoothemplik_kernelDensityMultiCPP(SEXP xSEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelDensityMultiCPP(x, xgrid, bw, gaussian));
    return rcpp_result_gen;
END_RCPP
}
// kernelSmoothMultiCPP
NumericVector kernelSmoothMultiCPP(NumericMatrix x, NumericVector y, NumericMatrix xgrid, NumericVector bw, bool gaussian, bool LOO);
RcppExport SEXP _smoothemplik_kernelSmoothMultiCPP(SEXP xSEXP, SEXP ySEXP, SEXP xgridSEXP, SEXP bwSEXP, SEXP gaussianSEXP, SEXP LOOSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xgrid(xgridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type bw(bwSEXP);
    Rcpp::traits::input_parameter< bool >::type gaussian(gaussianSEXP);
    Rcpp::traits::input_parameter< bool >::type LOO(LOOSEXP);
    rcpp_result_gen = Rcpp::wrap(kernelSmoothMultiCPP(x, y, xgrid, bw, gaussian, LOO));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_smoothemplik_kernelWeightsOneCPP", (DL_FUNC) &_smoothemplik_kernelWeightsOneCPP, 4},
    {"_smoothemplik_kernelDensityOneCPP", (DL_FUNC) &_smoothemplik_kernelDensityOneCPP, 4},
    {"_smoothemplik_kernelSmoothOneCPP", (DL_FUNC) &_smoothemplik_kernelSmoothOneCPP, 6},
    {"_smoothemplik_kernelWeightsMultiCPP", (DL_FUNC) &_smoothemplik_kernelWeightsMultiCPP, 4},
    {"_smoothemplik_kernelDensityMultiCPPold", (DL_FUNC) &_smoothemplik_kernelDensityMultiCPPold, 4},
    {"_smoothemplik_kernelDensityMultiCPP", (DL_FUNC) &_smoothemplik_kernelDensityMultiCPP, 4},
    {"_smoothemplik_kernelSmoothMultiCPP", (DL_FUNC) &_smoothemplik_kernelSmoothMultiCPP, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_smoothemplik(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
