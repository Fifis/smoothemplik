#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// svdlm: least-squares coefficients via SVD with Moore--Penrose inverse
// [[Rcpp::export]]
NumericVector svdlmCPP(const arma::mat& x, const arma::vec& y,
                        double rel_tol = 1e-9, double abs_tol = 1e-100) {
  if (x.n_rows != y.n_elem) stop("The number of rows in x must be equal to the length of y.");
  arma::mat U, V;
  arma::vec s;
  // economy-size SVD (U: m×r,  V: n×r,  s: r)
  bool ok = arma::svd_econ(U, s, V, x, "both", "std");
  if (!ok) stop("Armadillo SVD failed");

  double cut_off = std::max(rel_tol * s.max(), abs_tol);
  arma::vec s_inv = 1.0 / s;
  s_inv.elem( arma::find(s < cut_off) ).zeros();

  // Moore–Penrose pseudo-inverse  X+ = V * Sigma^-1 * U'
  // Dimensions:  (n.r) x (r.r) x (r.m)  ->  n.m
  //------------------------------------------------------------------
  arma::mat X_plus = V * arma::diagmat(s_inv) * U.t();
  arma::vec coef   = X_plus * y;

  return Rcpp::NumericVector(coef.begin(), coef.end());
}
