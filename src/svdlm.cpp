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

  // if (!ok) stop("Armadillo SVD for linear projection failed.");

  // (1) Trying the vanilla SVD
  if (ok) {  // SVD was successful
    double cut_off = std::max(rel_tol * s.max(), abs_tol);
    arma::vec s_inv = 1.0 / s;
    s_inv.elem( arma::find(s < cut_off) ).zeros();
    // Moore–Penrose pseudo-inverse  X+ = V * Sigma^-1 * U'
    // Dimensions:  (n.r) x (r.r) x (r.m)  ->  n.m
    arma::mat X_plus = V * arma::diagmat(s_inv) * U.t();
    arma::vec coef   = X_plus * y;
    return Rcpp::NumericVector(coef.begin(), coef.end());
  }

  // (2) Try pinv (non-throwing overload first)
  arma::mat X_plus;
  bool ok_pinv = false;
  try {
    ok_pinv = arma::pinv(X_plus, x);     // returns false instead of throwing
  } catch (...) {
    ok_pinv = false;  // defensive: handle the error
  }
  if (!ok_pinv) {  // Some versions may still throw in the returning overload:
    try {
      X_plus = arma::pinv(x);            // This one can throw
      ok_pinv = X_plus.is_finite();
    } catch (...) {
      ok_pinv = false;
    }
  }
  if (ok_pinv) {
    arma::vec coef = X_plus * y;
    if (coef.is_finite()) {
      return Rcpp::NumericVector(coef.begin(), coef.end());
    }
  }

  // If all else fails
  // (3) Levenberg-Marquardt fallback on normal equations:
  //    (X'X + lambda I) = X' y,  lambda = rel_tol * max(diag(X'X)) (>= abs_tol)
  arma::mat XtX = x.t() * x;
  arma::vec Xty = x.t() * y;

  double diag_max = 0.0;
  if (XtX.n_rows > 0) {
    arma::vec d = XtX.diag();
    diag_max = d.is_empty() ? 0.0 : d.max();
  }
  double lambda = std::max(rel_tol * diag_max, abs_tol);

  arma::vec coef_lm;
  bool solved = false;
  // Allow a handful of ridge increases if the system is stubborn
  for (int k = 0; k < 6 && !solved; ++k) {
    arma::mat M = XtX + lambda * arma::eye(XtX.n_rows, XtX.n_cols);
    try {
      coef_lm = arma::solve(M, Xty, arma::solve_opts::fast + arma::solve_opts::likely_sympd);
      solved = coef_lm.is_finite();
    } catch (...) {
      solved = false;
    }
    if (!solved) lambda *= 10.0; // strengthen the ridge and retry
  }

  if (!solved) {
    stop("svdlmCPP: SVD and pinv failed; LM ridge solve also failed.");
  }
  return Rcpp::NumericVector(coef_lm.begin(), coef_lm.end());
}
