#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declaration from svdlm.cpp
NumericVector svdlmCPP(const arma::mat& x, const arma::vec& y,
                       double rel_tol = 1e-9, double abs_tol = 1e-100);


// [[Rcpp::export]]
List dampedNewtonCPP(Function fn, NumericVector par,
                  double thresh = 1e-16, int itermax = 100, bool verbose = false,
                  double alpha = 0.3, double beta = 0.8, double backeps = 0.0) {
  const int d = par.size();
  arma::vec x(par.begin(), d);

  // Adapter that evaluates fn and returns (f, g, H)
  auto eval_fn = [&](const arma::vec& xx) {
    NumericVector xR(xx.begin(), xx.end());      // armadillo -> Rcpp
    List L = fn(xR);                             // call user function
    double f = as<double>(L["fn"]);
    NumericVector gR = L["gradient"];
    arma::vec g(gR.begin(), gR.size(), true);
    NumericMatrix Hnm = L["Hessian"];
    arma::mat H(Hnm.begin(), Hnm.nrow(), Hnm.ncol(), true);
    return std::tuple<double, arma::vec, arma::mat>(f, g, H);
  };

  // Initial state
  double f; arma::vec g; arma::mat h;
  std::tie(f, g, h) = eval_fn(x);

  bool converged = false;
  int iter = 0;
  double ndec = NA_REAL;
  double gradnorm = norm(g, 2);

  // Main loop
  while (!converged && iter < itermax) {
    ++iter;

    // Newton step  (fallback to pseudo-inverse if singular)
    arma::vec rhs = -g;
    NumericVector stepR = svdlmCPP(h, rhs);
    arma::vec step(stepR.begin(), d, false);

    // Newton decrement & gradient norm
    ndec      = std::sqrt( dot(g, -step) );
    gradnorm  = norm(g, 2);

    // Back-tracking line search
    double t = 1.0;
    double fnew; arma::vec gnew; arma::mat hnew;
    while (true) {
      arma::vec x_trial = x + t * step;
      std::tie(fnew, gnew, hnew) = eval_fn(x_trial);

      if (std::isfinite(fnew) && gnew.is_finite() && hnew.is_finite()) {
        double armijo = f + alpha * t * dot(g, step) + backeps;
        if (fnew <= armijo || t < 1e-4) {  // Accept
          // If all is good or all is bad, update the point
          x = x_trial;  f = fnew;  g = gnew;  h = hnew;
          // Could not find a satisfactory step -- fall back to the current point and leave the line search
          if (t < 1e-4 && verbose) Rcout << "Line search unsuccessful,  accepting a tiny step.\n";
          break;
        }
      }

      t *= beta;  // Shrink
    }

    if (verbose) Rcout << "Iter " << iter << "  f=" << f << "  |grad|=" <<gradnorm <<
      "  decr=" << ndec << "  shrink=" << t << "\n";

    converged = (ndec * ndec <= thresh);
  }

  return List::create(
    _["par"]       = NumericVector(x.begin(), x.end()),
    _["value"]     = f, _["counts"]    = iter, _["convergence"] = converged ? 0 : 1,
    _["ndec"]      = ndec, _["gradnorm"]  = gradnorm
  );
}
