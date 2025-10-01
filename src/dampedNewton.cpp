#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declaration from svdlm.cpp
NumericVector svdlmCPP(const arma::mat& x, const arma::vec& y,
                       double rel_tol = 1e-9, double abs_tol = 1e-100);

inline bool finite_all(const arma::vec& v) { return v.is_finite(); }
inline bool finite_all(const arma::mat& M) { return M.is_finite(); }

// [[Rcpp::export]]
List dampedNewtonCPP(Function fn, NumericVector par,
                  double thresh = 1e-16, int itermax = 100, bool verbose = false,
                  double alpha = 0.3, double beta = 0.8, double backeps = 0.0,
                  double grad_tol = 1e-12, double step_tol = 1e-12,
                  double f_tol = 1e-14, int stallmax = 5) {
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

  double f_prev = f;
  arma::vec x_prev = x;
  int stall = 0;

  // Main loop
  while (!converged && iter < itermax) {
    ++iter;

    // Newton step  (fallback to pseudo-inverse if singular)
    arma::vec rhs = -g;
    NumericVector stepR = svdlmCPP(h, rhs);
    arma::vec step(stepR.begin(), d, false);

    // Guarantee descent; fall back to steepest descent if needed
    double gTp = arma::dot(g, step);
    if (!step.is_finite() || !std::isfinite(gTp) || gTp >= 0.0) {
      step = -g;                      // steepest descent
      gTp  = -arma::dot(g, g);        // = -||g||^2
    }
    // Newton decrement & gradient norm (safe)
    ndec      = std::sqrt(std::max(0.0, -arma::dot(g, step)));
    gradnorm  = norm(g, 2);

    // Back-tracking line search
    double t = 1.0;
    double fnew; arma::vec gnew; arma::mat hnew;

    // Guards against infinite backtracking
    const double t_min = 1e-12;
    const int    max_backtracks = 200;
    int          n_back = 0;

    while (true) {
      arma::vec x_trial = x + t * step;
      std::tie(fnew, gnew, hnew) = eval_fn(x_trial);

      bool finite_trial = std::isfinite(fnew) && gnew.is_finite() && hnew.is_finite();
      if (finite_trial) {
        double armijo = f + alpha * t * arma::dot(g, step) + backeps;
        double df = f - fnew; // > 0 if decrease
        if (fnew <= armijo || df > 1e-16 || t < 1e-6) {
          // Accept trial update (either Armijo satisfied, or tiny measurable decrease, or tiny step)
          x = x_trial;  f = fnew;  g = gnew;  h = hnew;
          if (t < 1e-6 && verbose) Rcout << "Line search: tiny step " << t << " accepted.\n";
          break;
        }
      }

      t *= beta;  // Shrink the step
      ++n_back;

      // Global bail-out even when the trial is non-finite (may happen if logTaylor does not use Taylor)
      if (t < t_min || n_back >= max_backtracks) {
        if (verbose) Rcout << "Line search failed (non-finite or no Armijo). Proceeding without update.\n";
        // Plan B: take a tiny gradient step instead of Newton
        arma::vec gd_step = -g;
        x += t_min * gd_step / std::max(1.0, norm(gd_step, 2));
        std::tie(f, g, h) = eval_fn(x);
        break;
      }
    }

    double rel_step = arma::norm(x - x_prev, 2) / std::max(1.0, arma::norm(x, 2));
    double rel_f    = std::abs(f - f_prev) / std::max(1.0, std::abs(f));

    if (verbose) Rcout << "Iter " << iter << "  f=" << f << "  |grad|=" << gradnorm
                       << "  decr=" << ndec << "  shrink=" << t
                       << "  g'p=" << arma::dot(g, step)
                       << "  trace(H)=" << arma::trace(h)
                       << "  rel_step=" << rel_step
                       << "  rel_f=" << rel_f
                       << (gTp >= 0.0 ? "  (steepest)" : "")  // flag if fallback was used
                       << "\n";

    // Stall counting
    if (rel_f < f_tol && rel_step < step_tol) ++stall; else stall = 0;
    // Remember previous state
    f_prev = f;
    x_prev = x;

    // Multiple stop criteria
    converged = (ndec * ndec <= thresh) || (gradnorm <= grad_tol)  || (rel_step  <= step_tol) ||
      (rel_f     <= f_tol)    || (stall     >= stallmax);

  }

  return List::create(
    _["par"]       = NumericVector(x.begin(), x.end()),
    _["value"]     = f, _["counts"]    = iter, _["convergence"] = converged ? 0 : 1,
    _["ndec"]      = ndec, _["gradnorm"]  = gradnorm
  );
}
