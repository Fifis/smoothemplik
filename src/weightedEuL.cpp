#include <RcppArmadillo.h>

using namespace Rcpp;
using arma::mat;
using arma::vec;
using arma::uvec;
using arma::rowvec;

// weighted Euclidean likelihood: the C++ implementation is ~4 times faster than
// the R one

// [[Rcpp::export]]
Rcpp::List weightedEuLCPP(const arma::mat& z, arma::vec mu,
                          arma::vec ct, arma::vec shift,
                          const double n_orig, double weight_tolerance,
    const double trunc_to = 0.0, const bool SEL = true,
    const bool return_weights = false, const bool verbose = false,
    const bool chull_diag = false)
{
  const double me = std::numeric_limits<double>::epsilon();

  // These will be mutated
  mat zz = z;
  const std::size_t n = zz.n_rows;
  const std::size_t k = zz.n_cols;

  if (mu.n_elem == 1) mu.fill(mu[0]);
  if (mu.n_elem != k) stop("The length of mu must match the number of columns in z.");

  zz.each_row() -= mu.t(); // z <- z - mu

  if (ct.n_elem != n)
    stop("Length of ct must equal nrow(zz).");

  if (!zz.is_finite())
    stop("Non-finite observations (NA, NaN, Inf) are not welcome.");
  if (!ct.is_finite())
    stop("Non-finite weights (NA, NaN, Inf) are not welcome.");

  if (ct.min() < 0 && verbose)
    warning("Negative weights are present.");
  if (arma::accu(ct) <= 0.0)
    stop("The total sum of weights must be positive.");

  if (!SEL) stop("Only 'ct' representing weights (adding up to one) are supported (SEL = TRUE).");

  // Truncating tiny weights
  uvec tiny = arma::find(arma::abs(ct) < weight_tolerance);
  if (!tiny.empty()) {
    if (verbose) {
      Rcpp::warning("Counts closer to 0 than %1.2e have been replaced with %s",
                    weight_tolerance, (trunc_to == 0.0 ? "0." : "+-trunc.to of appropriate sign."));
    }
    ct.elem(tiny) = trunc_to * arma::sign(ct.elem(tiny));
  }

  // Discard zero-weight rows to save size
  uvec nonz = arma::find(ct != 0.0);
  zz   = zz.rows(nonz);
  ct  = ct.elem(nonz);
  shift = shift.elem(nonz);  // TODO: use later
  const std::size_t n_final = nonz.n_elem;

  if (SEL) ct /= arma::accu(ct);

  const double N = SEL ? static_cast<double>(n_orig) : arma::accu(ct);
  if (N <= 0) stop("Total weights after tolerance checks must be positive.");

  // Default output in case of failure
  double logelr      = -arma::datum::inf;
  vec    lam(k, arma::fill::value(arma::datum::inf));
  vec wts;
  if (return_weights) wts = vec(n, arma::fill::zeros);
  bool   converged   = false;
  int    iter        = 1;
  NumericVector bracket = NumericVector::create(R_NegInf, R_PosInf);
  double estim_prec  = NA_REAL;
  double f_root      = NA_REAL;
  int    exitcode    = 3;           // Default: failed to invert a matrix

  // Heart of the function: compute lambdas and ELR
  if (n_final < k) {               // Not enough points to solve the system
    exitcode = 2;
  } else {
    rowvec zbar = arma::mean(zz, 0);
    mat zc = zz.each_row() - zbar;
    vec m = zz.t() * ct;  // Count-weighted average of z
    mat zc_zc = zc.t() * zc;

    // Debug
    // Rcout << "zc: " << zc << "\n";
    // Rcout << "m: " << m << "\n";
    // Rcout << "zc'zc: " << zc_zc << "\n";

    bool solver_ok = true;
    try {
      lam = -arma::solve(zc_zc, m, arma::solve_opts::no_approx);
    } catch (const std::runtime_error& e) {
      solver_ok = false;
    }

    if (solver_ok) {
      vec wvec = ct + zc * lam;
      if (!SEL) wvec /= N;
      if (return_weights) wts.elem(nonz) = wvec;
      logelr     = -0.5 * arma::accu(arma::square(wvec - ct));
      converged  = true;
      estim_prec = me;
      // A measure of how close the derivative is to zero
      vec tmp = zz.t() * wvec;
      f_root  = arma::abs(tmp).max();

      exitcode = 0;

      // convex-hull diagnostic
      if (chull_diag) {
        bool outside = false;
        for (std::size_t j = 0; j < k; ++j) {
          double mn = zz.col(j).min(), mx = zz.col(j).max();
          if (mn * mx > 0) { outside = true; break; }
        }
        if (outside && exitcode == 0) exitcode = 1;
      }
    }
  }

  std::vector<std::string> msgs = {
    "successful convergence, mu is within the convex hull of z",
    "successful convergence, mu is certainly outside the convex hull of z",
    "at least ncol(z) points non-zero weights are needed for identification",
    "could not invert the matrix (probably columns of 'z' are collinear)"
  };
  if (verbose && exitcode > 0) warning(msgs[exitcode]);

  SEXP wts_out = R_NilValue;  // To handle NULL or vector
  if (return_weights) wts_out = Rcpp::NumericVector(wts.begin(), wts.end());

  return List::create(
    _["logelr"] = logelr, _["lam"] = Rcpp::NumericVector(lam.begin(), lam.end()),
    _["wts"] = wts_out,
    _["converged"] = converged, _["iter"] = iter, _["bracket"] = bracket,
    _["estim.prec"] = estim_prec, _["f.root"] = f_root, _["exitcode"] = exitcode,
    _["message"] = msgs[exitcode]
  );
}
