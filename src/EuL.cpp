#include <RcppArmadillo.h>

using namespace Rcpp;

// Weighted Euclidean likelihood: the C++ implementation is ~4 times faster than
// the R one

// [[Rcpp::export]]
List EuLCPP(const arma::mat& z, arma::vec mu,
            arma::vec ct, arma::vec vt, arma::vec shift,
            const double n_orig, const double weight_tolerance,
            const double trunc_to = 0.0, const bool SEL = true,
            const bool return_weights = false, const bool verbose = false,
            const bool chull_diag = false)
{
  const double me = std::numeric_limits<double>::epsilon();

  // These will be mutated
  arma::mat zz = z;
  const std::size_t n = zz.n_rows;
  const std::size_t k = zz.n_cols;

  if (mu.n_elem == 1) mu = arma::vec(k, arma::fill::value(mu[0]));
  if (mu.n_elem != k) stop("The length of mu must match the number of columns in z.");
  if (ct.n_elem != n) stop("Length of ct must equal nrow(zz).");
  if (vt.n_elem != n) stop("Length of ct must equal nrow(zz).");
  if (!zz.is_finite()) stop("Non-finite observations (NA, NaN, Inf) are not welcome.");
  if (!ct.is_finite()) stop("Non-finite weights (NA, NaN, Inf) are not welcome.");
  if (!vt.is_finite()) stop("Non-finite variance weights (NA, NaN, Inf) are not welcome.");
  if (ct.min() < 0 && verbose) warning("Negative weights are present.");
  if (vt.min() < 0) stop("Negative variance weights are not welcome.");
  if (arma::accu(ct) <= 0.0) stop("The total sum of weights must be positive.");
  if (arma::accu(vt) == 0.0) stop("All variances are zero -- nothing to optimise.");
  if (!SEL) stop("Only 'ct' representing weights (adding up to one) are supported (SEL = TRUE).");

  zz.each_row() -= mu.t(); // z <- z - mu

  // Truncating tiny SEuL weights
  arma::uvec tinyc = arma::find(arma::abs(ct) < weight_tolerance);
  if (!tinyc.empty()) {
    if (verbose) {
      Rcpp::warning("Counts closer to 0 than %1.2e have been replaced with %s",
                    weight_tolerance, (trunc_to == 0.0 ? "0." : "+-trunc.to of appropriate sign."));
    }
    ct.elem(tinyc) = trunc_to * arma::sign(ct.elem(tinyc));
  }
  // Truncating tiny variance weights
  arma::uvec tinyv = arma::find(arma::abs(vt) < weight_tolerance);
  if (!tinyv.empty()) {
    if (verbose) {
      Rcpp::warning("Variance weights closer to 0 than %1.2e have been replaced with %s",
                    weight_tolerance, (trunc_to == 0.0 ? "0." : "+-trunc.to of appropriate sign."));
    }
    vt.elem(tinyv).fill(trunc_to);
  }

  // Discard zero-weight rows to save size
  arma::uvec nonz = arma::find(vt != 0.0);
  zz   = zz.rows(nonz);
  ct   = ct.elem(nonz);
  vt   = vt.elem(nonz);
  shift = shift.elem(nonz);  // TODO: use later
  const std::size_t n_final = nonz.n_elem;

  // Re-normalising for SEL
  if (SEL) {
    ct /= arma::accu(ct);
  }
  vt /= arma::accu(vt);

  // Original sample size for normalisation
  const double N = SEL ? static_cast<double>(n_orig) : arma::accu(ct);
  if (N <= 0) stop("EuLCPP: Total weights after tolerance checks must be positive.");

  // Rcout << "zz: " << zz << "\n";
  // Rcout << "ct: " << ct << "\n";
  // Rcout << "vt: " << vt << "\n";
  // Rcout << "N: " << N << "\n";

  // Default output in case of failure
  double logelr = -arma::datum::inf;
  arma::vec lam(k, arma::fill::value(arma::datum::inf));
  arma::vec wts;
  if (return_weights) wts.zeros(n);
  bool converged   = false;
  int iter        = 1;
  NumericVector bracket = NumericVector::create(R_NegInf, R_PosInf);
  double estim_prec = NA_REAL;
  double f_root = NA_REAL;
  int exitcode = 3;           // Default: failed to invert a matrix

  // Heart of the function: compute lambdas and ELR
  if (n_final < k) {               // Not enough points to solve the system
    exitcode = 2;
  } else {
    arma::vec m = zz.t() * ct;  // Count-weighted average of z
    arma::mat zc = zz.each_row() - m.t();  // Z minus its conditional expectation for mean
    arma::vec mv = zz.t() * vt;  // Variance-weight-weighted average of z for de-meaning
    arma::mat zcv = zz.each_row() - mv.t();  // Z minus its conditional expectation for var
    arma::mat Omega = zcv.t() * (zcv.each_col() % vt);  // Weighted variance estimator

    // Rcout << "m: " << m << "\n";
    // Rcout << "zc: " << zc << "\n";
    // Rcout << "mv: " << mv << "\n";
    // Rcout << "zcv: " << zcv << "\n";
    // Rcout << "vt * zc'zc: " << Omega << "\n";

    bool solver_ok = true;
    try {
      lam = -arma::solve(Omega, m, arma::solve_opts::no_approx);
    } catch (const std::runtime_error& e) {
      solver_ok = false;
    }
    // Rcout << "lambda: " << lam << "\n";

    if (solver_ok) {
      arma::vec  zlam = zcv * lam;  // lambda' (z_j - mv_i)
      arma::vec  wvec = ct + vt % zlam;
      lam /= N;  // Because the true lambda is -1/n inv(Omega) m
      if (!SEL) wvec /= N;
      if (return_weights) wts.elem(nonz) = wvec;
      logelr     = -0.5 * arma::accu(arma::square(wvec - ct));
      converged  = true;
      estim_prec = me;
      // A measure of how close the derivative is to zero
      f_root  = arma::abs(zz.t() * wvec).max();

      exitcode = 0;

      // Convex-hull diagnostic
      // Zero cannot be in the convex hull of a column if its min and max
      // have the same non-zero sign (their product is strictly positive)
      if (chull_diag) {
        arma::rowvec colMax = arma::max(zz, /*dim=*/0);
        arma::rowvec colMin = arma::min(zz, /*dim=*/0);
        if (arma::any((colMin % colMax) > 0.0)) exitcode = 1;
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
