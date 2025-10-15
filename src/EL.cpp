#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// Forward declarations in other functions
SEXP logTaylorCPP(const NumericVector& x, NumericVector lower, NumericVector upper,
                  IntegerVector der, int order);

NumericVector svdlmCPP(const arma::mat& X, const arma::vec& y, double rel_tol = 1e-9, double abs_tol = 1e-100);

List dampedNewtonCPP(Function fn, NumericVector par, double thresh, int itermax,
                     bool verbose, double alpha, double beta, double backeps,
                     double grad_tol, double step_tol, double f_tol, int stallmax);

// Taylor-expanded log-likelihood and its derivatives
List wEL(const arma::vec& lambda,
         const arma::mat& Z,         // n*d centred data
         const arma::vec& ct,        // n weights
         const double shift,
         const NumericVector& lower,  // n
         const NumericVector& upper,  // n
         const int taylorOrd) {
  IntegerVector ders = IntegerVector::create(0, 1, 2);
  arma::vec s = 1.0 + Z * lambda + shift;
  NumericVector sR(s.begin(), s.end());

  NumericMatrix LT = as<NumericMatrix>(logTaylorCPP(sR, lower, upper, ders, taylorOrd));
  const arma::mat LTm( LT.begin(), LT.nrow(), LT.ncol(), /*copy_aux_mem*/ false );
  arma::vec v0 = -LTm.col(0);  // minus log and its derivatives
  arma::vec v1 = -LTm.col(1);  // -1/(1+l'Z)
  arma::vec v2 = -LTm.col(2);  //  1/(1+l'Z)^2

  double fn     = dot(ct, v0);                  // objective
  arma::vec g   = Z.t() * (ct % v1);            // gradient
  arma::mat h   = Z.t() * (Z.each_col() % (ct % v2));   // Hessian

  NumericVector gR(g.begin(), g.end());
  NumericMatrix hR(h.n_rows, h.n_cols);
  std::copy(h.begin(), h.end(), hR.begin());

  return List::create(_["fn"] = fn, _["gradient"] = gR, _["Hessian"]  = hR);
}


// Shared storage of file-wide static variables because the optimiser still expects
// a function of one argument (parameter vector)
static arma::mat      g_Z;
static arma::vec      g_ct;
static double         g_shift;
static NumericVector  g_lower;
static NumericVector  g_upper;
static int            g_order;

// Wrapper for the optimiser
static SEXP wELlambda(SEXP lambdaSEXP) {
  NumericVector lamR(lambdaSEXP);
  arma::vec lambda(lamR.begin(), lamR.size(), false);
  return wEL(lambda, g_Z, g_ct, g_shift, g_lower, g_upper, g_order);
}

// [[Rcpp::export]]
List ELCPP(NumericMatrix z, NumericVector ct, NumericVector mu, double shift,
           NumericVector lambda_init, bool return_weights, NumericVector lower, NumericVector upper,
           int order, double weight_tolerance, bool deriv = false,
           double thresh = 1e-16, int itermax = 100, bool verbose = false,
           double alpha = 0.3, double beta = 0.8, double backeps = 0.0,
           double grad_tol = 1e-12, double step_tol = 1e-12,
           double f_tol = 1e-14, int stallmax = 5) {
  const int n = z.nrow();
  const int d = z.ncol();

  // Computing the weighted mean centring -- used only for radial reduction
  arma::vec v(d);
  if (deriv) {
    double ct_sum = 0.0;
    for (int i = 0; i < n; ++i) ct_sum += ct[i];  // Positivity had been checked before
    NumericVector wm(d, 0.0);  // Weighted mean by row
    for (int i = 0; i < n; ++i) {
      const double w = ct[i];
      if (w > 0) {
        for (int j = 0; j < d; ++j) wm[j] += w * z(i, j);
      }
    }
    for (int j = 0; j < d; ++j) wm[j] /= ct_sum;
    // Direction of the ray: v = (mu - wm) / ||mu - wm||
    double nn = 0.0;
    for (int j = 0; j < d; ++j) {
      double diff = mu[j] - wm[j];
      v[j] = diff;
      nn += diff * diff;
      }
    if (nn > 0.0) v /= std::sqrt(nn); else v.zeros();   // If mu is wm, v=0
    // Rcpp::Rcout << "Weighted mean wm=" << wm << std::endl;
    // Rcpp::Rcout << "Vector v=" << v << std::endl;
  }

  g_Z  = arma::mat(z.begin(), n, d, /*copy_aux_mem =*/ true);  // Converting to Armadillo objects

  if (mu.size() == 1) mu = NumericVector(d, mu[0]);
  if (mu.size() != d) stop("ELCPP: The length of mu must match the number of columns in z.");

  // Centre by the hypothesised mean
  for (int j = 0; j < d; ++j) g_Z.col(j) -= mu[j];

  // Observation weights = counts
  if (min(ct) < 0) stop("ELCPP: Negative weights are not allowed.");
  for (double& w : ct) if (w > 0 && w < weight_tolerance) w = 0;
  if (sum(ct) == 0) stop("ELCPP: Total weight must be positive.");
  g_ct = arma::vec(ct.begin(), n,  /*copy_aux_mem =*/ true);

  g_shift = shift;

  // Cut-offs
  g_lower = as<NumericVector>(lower);
  if (g_lower.size() == 1) g_lower = NumericVector(n, g_lower[0]);
  g_upper = as<NumericVector>(upper);
  if (g_upper.size() == 1) g_upper = NumericVector(n, g_upper[0]);
  for (int i = 0; i < n; ++i) if (g_lower[i] > g_upper[i]) stop("ELCPP: lower > upper");

  g_order = order;   // Taylor order  (>=4)

  // Rcpp::Rcout << "Set the main variables" << std::endl;

  // Decide the starting lambda
  arma::vec lam = arma::vec(lambda_init.begin(), d);

  if (arma::any(lam != 0.0)) {  // Does user-supplied non-zero lambda yield a lower objective?
    arma::vec zerolam(d, fill::zeros);
    IntegerVector der0 = IntegerVector::create(0);
    // Testing 0: 1+lambda'Z = 1
    NumericVector one_l0(n, 1.0);
    arma::vec one_lambdaZ = 1.0 + g_Z * lam + shift;
    NumericVector one_lz(one_lambdaZ.begin(), one_lambdaZ.end());
    NumericVector lt0 = as<NumericVector>( logTaylorCPP(one_l0, g_lower, g_upper, der0, g_order));
    NumericVector lt1 = as<NumericVector>( logTaylorCPP(one_lz, g_lower, g_upper, der0, g_order));
    double f0 = -arma::dot(g_ct, arma::vec(lt0.begin(), n, false));
    double f1 = -arma::dot(g_ct, arma::vec(lt1.begin(), n, false));
    if (f0 < f1) lam.zeros();  // minus log-likelihood is lower -- start from 0 instead
  }

  // Rcpp::Rcout << "Chose lambda" << std::endl;

  // Damped Newton in lambda -- THE HEART OF EL
  Rcpp::InternalFunction objFunInternal(&wELlambda);
  Rcpp::Function objFun( (SEXP)objFunInternal );
  List opt = dampedNewtonCPP(objFun, NumericVector(lam.begin(), lam.end()), thresh,
                             itermax, verbose, alpha, beta, backeps,
                             grad_tol, step_tol, f_tol, stallmax);

  // Rcpp::Rcout << "Found lambda" << std::endl;

  double logelr = opt["value"];
  NumericVector par = opt["par"];  // This stays a NumericVector
  int it  = opt["counts"];
  arma::vec lambda(par.begin(), d, /*copy_aux_mem=*/ false);  // Here is an ARMA copy

  // Probabilities if asked for
  SEXP wts = R_NilValue;
  if (return_weights) {
    arma::vec wtsA = (g_ct / sum(g_ct)) / (1.0 + g_Z * lambda + g_shift);
    wts = NumericVector(wtsA.begin(), wtsA.end());
  }

  // Derivatives in the direction of mu if asked for (log-ELR scale)
  SEXP derivs = R_NilValue;
  if (deriv) {
    // Rcpp::Rcout << "v direction" << v << std::endl;
    arma::vec u      = 1.0 + g_Z * lambda + g_shift; // n x 1
    arma::vec invu   = 1.0 / u;
    arma::vec invu2  = invu % invu;
    arma::vec winvu2 = g_ct % invu2;  // weights / u^2

    const double S0 = arma::dot(g_ct, invu);
    const double S1 = arma::dot(g_ct, invu2);
    arma::vec  T1  = g_Z.t() * winvu2;                      // d x 1
    arma::mat  T2   = g_Z.t() * (g_Z.each_col() % winvu2);  // d x d

    arma::mat T2inv;  // Stable inverse
    bool ok = false;
    try { T2inv = arma::inv_sympd(T2); ok = true; } catch(...) { ok = false; }
    if (!ok) T2inv = arma::pinv(T2);

    const double lamv = arma::dot(lambda, v);
    arma::vec wv = S0 * v - lamv * T1;

    // log-ELR directional derivatives:
    // f'_v  =  S0 * (lambda' v)
    // f''_v = (wv' T2^{-1} wv) - (lambda' v)^2 S1
    double first = 0.0;
    if (arma::norm(v, 2) != 0.0) {  // mu == wm: zero direction
      first  = S0 * lamv;
    }
    double second = -(arma::as_scalar( wv.t() * T2inv * wv )) + (lamv * lamv) * S1;

    derivs = NumericVector::create(first, second);
  }

  // Final sanity check: if the convex-hull condition is violated, the procedure will never converge
  if (it >= itermax) {
    logelr = R_NegInf; // The exit code cannot be zero here
  }

  return List::create(
    _["logelr"] = logelr, _["lam"] = par, _["wts"] = wts,
    _["deriv"] = derivs,
    _["exitcode"] = opt["convergence"], _["iter"] = it,
    _["ndec"] = opt["ndec"], _["gradnorm"]  = opt["gradnorm"]
  );
}
