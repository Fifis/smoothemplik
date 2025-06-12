#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// Forward declarations in other functions
SEXP logTaylorCPP(const NumericVector& x, NumericVector lower, NumericVector upper,
                  IntegerVector der, int order);

NumericVector svdlmCPP(const arma::mat& X, const arma::vec& y, double rel_tol = 1e-9, double abs_tol = 1e-100);

List dampedNewtonCPP(Function fn, NumericVector par, double thresh, int itermax,
                     bool verbose, double alpha, double beta, double backeps);

// Taylor-expanded log-likelihood and its derivatives
List wEL(const arma::vec& lambda,
         const arma::mat& Z,         // n*d centred data
         const arma::vec& ct,        // n weights
         const NumericVector& lower,  // n
         const NumericVector& upper,  // n
         int taylorOrd) {
  IntegerVector ders = IntegerVector::create(0, 1, 2);
  arma::vec s = 1.0 + Z * lambda;
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
static NumericVector  g_lower;
static NumericVector  g_upper;
static int            g_order;

// Wrapper for the optimiser
static SEXP wELlambda(SEXP lambdaSEXP) {
  NumericVector lamR(lambdaSEXP);
  arma::vec lambda(lamR.begin(), lamR.size(), false);
  return wEL(lambda, g_Z, g_ct, g_lower, g_upper, g_order);
}

// [[Rcpp::export]]
List weightedELCPP(NumericMatrix z, NumericVector ct, NumericVector mu, NumericVector lambda_init,
                   bool return_weights, NumericVector lower, NumericVector upper,
                   int order, double weight_tolerance, double thresh, int itermax, bool verbose = false,
                   double alpha = 0.3, double beta = 0.8, double backeps = 0.0) {
  const int n = z.nrow();
  const int d = z.ncol();
  g_Z  = arma::mat(z.begin(), n, d, /*copy_aux_mem =*/ true);

  if (mu.size() == 1) mu = NumericVector(d, mu[0]);
  if (mu.size() != d) stop("The length of mu must match the number of columns in z.");

  // Centre by the hypothesised mean
  for (int j = 0; j < d; ++j) g_Z.col(j) -= mu[j];

  // Observation weights = counts
  if (min(ct) < 0) stop("Negative weights are not allowed.");
  for (double& w : ct) if (w > 0 && w < weight_tolerance) w = 0;
  if (sum(ct) == 0) stop("Total weight must be positive.");
  g_ct = arma::vec(ct.begin(), n,  /*copy_aux_mem =*/ true);

  // Cut-offs
  g_lower = as<NumericVector>(lower);
  if (g_lower.size() == 1) g_lower = NumericVector(n, g_lower[0]);
  g_upper = as<NumericVector>(upper);
  if (g_upper.size() == 1) g_upper = NumericVector(n, g_upper[0]);
  for (int i = 0; i < n; ++i) if (g_lower[i] > g_upper[i]) stop("weightedELCPP: lower > upper");

  g_order = order;   // Taylor order  (>=4)

  // Decide the starting lambda
  arma::vec lam = arma::vec(lambda_init.begin(), d);

  if (arma::any(lam != 0.0)) {  // user supplied non-zero lambda
    arma::vec zerolam(d, fill::zeros);
    IntegerVector der0 = IntegerVector::create(0);
    // Testing 0: 1+lambda'Z = 1
    NumericVector one_l0(n, 1.0);
    arma::vec one_lambdaZ = 1.0 + g_Z * lam;
    NumericVector one_lz(one_lambdaZ.begin(), one_lambdaZ.end());
    NumericVector lt0 = as<NumericVector>( logTaylorCPP(one_l0, g_lower, g_upper, der0, g_order));
    NumericVector lt1 = as<NumericVector>( logTaylorCPP(one_lz, g_lower, g_upper, der0, g_order));
    double f0 = -arma::dot(g_ct, arma::vec(lt0.begin(), n, false));
    double f1 = -arma::dot(g_ct, arma::vec(lt1.begin(), n, false));
    if (f0 < f1) lam.zeros();  // minus log-likelihood is lower -- start from 0 instead
  }

  // Damped Newton optimisation -- key line
  Rcpp::InternalFunction objFunInternal(&wELlambda);
  Rcpp::Function objFun( (SEXP)objFunInternal );
  List opt = dampedNewtonCPP(objFun, NumericVector(lam.begin(), lam.end()), thresh,
                             itermax, verbose, alpha, beta, backeps);

  double logelr = opt["value"];
  NumericVector par = opt["par"];
  int it  = opt["counts"];

  // Probabilities if asked for
  SEXP wts = R_NilValue;
  if (return_weights) {
    arma::vec lamA(par.begin(), d, false);
    arma::vec wtsA = (g_ct / sum(g_ct)) / (1.0 + g_Z * lamA);
    wts = NumericVector(wtsA.begin(), wtsA.end());
  }

  return List::create(
    _["logelr"] = logelr, _["lam"] = par, _["wts"] = wts,
    _["exitcode"] = opt["convergence"], _["iter"] = it,
    _["ndec"] = opt["ndec"], _["gradnorm"]  = opt["gradnorm"]
  );
}
