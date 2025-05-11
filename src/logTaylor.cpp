#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Functions related to the logarithm, its derivatives, its Taylor expansions,
// and the derivatives thereof

// Factorial for non-negative integers
inline double factorialCPP(int n) {
  if (n < 0)
    Rcpp::stop("factorial is not defined for negative integers.");

  static constexpr double table[16] = {
    1.0,
    1.0, 2.0, 6.0, 24.0, 120.0,  // 1..5
    720.0, 5040.0, 40320.0, 362880.0, 3628800.0,  // 6..10
    39916800.0, 479001600.0, 6227020800.0, 87178291200.0, 1307674368000.0 // 11..15
  };

  return n <= 15 ? table[n] : std::tgamma(static_cast<double>(n) + 1.0);
}


// Scalar d-th derivative of log(x); d = 0 gives log itself
// dlog(x, d) = (-1)^(d-1) * (d-1)! / x^d   for d ≥ 1
inline double dlog_scalarCPP(double x, int d) {
  if (d == 0) return std::log(x);
  if (d < 0)
    Rcpp::stop("dlog_scalarCPP: derivative order must be non-negative");
  // (d-1)!  --- very small d (<= 10) in practice, so looping is fine
  double fac  = factorialCPP(d - 1);
  double sign = ((d - 1) % 2 == 0) ? 1.0 : -1.0;
  return sign * fac / std::pow(x, d);
}


// Vectorised wrapper around the scalar version above

// [[Rcpp::export]]
NumericVector dlogCPP(const NumericVector& x, int d) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i)
    out[i] = dlog_scalarCPP(x[i], d);
  return out;
}


// Taylor expansion of the logarithm and its derivatives

// [[Rcpp::export]]
NumericVector tlogCPP(NumericVector x,
                      NumericVector a = NumericVector::create(1.0),
                      int k = 4, int d = 0) {
  int n_x = x.size();
  int len_a = a.size();
  if (!(len_a == 1 || len_a == n_x)) {
    stop("'a' must be either length 1 or the same length as 'x'.");
  }
  if (a.size() == 1)
    a = NumericVector(n_x, a[0]);

  if (k < 0 || d < 0)
    stop("The polynomial order 'k' and derivative order 'd' must be non-negative integer scalars.");

  if (d > k) {
    // If derivative order is greater than polynomial order, everything is zero
    // Return a zero vector of the same length as x
    return NumericVector(n_x, 0.0);
  }

  NumericVector out(n_x, 0.0);
  NumericVector xc = (x-a)/a;
  for (int i = 0; i < n_x; i++) {
    double accum = 0.0;

    for (int n = d; n <= k; n++) {
      if (n == d) {
        if (d == 0) {
          // Special case: n == d, lowest order: constant log(a)
          accum += std::log(a[i]);
        } else {
          double gamma_n = factorialCPP(n - 1);     // gamma(n)
          double denom   = std::pow(-a[i], (double)d);
          accum += -gamma_n / denom;
        }
      } else {
        double mult  = (d == 0) ? 1.0 : factorialCPP(n) / factorialCPP(n - d);
        double sign  = ((n - 1) % 2 == 0) ? 1.0 : -1.0;
        double denom = (double)n * std::pow(a[i], (double)d);
        double factor = sign * mult / denom;
        accum += factor * std::pow(xc[i], (double)(n - d));
      }
    }
    out[i] = accum;
  }

  return wrap(out);
}


// The big wrapper to do this expansion for vector inputs and cut-off points
// [[Rcpp::export]]
SEXP logTaylorCPP(const NumericVector&  x,
                  NumericVector lower,
                  NumericVector upper,
                  IntegerVector der  = IntegerVector::create(0),
                  int order  = 4     // pass NA_integer_ for 'use pure log / derivs, no Taylor'
) {
  const int n = x.size();
  const int m = der.size();  // number of derivative orders requested

  // Check lower and upper
  if (lower.size() == 1) lower = NumericVector(n, lower[0]);
  if (upper.size()   == 1) upper   = NumericVector(n, upper[0]);
  if (lower.size() != upper.size())
    stop("`lower` and `upper` must have the same length.");
  for (int i = 0; i < n; ++i) {
    if (lower[i] > upper[i]) stop("Thresholds not ordered: lower[i] must be <= upper[i].");
  }
  for (int d : der)
    if (d < 0) stop("All derivative orders must be non-negative.");

  NumericMatrix out(n, m);

  // Flag: should we *skip* Taylor and use plain log / derivatives only?
  bool noTaylorApprox = Rcpp::traits::is_na<INTSXP>(order) || order < 0;

  NumericVector x1(1), a1(1);  // Small buffers to avoid reallocating each loop

  // Main loop
  for (int j = 0; j < m; ++j) {
    const int d = der[j];  // Fix a column corresponding to a derivative order

    for (int i = 0; i < n; ++i) {   // Fix an observation

      const double xi   = x[i];
      const double loweri = lower[i];
      const double upperi   = upper[i];

      double val;

      if (noTaylorApprox) {  // plain derivative of log
        val = dlog_scalarCPP(xi, d);
      } else if (xi < loweri) { // left Taylor tail
        x1[0] = xi;  a1[0] = loweri;
        val = tlogCPP(x1, a1, order, d)[0];
      } else if (xi > upperi) { // right Taylor tail
        x1[0] = xi;  a1[0] = upperi;
        val = tlogCPP(x1, a1, order, d)[0];
      } else { // middle region: exact derivative
        val = dlog_scalarCPP(xi, d);
      }
      out(i, j) = val;
    }
  }

  // Drop dimensions to vector if requested
  if (m == 1) {
    NumericVector v(n);
    std::copy(out.begin(), out.begin() + n, v.begin());
    return v;
  }

  return out;
}

