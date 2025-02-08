#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Pre-computed factorials from 0! through 10!
static const double factorial_table[11] = {
  1.0, 1.0, 2.0, 6.0, 24.0, 120.0, // 0! ... 5!
  720.0, 5040.0, 40320.0, 362880.0, 3628800.0 // 6! ... 10!
};


// Factorial for non-negative integers
double factorialC(int m) {
  if (m < 0) {
    Rcpp::stop("factorialC is not defined for negative integers.");
  }

  if (m <= 10) {
    return factorial_table[m];
  }

  double result = factorial_table[10]; // 10!
  for (int i = 11; i <= m; i++) {
    result *= i;
  }
  return result;
}


// Taylor expansion of the logarithm

// [[Rcpp::export]]
NumericVector tlogCPP(NumericVector x,
                      NumericVector a = NumericVector::create(1.0),
                      int k = 4, int d = 0) {
  int n_x = x.size();
  int len_a = a.size();
  if (!(len_a == 1 || len_a == n_x)) {
    stop("'a' must be either length 1 or the same length as 'x'.");
  }
  if (a.size() == 1) {
    a = NumericVector(n_x, a[0]);
  }

  if (k < 0) {
    stop("The polynomial order 'k' must be a non-negative integer scalar.");
  }
  if (d < 0) {
    stop("The derivative order 'd' must be a non-negative integer scalar.");
  }
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
          double gamma_n = factorialC(n - 1);     // gamma(n)
          double denom   = std::pow(-a[i], (double)d);
          accum += -gamma_n / denom;
        }
      } else {
        double mult  = (d == 0) ? 1.0 : factorialC(n) / factorialC(n - d);
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
