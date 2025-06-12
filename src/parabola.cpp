#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector getParabola3CPP(const NumericVector& x, const NumericVector& y) {
  if (x.size() != 3 || y.size() != 3)
    Rcpp::stop("Input vectors x and y must each be of length 3.");

  const double me = std::numeric_limits<double>::epsilon();
  double x0 = x[0], x1 = x[1], x2 = x[2];
  double y0 = y[0], y1 = y[1], y2 = y[2];

  if (!std::isfinite(x0) || !std::isfinite(x1) || !std::isfinite(x2) ||
      !std::isfinite(y0) || !std::isfinite(y1) || !std::isfinite(y2)) {
      stop("Input values must be finite numbers (no NA/NaN/Inf).");
  }

  if (x0 == x1 || x0 == x2 || x1 == x2)
    stop("Input x values must be distinct.");

  double A1 = x0*x0 - x1*x1;
  double B1 = x0 - x1;
  double C1 = y0 - y1;
  double A2 = x1*x1 - x2*x2;
  double B2 = x1 - x2;
  double C2 = y1 - y2;
  // Compute the 2x2 determinant for the system solving a and b
  double det = A1 * B2 - A2 * B1;
  if (std::fabs(det) < me)
    warning("Poor quality of fitted parabola (the determinant in the denominator is less than the machine epsilon).");

  double a = (C1 * B2 - C2 * B1) / det;
  double b = (A1 * C2 - A2 * C1) / det;
  double c = y0 - a * x0 * x0 - b * x0;
  return NumericVector::create(a, b, c);
}


// [[Rcpp::export]]
NumericVector getParabolaCPP(double x, double f, double fp, double fpp) {
  if (!std::isfinite(x) || !std::isfinite(f) ||
      !std::isfinite(fp) || !std::isfinite(fpp))
      Rcpp::stop("Input values must be finite numbers (no NA/NaN/Inf).");
  // Compute coefficients for parabola matching f, f', f'' at x
  // a = fpp/2, b = fp - fpp*x, c = f - fp*x + 0.5*fpp*x^2
  double a = 0.5 * fpp;
  double b = fp - fpp * x;
  double c = f - fp * x + 0.5 * fpp * x * x;
  return NumericVector::create(a, b, c);
}
