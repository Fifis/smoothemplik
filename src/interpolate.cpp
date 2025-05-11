#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cmath>

using namespace Rcpp;


// Fritsch--Carlson tangents
static std::vector<double> fritschCarlsonTangents(const std::vector<double>& x,
                                                  const std::vector<double>& y) {
  const std::size_t n = x.size();
  std::vector<double> d(n - 1);          // secant slopes
  std::vector<double> m(n);              // tangents

  for (std::size_t i = 0; i < n - 1; ++i)
    d[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);

  m[0]       = d[0];
  m[n - 1]   = d[n - 2];

  for (std::size_t i = 1; i < n - 1; ++i) {
    if (d[i - 1] * d[i] <= 0.0) {
      m[i] = 0.0;  // local extremum -- flatten
    } else {
      const double h0 = x[i]     - x[i - 1];
      const double h1 = x[i + 1] - x[i];
      m[i] = (3.0 * (h0 + h1)) /
        ((2.0 * h1 + h0) / d[i - 1] + (h1 + 2.0 * h0) / d[i]);
    }
  }
  return m;
}

// Evaluate monotone cubic at one point
static double monoHermiteEval(double xq, const std::vector<double>& x,
                              const std::vector<double>& y, const std::vector<double>& m) {
    const std::size_t n = x.size();
    if (xq <= x[0])          return y[0];
    if (xq >= x[n - 1])      return y[n - 1];

    // locate interval  [x[k], x[k+1]]
    std::size_t k = std::upper_bound(x.begin(), x.end(), xq) - x.begin() - 1;

    const double h = x[k + 1] - x[k];
    const double t = (xq - x[k]) / h;

    const double h00 = (2.0 * t * t * t - 3.0 * t * t + 1.0);
    const double h10 = (t * t * t - 2.0 * t * t + t);
    const double h01 = (-2.0 * t * t * t + 3.0 * t * t);
    const double h11 = (t * t * t - t * t);

    return  h00 * y[k]  +  h10 * h * m[k]  +  h01 * y[k + 1]  +  h11 * h * m[k + 1];
  }

// Quadratic reference parabola
static inline NumericVector ref_parabola(const NumericVector& x, double mu, double var) {
  NumericVector out = clone(x);
  for (double& v : out) v = std::pow(v - mu, 2.0) / var;
  return out;
}


// Concatenate two NumericVectors ---------------------------------------------
static NumericVector mergevecs(const NumericVector& a, const NumericVector& b) {
  NumericVector res(a.size() + b.size());
  std::copy(a.begin(), a.end(), res.begin());
  std::copy(b.begin(), b.end(), res.begin() + a.size());
  return res;
}


// [[Rcpp::export]]
NumericVector interpToHigherCPP(NumericVector x, Function f, double mean, double var,
                                double at, double gap) {
  const bool incr = at > mean;  // Is the function increasing
  const std::size_t n = x.size();
  NumericVector fx = f(x);
  NumericVector out = clone(fx);
  NumericVector fax = ref_parabola(x, mean, var); // reference parabola

  NumericVector xleft(3), xright(3);
  if (incr) {
    xleft  = NumericVector::create(at, at + 0.05 * gap, at + 0.10 * gap);
    xright = NumericVector::create(at + 0.90 * gap, at + 0.95 * gap, at + 1.00 * gap);

    if (is_true(all(x >= max(xright)))) return fax;  // far right = pure parabola

    NumericVector fleft  = f(xleft);
    NumericVector fright = ref_parabola(xright, mean, var);
    NumericVector xp = mergevecs(xleft, xright);
    NumericVector yp = mergevecs(fleft, fright);

    // Build monotone spline coefficients
    std::vector<double> xs(xp.begin(), xp.end());
    std::vector<double> ys(yp.begin(), yp.end());
    std::vector<double> ms = fritschCarlsonTangents(xs, ys);

    for (std::size_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (xi <  at) out[i] = fx[i];
      else if (xi <  at + gap)
        out[i] = monoHermiteEval(xi, xs, ys, ms);
      else out[i] = fax[i];
    }
  } else {                                   //  at <= mean  (decreasing arm)
    xleft  = NumericVector::create(at - 0.90 * gap, at - 0.95 * gap, at - 1.00 * gap);
    xright = NumericVector::create(at - 0.10 * gap, at - 0.05 * gap, at);

    if (is_true(all(x <= min(xleft)))) return fax;  // far left = pure parabola

    NumericVector fleft  = ref_parabola(xleft, mean, var);
    NumericVector fright = f(xright);
    NumericVector xp = mergevecs(xleft, xright);
    NumericVector yp = mergevecs(fleft, fright);
    // sort xp/yp so that xp is strictly increasing
    const std::size_t nk = xp.size();  // Number of knots
    IntegerVector ord = seq(0, nk - 1);
    std::sort(ord.begin(), ord.end(), [&](int i, int j){ return xp[i] < xp[j]; });
    NumericVector xp_sorted(nk), yp_sorted(nk);
    for (std::size_t k = 0; k < nk; ++k) {
      xp_sorted[k] = xp[ord[k]];
      yp_sorted[k] = yp[ord[k]];
    }
    std::vector<double> xs(xp_sorted.begin(), xp_sorted.end());
    std::vector<double> ys(yp_sorted.begin(), yp_sorted.end());
    std::vector<double> ms = fritschCarlsonTangents(xs, ys);

    for (std::size_t i = 0; i < n; ++i) {
      const double xi = x[i];
      if (xi > at) out[i] = fx[i];
      else if (xi > at - gap)
        out[i] = monoHermiteEval(xi, xs, ys, ms);
      else out[i] = fax[i];
    }
  }
  return out;
}


// [[Rcpp::export]]
NumericVector interpToLowerCPP(NumericVector x, Function f,
                 double mean, double var, double at, double gap)
{
  const bool incr = at > mean;
  const std::size_t nx = x.size();
  NumericVector fx = f(x);
  NumericVector out = clone(fx);
  NumericVector fax = ref_parabola(x, mean, var);

  NumericVector xleft(3), xright(3);
  double at2;

  if (incr) {                                // right-hand branch
    xleft  = NumericVector::create(at, at + 0.01 * gap, at + 0.02 * gap);
    at2    = std::sqrt(as<double>(f(at)) * var) + mean;
    xright = NumericVector::create(at2 + 0.98 * gap, at2 + 0.99 * gap, at2 + 1.00 * gap);

    if (is_true(all(x >= max(xright)))) return fax;

    NumericVector fleft  = f(xleft);
    NumericVector fright = ref_parabola(xright, mean, var);
    NumericVector xp = mergevecs(xleft, xright);
    NumericVector yp = mergevecs(fleft, fright);
    std::vector<double> xs(xp.begin(), xp.end());
    std::vector<double> ys(yp.begin(), yp.end());
    std::vector<double> ms = fritschCarlsonTangents(xs, ys);

    for (std::size_t j = 0; j < nx; ++j) {
      const double xi = x[j];
      if (xi <  at) out[j] = fx[j];
      else if (xi <  at2 + gap)
        out[j] = monoHermiteEval(xi, xs, ys, ms);
      else out[j] = fax[j];
    }
  } else {  // Left branch
    xright = NumericVector::create(at - 0.02 * gap, at - 0.01 * gap, at);
    at2    = -std::sqrt(as<double>(f(at)) * var) + mean;
    xleft  = NumericVector::create(at2 - 0.98 * gap, at2 - 0.99 * gap, at2 - 1.00 * gap);

    if (is_true(all(x <= min(xleft))))
      return fax;

    NumericVector fleft  = ref_parabola(xleft, mean, var);
    NumericVector fright = f(xright);

    NumericVector xp = mergevecs(xleft, xright);
    NumericVector yp = mergevecs(fleft, fright);

    // sort xp/yp
    const std::size_t nk = xp.size();  // Number of knots
    IntegerVector ord = seq(0, nk - 1);
    std::sort(ord.begin(), ord.end(), [&](int i, int j){ return xp[i] < xp[j]; });
    NumericVector xp_sorted(nk), yp_sorted(nk);
    for (std::size_t k = 0; k < nk; ++k) {
      xp_sorted[k] = xp[ord[k]];
      yp_sorted[k] = yp[ord[k]];
    }
    std::vector<double> xs(xp_sorted.begin(), xp_sorted.end());
    std::vector<double> ys(yp_sorted.begin(), yp_sorted.end());
    std::vector<double> ms = fritschCarlsonTangents(xs, ys);

    for (std::size_t j = 0; j < nx; ++j) {
      const double xi = x[j];
      if (xi >  at) out[j] = fx[j];
      else if (xi >  at2 - gap)
        out[j] = monoHermiteEval(xi, xs, ys, ms);
      else out[j] = fax[j];
    }
  }
  return out;
}

