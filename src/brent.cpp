#include <Rcpp.h>
#include <cmath>
#include <cfloat>  // for DBL_EPSILON

using namespace Rcpp;

/*
 This is an adaptation of Brent's local minimization from John Burkardt's code
 (similar to "local_min" or "R_zeroin2"-style logic), but with the following additions:
 - We track the number of iterations.
 - We stop when the standard Brent criterion is met, or if we exceed a maximum iteration count.
 - We return a list with elements:
 $root       (location of the min)
 $f.root     (function value at that location)
 $iter          (total iteration count used)
 $estim.prec    (our estimate of the final bracket size)
 There are no preliminary iterations.
 
 The code stores the approximate final bracket width in $estim.prec, like in "uniroot".
 If the minimiser is pinned to an end point, we set $estim.prec = NA.
 */

static const double EPS = DBL_EPSILON;  // Machine precision

// [[Rcpp::export]]
Rcpp::List brentMin(
    Rcpp::Function f,  // The objective function to minimize
    double lower,          // Left endpoint of interval
    double upper,          // Right endpoint of interval
    double tol = 1e-8, // Tolerance (bracket width)
    int maxiter = 200   // Maximum allowed iterations
) {
  if (lower >= upper) {
    stop("brentMin: 'lower' must be strictly less than 'upper'.");
  }
  if (tol <= 0.0 || !R_finite(tol)) {
    stop("brentMin: 'tol' must be > 0 and finite.");
  }

  // Convert the R function into a C++ lambda that returns a double
  std::function<double(double)> fn = [&](double xVal){
    Rcpp::RObject valR = f(xVal);
    if (valR.isNULL()) {
      warning("Function returned NULL, returning +Inf.");
      return DBL_MAX;
    }
    NumericVector nv(valR);
    if (nv.size() < 1) {
      warning("Function returned empty vector, returning +Inf.");
      return DBL_MAX;
    }
    double val = nv[0];
    if (!R_finite(val)) {
      warning("Function returned NA/Inf, returning +Inf.");
      return DBL_MAX;
    }
    return val;
  };

  // c = the golden ratio
  const double c = 0.5 * (3.0 - std::sqrt(5.0));

  double sa = lower;
  double sb = upper;

  // Initial guess
  double x  = sa + c * (sb - sa);
  double w  = x;
  double v  = w;
  double e  = 0.0;  // storage for the step size from the prior iteration

  double fx = fn(x);
  double fw = fx;
  double fv = fw;

  int iterCount = 0;
  double estimPrec = NA_REAL; // Final bracket estimate

  for (; iterCount < maxiter; iterCount++) {
    double m   = 0.5 * (sa + sb);
    double tol1= std::sqrt(EPS)*std::fabs(x) + tol;
    double t2  = 2.0 * tol1;

    // Stopping criterion: if the interval is sufficiently small
    if (std::fabs(x - m) <= (t2 - 0.5 * (sb - sa))) {
      // Convergence; store the approximate bracket size, like uniroot's fabs(c - upper)
      // If the solution is at an endpoint, set NA to match uniroot's style.
      if (x == sa || x == sb) {
        estimPrec = NA_REAL;
      } else {
        estimPrec = std::fabs(sb - sa);
      }
      break;
    }

    // Try a parabolic fit or fallback to golden section
    double r = 0.0, q = 0.0, p = 0.0;
    if (std::fabs(e) > tol1) {
      r = (x - w)*(fx - fv);
      q = (x - v)*(fx - fw);
      p = (x - v)*q - (x - w)*r;
      q = 2.0*(q - r);
      if (q > 0.0) p = -p; else q = -q;
      r = e;
      e = 0.0;

      // Check if the parabolic step is acceptable
      if ( (std::fabs(p) < std::fabs(0.5*q*r)) &&
        (p > q*(sa - x)) &&
        (p < q*(sb - x)) ) {
        e = p/q; // Parabolic step
      } else {
        if (x < m) e = sb - x; else e = sa - x; // Golden section
        e *= c;
      }
    } else {
      if (x < m) e = sb - x; else e = sa - x; // Golden section
      e *= c;
    }

    // If e is too small, push it to at least +/- tol1
    double d;
    if (std::fabs(e) >= tol1) {
      d = e;
    } else {
      d = (e > 0.0) ? tol1 : -tol1;
    }

    double u = x + d;
    double fu = fn(u);

    // Update the bracket [sa, sb], plus v, w, x
    if (fu <= fx) {
      if (u < x) {
        sb = x;
      } else {
        sa = x;
      }
      v  = w;   fv = fw;
      w  = x;   fw = fx;
      x  = u;   fx = fu;
    } else {
      if (u < x) {
        sa = u;
      } else {
        sb = u;
      }
      if ( (fu <= fw) || (w == x) ) {
        v  = w;   fv = fw;
        w  = u;   fw = fu;
      } else if ( (fu <= fv) || (v == x) || (v == w) ) {
        v  = u;   fv = fu;
      }
    }
  }

  // If we exit because of maxiter, we still set the bracket size:
  if (iterCount >= maxiter) {
    warning("brentMin: Maximum iteration limit reached.");
    // The final "x" is the best we have, so define a bracket-based precision:
    if (x == sa || x == sb) {
      estimPrec = NA_REAL;
    } else {
      estimPrec = std::fabs(sb - sa);
    }
  }

  return Rcpp::List::create(
    _["root"]     = x,
    _["f.root"]   = fx,
    _["iter"]        = iterCount,
	_["init.it"]        = NA_REAL,
    _["estim.prec"]  = estimPrec
  );
}


/*
 This is an adaptation of Brent's root search from John Burkardt's code
   1. Take an R function f + a bracket [a, upper].
   2. Assume that f(lower) and f(upper) have opposite signs (so there's a root).
   3. Return a list analogous to uniroot's output in base R:
      $root        the final solution
      $f.root      the function value at that solution
      $iter        the number of iterations used
      $estim.prec  the approximate final bracket width
                   (NA if the solution is at an endpoint)
*/

// [[Rcpp::export]]
Rcpp::List brentZero(
    Rcpp::Function f,  // R function with one argument -> returns double
    double lower,          // left endpoint (f(lower) <> 0)
    double upper,          // right endpoint (f(upper) >< 0)
    double tol = 1e-8, // absolute error tolerance
    int maxiter = 100   // maximum number of iterations
) {
  if (lower >= upper) {
    stop("brentZero: 'lower' must be strictly less than 'upper'.");
  }
  if (tol <= 0.0 || !R_finite(tol)) {
    stop("brentZero: 'tol' must be > 0 and finite.");
  }

  std::function<double(double)> fn = [&](double xVal) {
    Rcpp::RObject valR = f(xVal);
    if (valR.isNULL()) {
      Rcpp::warning("Function returned NULL, returning NA for safety.");
      return NA_REAL;
    }
    NumericVector nv(valR);
    if (nv.size() < 1) {
      Rcpp::warning("Function returned empty vector, returning NA for safety.");
      return NA_REAL;
    }
    double val = nv[0];
    if (!R_finite(val)) {
      Rcpp::warning("Function returned NA/Inf, returning NA for safety.");
      return NA_REAL;
    }
    return val;
  };

  const double macheps = DBL_EPSILON;

  double sa = lower;
  double sb = upper;
  double fa = fn(sa);
  double fb = fn(sb);

  // Check for sign difference
  if (fa * fb > 0.0) {
    Rcpp::stop("f(lower) and f(upper) must have opposite signs.");
  }

  double c  = sa;
  double fc = fa;
  double e  = sb - sa;
  double d  = e; // step size in iteration

  int iterCount = 0;
  double root   = sb; // To store the final solution
  double froot  = fb;

  if (fa == 0.0) {
    // Found lower root at 'lower'
    // iterations used = 0, tolerance = 0
    return List::create(
      _["root"]       = lower,
      _["f.root"]     = fa,
      _["iter"]       = 0,
      _["init.it"]    = NA_REAL,
      _["estim.prec"] = 0.0
    );
  }
  if (fb == 0.0) {
    return List::create(
      _["root"]       = upper,
      _["f.root"]     = fb,
      _["iter"]       = 0,
      _["init.it"]    = NA_REAL,
      _["estim.prec"] = 0.0
    );
  }

  for (; iterCount < maxiter; iterCount++) {

    // If fc is smaller in magnitude than fb, swap roles
    if (std::fabs(fc) < std::fabs(fb)) {
      sa = sb;
      sb = c;
      c  = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    double tol1 = 2.0 * macheps * std::fabs(sb) + tol;
    double m    = 0.5 * (c - sb);

    // Stopping criterion: if bracket is small or fb == 0
    if ( (std::fabs(m) <= tol1) || (fb == 0.0) ) {
      // Converged
      root  = sb;
      froot = fb;
      break;
    }

    if ( (std::fabs(e) < tol1) || (std::fabs(fa) <= std::fabs(fb)) ) {
      e = m;
      d = e;
    } else {
      double s = fb / fa;
      double p, q, r;
      if (sa == c) {
        // linear interpolation
        p = 2.0 * m * s;
        q = 1.0 - s;
      } else {
        // inverse quadratic
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }

      if (p > 0.0) {
        q = -q;
      } else {
        p = -p;
      }

      double sOld = e;
      e = d;

      // Accept interpolation if it falls well within the bracket
      if ( (2.0 * p < 3.0 * m * q - std::fabs(tol1 * q)) &&
           (p < std::fabs(0.5 * sOld * q)) ) {
        d = p / q;
      } else {
        // otherwise do a bisection
        e = m;
        d = e;
      }
    }

    sa = sb;
    fa = fb;

    if (std::fabs(d) > tol1) {
      sb = sb + d;
    } else if (m > 0.0) {
      sb = sb + tol1;
    } else {
      sb = sb - tol1;
    }

    fb = fn(sb);

    // If the sign of fb changed relative to fc, update
    if ((fb > 0.0 && fc > 0.0) || (fb <= 0.0 && fc <= 0.0)) {
      c  = sa;
      fc = fa;
      e  = sb - sa;
      d  = e;
    }
  }

  // If we exited because max iterations were used, we still have our best guess:
  if (iterCount >= maxiter) {
    Rcpp::warning("brentZero: maximum iteration limit reached.");
    // last computed 'sb' is the best we have
    root  = sb;
    froot = fb;
  }

  // Estimate precision from the last bracket: typically |c - sb|
  // or if root is pinned at an endpoint, set NA similarly to uniroot
  double estimPrec = std::fabs(c - sb);
  if (root == lower || root == upper) {
    estimPrec = NA_REAL;
  }

  // Build the result
  return Rcpp::List::create(
    _["root"]       = root,
    _["f.root"]     = froot,
    _["iter"]       = iterCount,
    _["init.it"]    = NA_REAL,
    _["estim.prec"] = estimPrec
  );
}
