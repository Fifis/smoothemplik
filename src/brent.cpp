#include <Rcpp.h>
#include <cmath>
#include <cfloat>  // for DBL_EPSILON
#include <string>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List brentMinCPP(
    Function f,
    NumericVector interval = NumericVector(),
    double lower = NA_REAL,
    double upper = NA_REAL,
    double tol = 1e-8,
    int maxiter = 500,
    int trace = 0
) {

  // If the user gave an interval, replicate "if (!missing(interval) && length(interval)==2)"
  if (interval.size() == 2) {
    lower = std::min(interval[0], interval[1]);
    upper = std::max(interval[0], interval[1]);
  } else if (interval.size() != 0 && interval.size() != 2) {
    stop("brentMin: 'interval' must be a vector of length 2 if given, else empty.");
  }
  if (!R_finite(lower) || !R_finite(upper) || lower >= upper) {
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

  if (trace >= 1) {
    Rcpp::Rcout << "brentMin: Starting search in ["
    << sa << ", " << sb << "], initial x=" << x
    << ", fx=" << fx << "\n";
  }

  for (; iterCount < maxiter; iterCount++) {
    double m   = 0.5 * (sa + sb);
    double tol1= std::sqrt(DBL_EPSILON) * std::fabs(x) + tol;
    double t2  = 2.0 * tol1;

    if (trace >= 2) {
      Rcpp::Rcout << "Iter " << iterCount
      << ": bracket=[" << sa << ", " << sb << "], x=" << x
      << ", fx=" << fx << "\n";
    }

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

  if (trace >= 1) {
    Rcpp::Rcout << "brentMin: Done in " << iterCount << " iterations, x="
    << x << ", fx=" << fx << ", bracket=[" << sa << ", " << sb
    << "]\n";
  }

  return Rcpp::List::create(
    _["root"] = x, _["f.root"]   = fx,
    _["iter"] = iterCount, _["init.it"] = NA_REAL, _["estim.prec"]  = estimPrec
  );
}


// Helper to compute sign(x); sign(0) = 0
static inline double getsign(double x) {
  return (x > 0.0) ? 1.0 : (x < 0.0 ? -1.0 : 0.0);
}

// The "extendInt" logic from R's uniroot.
// match.arg(extendInt) in R is one of {"no","yes","downX","upX"}.
enum ExtendIntMode {
  EXT_NO = 0,
  EXT_YES,
  EXT_DOWN,
  EXT_UP,
  EXT_RIGHT,
  EXT_LEFT
};

static ExtendIntMode parseExtendInt(const std::string &str) {
  if (str == "no")    return EXT_NO;
  if (str == "yes")   return EXT_YES;
  if (str == "downX") return EXT_DOWN;
  if (str == "upX")   return EXT_UP;
  if (str == "right") return EXT_RIGHT;
  if (str == "left")  return EXT_LEFT;
  if (str == "upper") return EXT_RIGHT;  // Safe catching of popular mis-spellings
  if (str == "lower") return EXT_LEFT;
  Rcpp::stop("Invalid 'extendInt' argument; must be one of 'no', 'yes', 'downX', 'upX', 'right', 'left'.");
  return EXT_NO; // not reached
}

// Clamp to machine max
static inline double truncateBig(double x) {
  // like R's pmax.int(pmin(x, .Machine$double.xmax), -.Machine$double.xmax)
  const double big = DBL_MAX;
  if (x >  big) return  big;
  if (x < -big) return -big;
  return x;
}

// Delta(u) = 0.01 * pmax(1e-4, abs(u))
static inline double Delta(double u) {
  return 0.01 * std::max(1e-4, std::fabs(u));
}

// ----------------------------------------------------------------------
// This function extends [lower,upper] until f(lower)*f(upper) <= 0
// or until we hit the maximum allowed extension steps.
// - If mode=YES, we expand in both directions as needed (like uniroot "yes").
// - If mode=DOWN, we only expand downward (like "downX") and then possibly upward if needed.
// - If mode=UP,   we only expand upward (like "upX") and then possibly downward if needed.
// - If mode=NO,   we do nothing.
//
// This extension step count is returned in 'initSteps', and the final
// bracket is returned in [lo,hi], with f(lo), f(hi) returned in fLo, fHi.
//
// Instead of throwing when we exceed maxSteps or fail to find a sign change, we set softFail=true
// and return the current "best" endpoint (smaller |f|) in candX/candFX.

void doExtendInterval(
  std::function<double(double)> fn, // the function
                      double &lo, double &hi,           // in/out bracket endpoints
                      double &fLo, double &fHi,         // in/out f-values at endpoints
                      int &initSteps,                   // out: how many expansions we did
                      int maxSteps,                     // maximum allowed expansions
                      ExtendIntMode mode,               // which mode
                      int trace,                        // printing level
                      bool &softFail,                   // out: did we hit a soft failure?
                      double &candX, double &candFX     // out: best candidate so far
) {
  initSteps = 0;
  softFail  = false;
  candX     = NA_REAL;
  candFX    = NA_REAL;

  // truncate the endpoints to avoid overflows as R does:
  double fLoTrunc = truncateBig(fLo);
  double fHiTrunc = truncateBig(fHi);

  // doX logic:
  // If mode = yes (EXT_YES), do expansions if fLoTrunc*fHiTrunc > 0
  // If mode = no (EXT_NO), skip expansions
  // If mode = downX => sign(...) * fLoTrunc>0 or sign(...) * fHiTrunc<0
  // If mode = upX => similarly
  // If mode = left/right => extend only that side if no sign change yet
  bool doX = false;
  if (mode == EXT_YES) {
    if (fLoTrunc * fHiTrunc > 0) doX = true;
  } else if (mode == EXT_NO) {
    doX = false;
  } else if (mode == EXT_LEFT || mode == EXT_RIGHT) {
    if (fLoTrunc * fHiTrunc > 0) doX = true; // Only if we still need a sign change
  } else {
    // mode is EXT_DOWN or EXT_UP => Sig = -1 or +1
    double Sig = (mode == EXT_DOWN ? -1.0 : 1.0);
    if (Sig * fLoTrunc > 0 || Sig * fHiTrunc < 0) doX = true;
  }

  if (!doX) {
    return; // no extension needed
  }

  if (trace >= 1) {
    Rcpp::Rcout << "Extending interval from [" << lo << "," << hi << "] ...\n";
  }

  auto setSoftFailWithBest = [&](){
    // choose endpoint with smaller |f|
    if (std::fabs(fLo) <= std::fabs(fHi)) { candX = lo; candFX = fLo; }
    else { candX = hi; candFX = fHi; }
    softFail = true;
  };

  if (mode == EXT_YES) {
    // Expand both ends if fLo*fHi > 0
    // This matches uniroot's 'while(f.lower*f.upper>0)'
    // doubling the step each time
    std::vector<double> d(2);
    d[0] = Delta(lo);
    d[1] = Delta(hi);

    while ((fLo * fHi > 0.0) && (std::isfinite(lo) || std::isfinite(hi))) {
      initSteps++;
      if (initSteps > maxSteps) {
        Rcpp::warning("Extension step limit hit after %d steps (soft stop).\n", initSteps-1);
        setSoftFailWithBest();
        return;
      }
      // Try to move lower downward
      if (std::isfinite(lo)) {
        double oldLo = lo;
        double oldf  = fLo;
        lo = lo - d[0];
        double tmp = fn(lo);
        if (Rcpp::NumericVector::is_na(tmp)) {
          // revert
          lo   = oldLo;
          fLo  = oldf;
          d[0] = d[0] / 4.0;
        } else {
          fLo = tmp;
        }
      }

      // Try to move upper upward
      if (std::isfinite(hi)) {
        double oldHi = hi;
        double oldf  = fHi;
        hi = hi + d[1];
        double tmp = fn(hi);
        if (Rcpp::NumericVector::is_na(tmp)) {
          // revert
          hi   = oldHi;
          fHi  = oldf;
          d[1] = d[1] / 4.0;
        } else {
          fHi = tmp;
        }
      }
      if (trace >= 2) {
        Rcpp::Rcout << " .. step " << initSteps
        << " => bracket=[" << lo << "," << hi << "], f(lo)="
        << fLo << ", f(hi)=" << fHi << "\n";
      }
      // double the step
      d[0] *= 2.0;
      d[1] *= 2.0;
    }
  } else if (mode == EXT_DOWN || mode == EXT_UP) {
    // Sig is -1 or +1
    double Sig = (mode == EXT_DOWN ? -1.0 : 1.0);

    // Expand lower if Sig*fLo>0
    double dd = Delta(lo);
    while ((Sig * fLo > 0.0)) {
      initSteps++;
      if (initSteps > maxSteps) {
        Rcpp::warning("Lower extension step limit hit after %d steps (soft stop).\n", initSteps-1);
        // We were extending the lower endpoint
        candX  = lo;
        candFX = fLo;
        softFail = true;
        return;
      }
      // Move lower endpoint
      if (mode == EXT_DOWN) {
        lo -= std::fabs(dd);  // ensure we go downward
      } else {
        // upX => lo stays put, or we can do the symmetrical logic
        // But uniroot does "Sig * f.lower > 0 => we shift the lower down if Sig=-1
        // or shift it up if Sig=+1. It's a bit unusual, but let's keep it.
        lo -= Sig*dd;
      }

      fLo = fn(lo);
      if (trace >= 2) {
        Rcpp::Rcout << " .. extended lower => lo=" << lo << ", f(lo)=" << fLo
        << ", step=" << dd << "\n";
      }
      dd *= 2.0;
    }

    // Expand upper if Sig*fHi<0
    dd = Delta(hi);
    while ((Sig * fHi < 0.0)) {
      initSteps++;
      if (initSteps > maxSteps) {
        Rcpp::warning("Upper extension step limit hit after %d steps (soft stop).\n", initSteps-1);
        // We were extending the upper endpoint
        candX  = hi;
        candFX = fHi;
        softFail = true;
        return;
      }
      if (mode == EXT_UP) {
        hi += std::fabs(dd); // push hi upward
      } else {
        hi -= Sig*dd;
      }

      fHi = fn(hi);
      if (trace >= 2) {
        Rcpp::Rcout << " .. extended upper => hi=" << hi << ", f(hi)=" << fHi
        << ", step=" << dd << "\n";
      }
      dd *= 2.0;
    }
  } else if (mode == EXT_LEFT) {
    // Extend lower (left) end only until sign change or we hit maxSteps
    double dd = Delta(lo);
    bool   seen_bad = false;
    double bad = NA_REAL;             // first non-finite to the left
    double good = lo;                 // last finite to the left
    double fgood = fLo;

    while ((getsign(fLo) * getsign(fHi) > 0.0)) {
      initSteps++;
      if (initSteps > maxSteps) {
        Rcpp::warning("Left-only extension step limit hit after %d steps (soft stop).\n", initSteps-1);
        if (std::fabs(fLo) <= std::fabs(fHi)) { candX = lo; candFX = fLo; }
        else { candX = hi; candFX = fHi; }
        softFail = true;
        return;
      }

      double trial = lo - dd;
      double tmp = fn(trial);

      if (!R_finite(tmp)) {
        // Hit a non-finite barrier on the left; remember it and backtrack once.
        if (!seen_bad) { seen_bad = true; bad = trial; }
        if (trace >= 2) Rcpp::Rcout << " .. left-only backtrack from non-finite\n";

        // Binary search in (trial,best_finite) to get the closest finite "good" below the barrier
        double right = good,  f_right = fgood;  // finite
        double left  = trial;                   // non-finite
        int    tries = 0;
        while (tries++ < maxSteps) {
          double mid  = 0.5*(left + right);
          double fmid = fn(mid);
          if (!R_finite(fmid)) {
            left = mid;           // still non-finite
          } else {
            right = mid; f_right = fmid;  // finite, move toward barrier
          }
          if (std::fabs(right - left) <= 2.0*DBL_EPSILON*std::fabs(right)) break;
        }
        lo = right; fLo = f_right; good = lo; fgood = fLo;

        // From now on, stay inside (good,bad): shrink step to the gap
        dd = 0.5 * (std::fabs(bad - good));
        if (!(dd > 0.0)) dd = Delta(lo);

      } else {
        // finite step
        lo = trial; fLo = tmp; good = lo; fgood = fLo;
        if (!seen_bad) dd *= 2.0;     // still exploring
        else {
          dd = 0.5 * (std::fabs(bad - good));   // approach barrier from inside
          if (!(dd > 0.0)) dd = Delta(lo);
        }
      }

      if (trace >= 2) {
        Rcpp::Rcout << " .. left-only step " << initSteps
                    << " => lo=" << lo << ", f(lo)=" << fLo << "\n";
      }

      if (getsign(fLo) * getsign(fHi) <= 0.0) break; // got a sign change

      // If we know the barrier and the window collapsed, stop with soft fail
      if (seen_bad && std::fabs(bad - good) <= 2.0*DBL_EPSILON*std::fabs(good)) {
        candX  = (std::fabs(fLo) <= std::fabs(fHi)) ? lo : hi;
        candFX = (candX == lo) ? fLo : fHi;
        softFail = true; return;
      }
    }

  } else if (mode == EXT_RIGHT) {
    // Extend upper (right) end only until sign change or we hit maxSteps
    double dd = Delta(hi);
    bool   seen_bad = false;
    double bad = NA_REAL;             // first non-finite to the right
    double good = hi;                 // last finite to the right
    double fgood = fHi;

    while ((getsign(fLo) * getsign(fHi) > 0.0)) {
      initSteps++;
      if (initSteps > maxSteps) {
        Rcpp::warning("Right-only extension step limit hit after %d steps (soft stop).\n", initSteps-1);
        if (std::fabs(fLo) <= std::fabs(fHi)) { candX = lo; candFX = fLo; }
        else { candX = hi; candFX = fHi; }
        softFail = true;
        return;
      }

      double trial = hi + dd;
      double tmp   = fn(trial);

      if (!R_finite(tmp)) {
        // Hit a non-finite barrier on the right; remember it and backtrack once.
        if (!seen_bad) { seen_bad = true; bad = trial; }
        if (trace >= 2) Rcpp::Rcout << " .. right-only backtrack from non-finite\n";

        // Binary search in (best_finite, trial) to get the closest finite "good" below the barrier
        double left  = good,  f_left = fgood;   // finite
        double right = trial;                   // non-finite
        int    tries = 0;
        while (tries++ < maxSteps) {
          double mid  = 0.5*(left + right);
          double fmid = fn(mid);
          if (!R_finite(fmid)) {
            right = mid;          // still non-finite
          } else {
            left = mid; f_left = fmid;  // finite, move toward barrier
          }
          if (std::fabs(right - left) <= 2.0*DBL_EPSILON*std::fabs(left)) break;
        }
        hi = left; fHi = f_left; good = hi; fgood = fHi;

        // From now on, stay inside (good,bad): shrink step to the gap
        dd = 0.5 * (std::fabs(bad - good));
        if (!(dd > 0.0)) dd = Delta(hi);
      } else {
        // finite step
        hi = trial; fHi = tmp; good = hi; fgood = fHi;
        if (!seen_bad) dd *= 2.0;     // still exploring
        else {
          dd = 0.5 * (std::fabs(bad - good));   // approach barrier from inside
          if (!(dd > 0.0)) dd = Delta(hi);
        }
      }

      if (trace >= 2) {
        Rcpp::Rcout << " .. right-only step " << initSteps
                    << " => hi=" << hi << ", f(hi)=" << fHi << "\n";
      }

      if (getsign(fLo) * getsign(fHi) <= 0.0) break; // got a sign change

      // If we know the barrier and the window collapsed, stop with soft fail
      if (seen_bad && std::fabs(bad - good) <= 2.0*DBL_EPSILON*std::fabs(good)) {
        candX  = (std::fabs(fLo) <= std::fabs(fHi)) ? lo : hi;
        candFX = (candX == lo) ? fLo : fHi;
        softFail = true; return;
      }
    }
  }

  // Possibly a short message
  if (trace > 0 && initSteps > 0) {
    Rcpp::Rcout << "Extended to [" << lo << "," << hi << "] in " << initSteps << " steps\n";
  }

  // If still no sign change, soft stop with best endpoint
  if (getsign(fLo) * getsign(fHi) > 0.0) {
    Rcpp::warning("No sign change found during extension in 'yes' mode after %d steps.", initSteps-1);
    setSoftFailWithBest();
    return;
  }
}


// This is an adaptation of Brent's root search from John Burkardt's code
//   1. Take an R function f + a bracket [a, upper].
//   2. Assume that f(lower) and f(upper) have opposite signs (so there's a root).
//   3. Return a list analogous to uniroot's output in base R:
//      $root        the final solution
//      $f.root      the function value at that solution
//      $iter        the number of iterations used
//      $init.it     the number of initial iterations to find the sign change
//      $estim.prec  the approximate final bracket width
//                   (NA if the solution is at an endpoint)
//      $exitcode    0 ok, 1 extension max steps, 2 main-loop maxiter

// [[Rcpp::export]]
List brentZeroCPP(
  Function f,
  NumericVector interval = NumericVector(),
  double lower = NA_REAL,
  double upper = NA_REAL,
  Nullable<NumericVector> f_lower = R_NilValue,
  Nullable<NumericVector> f_upper = R_NilValue,
  std::string extendInt = "no",
  double tol = 1e-8,
  int maxiter = 500,
  int trace = 0
) {

  // Count all function evaluations (nfeval) and use a single wrapper
  long long nfeval = 0;
  auto eval = [&](double x)->double {
    double y = as<double>(f(x));
    ++nfeval;
    return y;
  };

  // If the user gave an interval, replicate "if (!missing(interval) && length(interval)==2)"
  if (interval.size() == 2) {
    lower = std::min(interval[0], interval[1]);
    upper = std::max(interval[0], interval[1]);
  } else if (interval.size() != 0 && interval.size() != 2) {
    stop("brentZero: 'interval' must be a vector of length 2 if given, else empty.");
  }
  if (!R_finite(lower) || !R_finite(upper) || lower >= upper) {
    stop("brentZero: 'lower' must be strictly less than 'upper'.");
  }

  // Evaluate f at lower, upper if not supplied
  double fa, fb;
  if (f_lower.isNotNull()) {
    NumericVector tmp(f_lower);
    if (tmp.size() < 1) stop("f.lower invalid");
    fa = tmp[0];
  } else {
    fa = eval(lower);
  }
  if (f_upper.isNotNull()) {
    NumericVector tmp(f_upper);
    if (tmp.size() < 1) stop("f.upper invalid");
    fb = tmp[0];
  } else {
    fb = eval(upper);
  }

  if (ISNAN(fa)) stop("f.lower = f(lower) is NA");
  if (ISNAN(fb)) stop("f.upper = f(upper) is NA");

  ExtendIntMode mode = parseExtendInt(extendInt);

  // Preliminary extension pass if necessary
  int initSteps = 0;
  bool softFail = false;
  double candX = NA_REAL, candFX = NA_REAL;

  doExtendInterval(
    [&](double xx){ return eval(xx); },
    lower, upper, fa, fb, initSteps, maxiter, mode, trace, softFail, candX, candFX
  );

  if (softFail) {  // Soft stop: return the current best endpoint
    return List::create(
      _["root"]       = candX,
      _["f.root"]     = candFX,
      _["iter"]       = initSteps,  // only extension iters
      _["init.it"]    = (initSteps > 0 ? initSteps : R_NaInt),
      _["estim.prec"] = NA_REAL,
      _["exitcode"]   = 1
    );
  }

  // Now [lower, upper] should have opposite signs
  if (getsign(fa) * getsign(fb) > 0) {
    // Shouldn't happen due to soft handling; safeguard any way
    double root = (std::fabs(fa) <= std::fabs(fb)) ? lower : upper;
    double froot = (std::fabs(fa) <= std::fabs(fb)) ? fa : fb;
    return List::create(
      _["root"]       = root,
      _["f.root"]     = froot,
      _["iter"]       = initSteps,
      _["init.it"]    = (initSteps > 0 ? initSteps : R_NaInt),
      _["estim.prec"] = NA_REAL,
      _["exitcode"]   = 1
    );
  }

  // Ensure BOTH end-points are finite; if one end is not, shrink it inward
  auto shrink_right_to_finite = [&]() {
    int tries = 0;
    if (trace >= 1 && !R_finite(fb))
      Rcpp::Rcout << "Shrinking right endpoint to obtain finite f(upper) ...\n";
    while (!R_finite(fb) && tries < maxiter) {
      upper = 0.5 * (upper + lower);     // move inward
      fb = eval(upper)  ;
      tries++;
      if (trace >= 2)
        Rcpp::Rcout << " .. right => upper=" << upper << ", f(upper)=" << fb << "\n";
      double tol1 = 2.0*DBL_EPSILON*std::fabs(upper) + tol;
      if (std::fabs(upper - lower) <= tol1) break;
    }
  };

  auto shrink_left_to_finite = [&]() {
    int tries = 0;
    if (trace >= 1 && !R_finite(fa))
      Rcpp::Rcout << "Shrinking left endpoint to obtain finite f(lower) ...\n";
    while (!R_finite(fa) && tries < maxiter) {
      lower = 0.5 * (lower + upper);     // move inward
      fa = eval(lower);
      tries++;
      if (trace >= 2)
        Rcpp::Rcout << " .. left => lower=" << lower << ", f(lower)=" << fa << "\n";
      double tol1 = 2.0*DBL_EPSILON*std::fabs(lower) + tol;
      if (std::fabs(upper - lower) <= tol1) break;
    }
  };

  // Apply shrinkers when appropriate
  if (R_finite(fa) && !R_finite(fb)) {
    if (mode == EXT_UP || mode == EXT_RIGHT || mode == EXT_YES || mode == EXT_NO) shrink_right_to_finite();
  } else if (!R_finite(fa) && R_finite(fb)) {
    if (mode == EXT_DOWN || mode == EXT_LEFT || mode == EXT_YES || mode == EXT_NO) shrink_left_to_finite();
  }

  // After shrinking, if we STILL do not have opposite signs, try to (re)extend
  // on the allowed side to create a valid bracket before entering Brent.
  if (getsign(fa) * getsign(fb) > 0) {
    int addInit = 0; bool sf2 = false; double cx = NA_REAL, cfx = NA_REAL;
    double lo2 = lower, hi2 = upper, fLo2 = fa, fHi2 = fb;

    // Choose which way to re-extend
    ExtendIntMode reMode = mode;
    if (mode == EXT_RIGHT || mode == EXT_UP)      reMode = EXT_RIGHT;
    else if (mode == EXT_LEFT || mode == EXT_DOWN) reMode = EXT_LEFT;
    else if (mode == EXT_NO)                       reMode = EXT_YES;  // allow either side
    // EXT_YES stays EXT_YES

    doExtendInterval(
      [&](double xx){ return eval(xx); },
      lo2, hi2, fLo2, fHi2, addInit, maxiter, reMode, trace, sf2, cx, cfx
    );

    if (sf2) {
      // Could not recover a bracket: return best endpoint
      return List::create(
        _["root"]       = cx,
        _["f.root"]     = cfx,
        _["iter"]       = initSteps + addInit,
        _["init.it"]    = (initSteps + addInit > 0 ? initSteps + addInit : R_NaInt),
        _["estim.prec"] = NA_REAL,
        _["exitcode"]   = 1
      );
    }

    // Adopt the recovered bracket
    lower = lo2; upper = hi2; fa = fLo2; fb = fHi2;
    initSteps += addInit;
  }

  // Final safety: if still same sign, bail with best endpoint
  if (getsign(fa) * getsign(fb) > 0) {
    double root = (std::fabs(fa) <= std::fabs(fb)) ? lower : upper;
    double froot = (root == lower) ? fa : fb;
    return List::create(
      _["root"]       = root,
      _["f.root"]     = froot,
      _["iter"]       = initSteps,
      _["init.it"]    = (initSteps ? initSteps : R_NaInt),
      _["estim.prec"] = NA_REAL,
      _["exitcode"]   = 1
    );
  }


  //////////////////////////////////////////////////////////////////////////////
  // Standard Brent iteration. 'iter' is tracked in local code.
  double sa = lower;
  double sb = upper;
  double fsa = fa;
  double fsb = fb;

  // If f=0 at an endpoint, done
  if (fsa == 0.0) {
    int finalIter = 0 + initSteps;
    return List::create(
      _["root"]       = sa, _["f.root"] = fsa,
      _["iter"]       = finalIter, _["init.it"] = initSteps,
      _["estim.prec"] = 0.0, _["exitcode"]   = 0
    );
  }
  if (fsb == 0.0) {
    int finalIter = 0 + initSteps;
    return List::create(
      _["root"]       = sb, _["f.root"] = fsb,
      _["iter"]       = finalIter, _["init.it"] = initSteps,
      _["estim.prec"] = 0.0, _["exitcode"]   = 0
    );
  }

  // Standard brent loop
  double c  = sa;
  double fc = fsa;
  double e  = sb - sa;
  double d  = e;
  int    iterCount = 0;
  double root   = sb;
  double froot  = fsb;

  if (trace >= 1) {
    Rcpp::Rcout << "brentZero: Searching for the root in [" << sa << "," << sb
    << "], f(sa)=" << fsa << ", f(sb)=" << fsb
    << ", init.it=" << initSteps << "\n";
  }

  for (; iterCount < maxiter; iterCount++) {
    if (trace >= 2) {
      double bl = std::min(sb, c), br = std::max(sb, c);
      Rcpp::Rcout << "   Iter " << iterCount
                  << ": bracket=[" << bl << ", " << br << "], fb=" << fsb
                  << ", fc=" << fc << "\n";
    }
    // If fc is smaller in magnitude than fb, swap roles
    if (std::fabs(fc) < std::fabs(fsb)) {
      std::swap(sb, c);
      std::swap(fsb, fc);
      // Do NOT swap sa and fsa here.
    }
    double tol1 = 2.0*DBL_EPSILON*std::fabs(sb) + tol;
    double m    = 0.5*(c - sb);
    if ((std::fabs(m) <= tol1) || (fsb == 0.0)) {
      root  = sb;
      froot = fsb;
      break;
    }
    if ((std::fabs(e) < tol1) || (std::fabs(fsa) <= std::fabs(fsb))) {
      e = m;
      d = e;
    } else {
      double s = fsb / fsa;
      double p, q, r;
      if (sa == c) {
        // linear interpolation -- secant
        p = 2.0*m*s;
        q = 1.0 - s;
      } else {
        // inverse quadratic interpolation
        q = fsa/fc;
        r = fsb/fc;
        p = s*(2.0*m*q*(q-r) - (sb - sa)*(r - 1.0));
        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
      }
      if (p > 0.0) q = -q; else p = -p;
      double sOld = e;
      e = d;
      // Accept interpolation if it falls well within the bracket
      if ((2.0*p < 3.0*m*q - std::fabs(tol1*q)) && (p < std::fabs(0.5*sOld*q))) {
        d = p / q;
      } else {
        e = m;
        d = e;
      }
    }
    sa  = sb;
    fsa = fsb;

    if (std::fabs(d) > tol1) {
      sb += d;
    } else if (m > 0.0) {
      sb += tol1;
    } else {
      sb -= tol1;
    }
    fsb = eval(sb);
    // if (!R_finite(fsb)) fsb=0.0; // or handle NA as 0, etc.
    // Reject non-finite trial points; bisect inside the bracket
    if (!R_finite(fsb)) {
      // Fall back to bisection toward the interior of [sb,c]
      // (sb and c carry opposite signs by construction)
      do {
        sb = 0.5 * (sb + c);
        fsb = eval(fsb);
        tol1 = 2.0*DBL_EPSILON*std::fabs(sb) + tol;
      } while (!R_finite(fsb) && std::fabs(c - sb) > tol1);

      // If even repeated bisection cannot find a finite value, bail with best finite endpoint.
      if (!R_finite(fsb)) {
        double rootCand = (std::isfinite(fsa) && (std::fabs(fsa) <= std::fabs(fc))) ? sa : c;
        double fCand    = (rootCand == sa) ? fsa : fc;
        return List::create(
          _["root"] = rootCand, _["f.root"] = fCand,
          _["iter"] = iterCount + initSteps,
          _["init.it"] = (initSteps ? initSteps : R_NaInt),
          _["estim.prec"] = std::fabs(c - sb),
          _["exitcode"] = 2
        );
      }
    }

    // Maintain the bracket: ensure sign(fsb) != sign(fc)
    if (getsign(fsb) == getsign(fc)) {
      c  = sa;
      fc = fsa;
      e  = sb - sa;
      d  = e;
    }
  }

  // If we exit due to iteration limit
  int  exitcode = 0;
  if (iterCount >= maxiter) {
    exitcode = 2; // main-loop maxiter
    Rcpp::warning("brentZero: max iteration limit reached, no convergence.");
    // last best guess:
    root  = sb;
    froot = fsb;
  }

  double estimPrec = std::fabs(c - sb);
  // Estimate precision from the last bracket: typically |c - sb|
  // or if root is pinned at an endpoint, set NA similarly to uniroot
  if (root == lower || root == upper) {
    estimPrec = NA_REAL;
  }

  int finalIter = iterCount;

  // The total iteration is "main brent iteration + initSteps"
  int totalIter = finalIter + initSteps;

  // Return
  return List::create(
    _["root"] = root, _["f.root"] = froot,
    _["iter"] = totalIter, _["init.it"] = (initSteps > 0 ? initSteps : R_NaInt),
    _["estim.prec"] = estimPrec, _["exitcode"] = exitcode
  );
}
