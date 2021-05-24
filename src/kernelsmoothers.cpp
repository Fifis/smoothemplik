#include <Rcpp.h>
using namespace Rcpp;

double dstdnorm(double x) {
  return 0.39894228040143267794 * std::exp(-0.5 * x * x);
}

// Epanechnikov divided by sqrt(5)
double depan(double x) {
  if (std::abs(x) < 2.23606797749978969641) {
    return 0.33541019662496845446 - 0.067082039324993690892 * x * x;
  } else {
    return 0;
  }
}

// In the innermost loop, using ksum += dnorm((xgrid[i] - x[j]) / bw, 0, 1, 0) gives terrible performance
// [[Rcpp::export]]
NumericMatrix kernelWeightsOneCPP(NumericVector x, NumericVector xgrid, double bw, bool gaussian = true) {
  int i, j;
  int n = xgrid.size();
  int nx = x.size();
  NumericMatrix kw(n, nx);

  if (gaussian) {
    for (i=0; i < n; i++) {
      NumericVector thisrow(nx);
      for (j=0; j < nx; j++) {
        thisrow[j] = dstdnorm((xgrid[i] - x[j]) / bw);
      }
      kw(i, _) = thisrow;
    }
  } else {
    for (i=0; i < n; i++) {
      NumericVector thisrow(nx);
      for (j=0; j < nx; j++) {
        thisrow[j] = depan((xgrid[i] - x[j]) / bw);
      }
      kw(i, _) = thisrow;
    }
  }

  return kw;
}

// [[Rcpp::export]]
NumericVector kernelDensityOneCPP(NumericVector x, NumericVector xgrid, double bw, bool gaussian = true) {
  int i, j;
  int n = xgrid.size();
  int nx = x.size();
  NumericVector out(n);

  double nb = (double)nx * bw;

  NumericMatrix kw = kernelWeightsOneCPP(x, xgrid, bw, gaussian);

  for (i=0; i < n; i++) {
    double ksum = 0;
    NumericVector krow = kw(i, _);
    for (j=0; j < nx; j++) {
      ksum += krow[j];
    }
    out[i] = ksum / nb;
  }

  return out;
}

// [[Rcpp::export]]
NumericVector kernelSmoothOneCPP(NumericVector x, NumericVector y, NumericVector xgrid, double bw, bool gaussian = true, bool LOO = false) {
  int i, j;
  int n = xgrid.size();
  int nx = x.size();
  NumericVector out(n);
  NumericMatrix kw = kernelWeightsOneCPP(x, xgrid, bw, gaussian);

  // LOO: setting diagonal elements to zero, assuming x = xgrid! (R makes sure it happens, though.)
  if (LOO) {
    for (i=0; i < n; i++) {
      kw(i, i) = 0;
    }
  }

  for (i=0; i < n; i++) {
    double ysum = 0;
    double ksum = 0;
    NumericVector krow = kw(i, _);
    for (j=0; j < nx; j++) {
      ysum += y[j]*krow[j];
      ksum += krow[j];
    }
    out[i] = ysum / ksum;
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix kernelWeightsMultiCPP(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian = true) {
  int i, j, k;
  int nx = x.nrow(), d = x.ncol(), n = xgrid.nrow();

  NumericMatrix pk(n, nx); // Product kernel matrix
  NumericVector out(n);
  pk.fill(1.0);

  double nb = (double)nx;
  for (k=0; k < d; k++) {
    nb *= bw[k];
  }

  if (gaussian) {
    for (k=0; k < d; k++) {
      for (i=0; i < n; i++) {
        for (j=0; j < nx; j++) {
          pk(i, j) *= dstdnorm((xgrid(i, k) - x(j, k)) / bw[k]);
        }
      }
    }
  } else {
    for (k=0; k < d; k++) {
      for (i=0; i < n; i++) {
        for (j=0; j < nx; j++) {
          pk(i, j) *= depan((xgrid(i, k) - x(j, k)) / bw[k]);
        }
      }
    }
  }

  return pk;
}

// [[Rcpp::export]]
NumericVector kernelDensityMultiCPPold(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian = true) {
  int i, j, k;
  int nx = x.nrow(), d = x.ncol(), n = xgrid.nrow();
  NumericVector out(n);

  double nbs = (double)nx;
  for (k=0; k < d; k++) {
    nbs *= bw[k];
  }

  if (gaussian) {
    for (i=0; i < n; i++) {
      double ksum = 0;
      for (j=0; j < nx; j++) {
        NumericVector dx = xgrid(i,_) - x(j,_);
        double kprod = 1.0;
        for (k=0; k < d; k++) {
          kprod *= dstdnorm(dx[k] / bw[k]);
        }
        ksum += kprod;
      }
      out[i] = ksum / nbs;
    }
  } else {
    for (i=0; i < n; i++) {
      double ksum = 0;
      for (j=0; j < nx; j++) {
        NumericVector dx = xgrid(i,_) - x(j,_);
        double kprod = 1.0;
        for (k=0; k < d; k++) {
          kprod *= depan(dx[k] / bw[k]);
        }
        ksum += kprod;
      }
      out[i] = ksum / nbs;
    }
  }

  return out;
}

// [[Rcpp::export]]
NumericVector kernelDensityMultiCPP(NumericMatrix x, NumericMatrix xgrid, NumericVector bw, bool gaussian = true) {
  int i, j, k;
  int nx = x.nrow(), d = x.ncol(), n = xgrid.nrow();

  NumericVector out(n);

  double nb = (double)nx;
  for (k=0; k < d; k++) {
    nb *= bw[k];
  }

  NumericMatrix kw = kernelWeightsMultiCPP(x, xgrid, bw, gaussian);

  for (i=0; i < n; i++) {
    double ksum = 0;
    for (j=0; j < nx; j++) {
      ksum += kw(i, j);
    }
    out[i] = ksum / nb;
  }

  return out;
}


// [[Rcpp::export]]
NumericVector kernelSmoothMultiCPP(NumericMatrix x, NumericVector y, NumericMatrix xgrid, NumericVector bw, bool gaussian = true, bool LOO = false) {
  int i, j;
  int nx = x.nrow(), n = xgrid.nrow();
  NumericVector out(n);
  NumericMatrix kw = kernelWeightsMultiCPP(x, xgrid, bw, gaussian);

  // LOO: setting diagonal elements to zero, assuming x = xgrid! (R makes sure it happens, though.)
  if (LOO) {
    for (i=0; i < n; i++) {
      kw(i, i) = 0;
    }
  }

  for (i=0; i < n; i++) {
    double ysum = 0;
    double ksum = 0;
    NumericVector krow = kw(i, _);
    for (j=0; j < nx; j++) {
      ysum += y[j]*krow[j];
      ksum += krow[j];
    }
    out[i] = ysum / ksum;
  }

  return out;
}

