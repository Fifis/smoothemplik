#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// Uniform kernel, unscaled (with roughness = int x^2 f(x) dx = 1/3)
arma::vec kuniform(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (int i=0; i < ax.n_elem; i++) {
    if (ax[i] < 1.0) {
      ax[i] = 0.5;
    } else {
      ax[i] = 0;
    }
  }
  return ax;
}

// Triangular kernel, unscaled (with roughness = int x^2 f(x) dx = 1/6)
arma::vec ktriangular(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (int i=0; i < ax.n_elem; i++) {
    if (ax[i] < 1.0) {
      ax[i] = 1.0 - ax[i];
    } else {
      ax[i] = 0;
    }
  }
  return ax;
}

// Epanechnikov kernel, unscaled (with roughness = int x^2 f(x) dx = 1/5)
arma::vec kepanechnikov(arma::vec x) {
  arma::vec ax = arma::abs(x);
  for (int i=0; i < ax.n_elem; i++) {
    if (ax[i] < 1.0) {
      ax[i] = 0.75 * (1.0 - ax[i] * ax[i]);
    } else {
      ax[i] = 0;
    }
  }
  return ax;
}

// Gaussian kernel unscaled (with roughness 1)
arma::vec kgaussian(arma::vec x) {
  return 0.398942280401432678 * arma::exp(-0.5 * (x % x));
}

// In the innermost loop, using ksum += dnorm((xgrid[i] - x[j]) / bw, 0, 1, 0) gives terrible performance
// [[Rcpp::export]]
arma::mat kernelWeightsOneCPP(arma::vec x, arma::vec xgrid, double bw, std::string kernel = "gaussian") {
  int i;
  int n = xgrid.n_elem;
  int nx = x.n_elem;
  arma::mat kw(n, nx);
  arma::vec xs = x / bw;
  arma::vec gs = xgrid / bw; // Scaling by the bandwidth

  if (kernel == "gaussian") {
    for (i=0; i < n; i++) {
      kw.row(i) = kgaussian(gs[i] - xs).t();
    }
  } else if (kernel == "uniform") {
    for (i=0; i < n; i++) {
      kw.row(i) = kuniform(gs[i] - xs).t();
    }
  } else if (kernel == "triangular") {
    for (i=0; i < n; i++) {
      kw.row(i) = ktriangular(gs[i] - xs).t();
    }
  } else if (kernel == "epanechnikov") {
    for (i=0; i < n; i++) {
      kw.row(i) = kepanechnikov(gs[i] - xs).t();
    }
  }

  return kw;
}

// [[Rcpp::export]]
arma::mat kernelWeightsCPP(arma::mat x, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian") {
  int d = x.n_cols;
  // The product kernel matrix starts with the first dimension (there is at least one column or row)
  arma::mat pk = kernelWeightsOneCPP(x.col(0), xgrid.col(0), bw[0], kernel);
  if (d > 1) { // We need to compute the product kernel starting from the 2nd till the last dimension
    for (int k=1; k < d; k++) {
      pk %= kernelWeightsOneCPP(x.col(k), xgrid.col(k), bw[k], kernel);
    }
  }
  return pk;
}

// [[Rcpp::export]]
Rcpp::NumericVector kernelDensityCPP(arma::mat x, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian") {
  int d = x.n_cols;

  double nb = (double)x.n_rows; // n*prod(b) in the denominator
  for (int k=0; k < d; k++) {
    nb *= bw[k];
  }

  arma::mat kw = kernelWeightsCPP(x, xgrid, bw, kernel);
  arma::vec out = sum(kw, 1) / nb;
  Rcpp::NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

// [[Rcpp::export]]
Rcpp::NumericVector kernelSmoothCPP(arma::mat x, arma::vec y, arma::mat xgrid, arma::vec bw, std::string kernel = "gaussian", bool LOO = false) {
  arma::mat kw = kernelWeightsCPP(x, xgrid, bw, kernel);

  // LOO: setting diagonal elements to zero, assuming x = xgrid! (R makes sure it happens, though.)
  if (LOO) {
    for (int i=0; i < x.n_rows; i++) {
      kw(i, i) = 0;
    }
  }

  arma::vec ksum = sum(kw, 1); // Nadaraya--Watson denominator
  kw.each_row() %= y.t(); // Nadaraya--Watson numerator: y_i * w_ij (in place to save memory)
  arma::vec ysum = sum(kw, 1);
  arma::vec out = ysum / ksum;
  Rcpp::NumericVector rout = Rcpp::NumericVector(out.begin(), out.end());
  return rout;
}

