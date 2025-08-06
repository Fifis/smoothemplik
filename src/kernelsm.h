#ifndef KERNEL_DENSITY_H
#define KERNEL_DENSITY_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

arma::mat kernelWeightsCPP(arma::mat x, arma::mat xout, arma::mat bw, std::string kernel, int order, bool convolution);

struct KernelDensityWorker : public RcppParallel::Worker {
  const arma::mat x; const arma::mat xout;
  const arma::vec weights;
  const arma::mat bw;
  const std::string kernel; const int order; const bool convolution; const arma::uvec starts;
  const arma::uvec ends;
  arma::vec& out;
  const arma::vec nb;

  KernelDensityWorker(const arma::mat x, const arma::mat xout, const arma::vec weights,
                      const arma::mat bw, const std::string kernel,
                      const int order, const bool convolution,
                      const arma::uvec starts, const arma::uvec ends,
                      arma::vec& out, const arma::vec nb)
    : x(x), xout(xout), weights(weights), bw(bw), kernel(kernel), order(order), convolution(convolution), starts(starts), ends(ends), out(out), nb(nb) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      arma::mat xsubgrid = xout.rows(starts[j], ends[j]);
      arma::mat bwsub = bw.rows(starts[j], ends[j]);
      arma::mat kw = kernelWeightsCPP(x, xsubgrid, bwsub, kernel, order, convolution);
      kw.each_row() %= weights.t();
      out.subvec(starts[j], ends[j]) = arma::sum(kw, 1) / nb.subvec(starts[j], ends[j]);
    }
  }
};

Rcpp::NumericVector kernelDensityCPP(arma::mat x, arma::mat xout, arma::vec weights, arma::mat bw, std::string kernel, int order, bool convolution, int chunks);
#endif // KERNEL_DENSITY_H

#ifndef KERNEL_SMOOTH_H
#define KERNEL_SMOOTH_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>

// arma::mat kernelWeightsCPP(arma::mat x, arma::mat xout, arma::vec bw, std::string kernel, int order, bool convolution);

struct KernelSmoothWorker : public RcppParallel::Worker {
  const arma::mat x;
  const arma::vec y;
  const arma::mat xout;
  const arma::vec weights;
  const arma::mat bw;
  const std::string kernel;
  const int order;
  const bool convolution;
  const bool LOO;
  const arma::uvec starts;
  const arma::uvec ends;
  arma::vec& ksum;
  arma::vec& ysum;

  KernelSmoothWorker(const arma::mat x, const arma::vec y, const arma::mat xout,
                     const arma::vec weights, const arma::mat bw, const std::string kernel,
                     const int order, const bool convolution, const bool LOO,
                     const arma::uvec starts, const arma::uvec ends,
                     arma::vec& ksum, arma::vec& ysum)
    : x(x), y(y), xout(xout), weights(weights), bw(bw), kernel(kernel), order(order), convolution(convolution), LOO(LOO), starts(starts), ends(ends), ksum(ksum), ysum(ysum) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; j++) {
      arma::mat xsubgrid = xout.rows(starts[j], ends[j]);
      arma::mat bwsub = bw.rows(starts[j], ends[j]);
      arma::mat kw = kernelWeightsCPP(x, xsubgrid, bwsub, kernel, order, convolution);
      kw.each_row() %= weights.t();
      if (LOO) { // LOO: being careful because this matrix is not square
        for (arma::uword k=starts[j]; k <= ends[j]; k++) { // Note the <= !
          kw(k-starts[j], k) = 0;
        }
      }
      ksum.subvec(starts[j], ends[j]) = arma::sum(kw, 1);
      kw.each_row() %= y.t();
      ysum.subvec(starts[j], ends[j]) = arma::sum(kw, 1);
    }
  }
};

Rcpp::NumericVector kernelSmoothCPP(arma::mat x, arma::vec y, arma::mat xout, arma::vec weights, arma::mat bw, std::string kernel, int order, bool convolution, int chunks);
#endif // KERNEL_SMOOTH_H

