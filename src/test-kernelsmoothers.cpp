#include <testthat.h>
#include <RcppArmadillo.h>

extern arma::vec kernelFunCPP(arma::vec, std::string, int, bool);

context("C++ kernel functions") {
  test_that("kernel weights take correct values") {
    arma::vec x = arma::linspace(-1.25, 1.25, 11);  // Step 0.25
    std::string kernel = "epanechnikov";
    int order = 2;
    arma::vec k = kernelFunCPP(x, kernel, order, false);
    arma::vec ktrue = {0, 0, 21, 36, 45, 48, 45, 36, 21,  0,  0};
    ktrue /= 64.0;
    arma::vec diff = arma::abs(ktrue - k);
    expect_true(diff.max() == 0);
  }
}
