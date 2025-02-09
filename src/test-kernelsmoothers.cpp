#include <testthat.h>
#include <RcppArmadillo.h>

extern arma::vec kernelFunCPP(arma::vec, std::string, int, bool);

context("C++ kernel functions") {
  test_that("2nd-order kernel take correct values") {
    double macheps = std::numeric_limits<double>::epsilon();
    arma::vec x = arma::linspace(-1.25, 1.25, 11);  // Step 0.25
    arma::vec kvals, ktrue;


    kvals = kernelFunCPP(x, "uniform", 2, false);
    ktrue = {0, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 5e-1, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "triangular", 2, false);
    ktrue = {0, 0, 2.5e-1, 5e-1, 7.5e-1, 1e0, 7.5e-1, 5e-1, 2.5e-1, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "epanechnikov", 2, false);
    ktrue = {0, 0, 3.28125e-1, 5.625e-1, 7.03125e-1, 7.5e-1, 7.03125e-1, 5.625e-1, 3.28125e-1, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "quartic", 2, false);
    ktrue = {0, 0, 1.79443359375e-1, 5.2734375e-1, 8.23974609375e-1, 9.375e-1, 8.23974609375e-1, 5.2734375e-1, 1.79443359375e-1, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "gaussian", 2, false);
    ktrue = {1.826490853890219e-1, 2.4197072451914335e-1, 3.011374321548044e-1, 3.5206532676429948e-1, 3.8666811680284921e-1, 3.9894228040143268e-1, 3.8666811680284921e-1, 3.5206532676429948e-1, 3.011374321548044e-1, 2.4197072451914335e-1, 1.826490853890219e-1};
     expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);
  }

  test_that("2nd-order convolution kernels take correct values") {
    double macheps = std::numeric_limits<double>::epsilon();
    arma::vec x = arma::linspace(-2.5, 2.5, 11);  // Step 0.5
    arma::vec kvals, ktrue;


    kvals = kernelFunCPP(x, "uniform", 2, true);
    ktrue = {0, 0, 1.25e-1, 2.5e-1, 3.75e-1, 5e-1, 3.75e-1, 2.5e-1, 1.25e-1, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "triangular", 2, true);
    ktrue = {0, 0, 2.0833333333333333e-2, 1.6666666666666667e-1, 4.7916666666666667e-1,
             6.6666666666666667e-1, 4.7916666666666667e-1, 1.6666666666666667e-1, 2.0833333333333333e-2, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "epanechnikov", 2, true);
    ktrue = {0, 0, 3.57421875e-2, 2.0625e-1, 4.587890625e-1, 6e-1, 4.587890625e-1,
             2.0625e-1, 3.57421875e-2, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "quartic", 2, true);
    ktrue = {0, 0, 8.5367475237165179e-3, 1.4369419642857143e-1, 4.906327383858817e-1,
             7.1428571428571429e-1, 4.906327383858817e-1, 1.4369419642857143e-1,
             8.5367475237165179e-3, 0, 0};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);

    kvals = kernelFunCPP(x, "gaussian", 2, true);
    ktrue = {5.9130280611822697e-2, 1.0377687435514868e-1, 1.6073276729880183e-1, 2.196956447338612e-1,
             2.6500353234402856e-1, 2.8209479177387814e-1, 2.6500353234402856e-1, 2.196956447338612e-1,
             1.6073276729880183e-1, 1.0377687435514868e-1, 5.9130280611822697e-2};
    expect_true(arma::abs(ktrue - kvals).max() <= 2 * macheps);
  }
}
