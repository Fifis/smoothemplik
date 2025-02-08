test_that("the rule-of-thumb bandwidth is invoked with a warning", {
  expect_warning(kernelWeights(1:10, PIT = TRUE), "No bandwidth supplied")
})

test_that("sparse matrices are identical to dense ones", {
  x   <- seq(-5, 5, 0.02)
  g   <- seq(-10, 10, 0.1)
  w1   <- kernelWeights(x, g, bw = 2, kernel = "triangular")
  w2   <- kernelWeights(x, g, bw = 2, kernel = "triangular", sparse = TRUE)
  expect_gt(object.size(w1), 3 * object.size(w2))
  expect_equal(as.vector(w1), as.vector(w2), tolerance = 1e-15)
})

test_that("the kernels have expected properties", {
  x <- seq(-1.50, 1.50, 0.05)
  expect_true(all(kernelWeights(x, xout = 0, bw = 1, kernel = "uniform", order = 2) >= 0))
  expect_true(all(kernelWeights(x, xout = 0, bw = 1, kernel = "uniform", order = 2, convolution = TRUE) >= 0))

  expect_true(all(range(kernelWeights(x, xout = 0, bw = 1, kernel = "uniform", order = 4)) == c(-0.4, 0.95)))
  r <- range(kernelWeights(x, xout = 0, bw = 1, kernel = "epanechnikov", order = 4))
  expect_true(r[1] < 0 & r[2] > 0)
  r <- range(kernelWeights(x, xout = 0, bw = 1, kernel = "quartic", order = 4))
  expect_true(r[1] < 0 & r[2] > 0)
  x <- seq(-5, 5, 0.1)
  r <- range(kernelWeights(x, xout = 0, bw = 1, kernel = "gaussian", order = 4))
  expect_true(r[1] < 0 & r[2] > 0)
})

