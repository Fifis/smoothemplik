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
