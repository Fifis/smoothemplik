test_that("the SEL weights for various kernels", {
  expect_equal(getSELWeights((1:5)/2, bw = 1, kernel = "uniform")[[2]],
               list(idx = 1:4, ct = rep(1, 4)/4))
  expect_equal(getSELWeights((1:5)/2, bw = 1, kernel = "triangular")[[2]],
               list(idx = 1:3, ct = c(1, 2, 1)/4))
  expect_equal(getSELWeights((1:5)/2, bw = 1, kernel = "epanechnikov")[[2]],
               list(idx = 1:3, ct = c(0.3, 0.4, 0.3)))
})
