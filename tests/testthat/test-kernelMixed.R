test_that("kernelMixedSmooth and kernelMixedDensity de-duplicate correctly", {
  set.seed(1)
  x  <- rnorm(1000)
  xout <- rnorm(500)
  xr <- round(x)
  xrout <- round(xout)
  w <- runif(1000, 1, 3)
  y  <- 1 + x^2 + rnorm(1000)
  by <- rep(1:4, each = 250)
  byout <- rep(1:4, each = 125)
  kd1 <- kernelMixedDensity(x = xr, by = by, weights = w,
                            xout = xrout, byout = byout, bw = 2)
  kd2 <- kernelMixedDensity(x = xr, by = by, weights = w,
                            xout = xrout, byout = byout, bw = 2, no.dedup = TRUE)
  expect_equal(kd1, kd2)
  km1 <- kernelMixedSmooth(x = xr, y = y, by = by, weights = w,
                           xout = xrout, byout = byout, bw = 2)
  km2 <- kernelMixedSmooth(x = xr, y = y, by = by, weights = w,
                           xout = xrout, byout = byout, bw = 2, no.dedup = TRUE)
  expect_equal(km1, km2)
})
