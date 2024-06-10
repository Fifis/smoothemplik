test_that("kernels integrate to 1, are symmetrical, and the order is correct", {
  ks <- c("gaussian", "uniform", "triangular", "epanechnikov", "quartic")
  os <- c(2, 4, 6)
  convs <- c(FALSE, TRUE)
  par.grid <- expand.grid(kernel = ks, order = os, convolution = convs, stringsAsFactors = FALSE)
  m <- t(sapply(seq_len(nrow(par.grid)), function(i) {
    f <- function(x) kernelFun(x, kernel = par.grid$kernel[i],
                               order = par.grid$order[i], convolution = par.grid$convolution[i], rescale = FALSE)
    lims <- if (par.grid$convolution[i]) c(-2, 2) else c(-1, 1)
    if (par.grid$kernel[i] == "gaussian") lims <- c(-Inf, Inf)
    rt <- 2e-13
    m0 <- integrate(f, lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m1 <- integrate(function(x) x*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m2 <- integrate(function(x) x^2*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m3 <- integrate(function(x) x^3*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m4 <- integrate(function(x) x^4*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m5 <- integrate(function(x) x^5*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    m6 <- integrate(function(x) x^6*f(x), lims[1], lims[2], rel.tol = rt, subdivisions = 200)$value
    c(m0, m1, m2, m3, m4, m5, m6)
  }))
  expect_equal(m[, 1], rep(1, nrow(m)), tolerance = 1e-14)
  expect_equal(m[, 2], rep(0, nrow(m)), tolerance = 1e-14)
  expect_equal(m[, 4], rep(0, nrow(m)), tolerance = 1e-14)
  expect_equal(m[, 6], rep(0, nrow(m)), tolerance = 1e-14)
  # Higher-order kernels have zero variance
  expect_equal(m[par.grid$order > 2, 3], rep(0, sum(par.grid$order > 2)), tolerance = 1e-14)
  expect_equal(m[par.grid$order > 4, 5], rep(0, sum(par.grid$order > 4)), tolerance = 1e-13)
  expect_true(all(m[par.grid$order <= 2, 3] > 1e-8))
  expect_true(all(abs(m[par.grid$order <= 4, 5]) > 1e-8))
  expect_true(all(abs(m[, 7]) > 1e-8))
})
