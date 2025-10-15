test_that("Extrapolation works for uni-variate data", {
  z <- c(1, 4, 5, 5, 6, 6)
  ct <- 1:6
  e0 <- ExEL1(z, 0.5)
  e1 <- ExEL1(z, 0.5, ct = ct)
  e2 <- ExEL2(z, 0.5, ct = ct)
  expect_lt(e1, 0)
  expect_lt(e1, e2)  # Higher weights on larger obs --> worse for the opposite hypothesis
  expect_lt(e2, 0)

  # The values in the middle must coincide
  xseq <- seq(0, 7, 0.2)
  ctrl0 <- list(xlim = c(-1, 8))  # No extrapolation
  ctrl1 <- list(xlim = c(2.5, 5.5))
  e0 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl0)
  e1 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl1)
  expect_identical(sum(e0 == e1), 15L)

  # Root searching
  ctrl2 <- list(fmax = qchisq(0.999, 1))
  e2 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl2)
  xlim <- attr(e2, "xlim")
  cond <- xseq > xlim[1] & xseq < xlim[2]
  expect_identical(e0[cond], e2[cond])

})

test_that("Extrapolation works with EL0 and EL1 similarly", {
  z <- c(1, 4, 5, 5, 6, 6)
  ct <- 1:6
  xseq <- seq(0, 7, 0.2)
  # e0 <- ExEL1(z, xseq, ct = ct, exel.control = list(xlim = c(0, 7)))
  e1 <- ExEL1(z, xseq, ct = ct, type = "EL0")
  e2 <- ExEL1(z, xseq, ct = ct, type = "EL1")
  expect_lt(max(abs(e1 - e2)), sqrt(.Machine$double.eps))
})

test_that("Errors are caught in ExEL1", {
  set.seed(1)
  X <- matrix(rnorm(15), ncol = 3)
  expect_error(ExEL1(X, exel.control = list(xlim = c(-1, 1))), "only xlim='auto'")
  expect_error(ExEL1(X, mu = matrix(0, ncol = 1)), "mu must have d")
  expect_lt(ExEL1(X, mu = matrix(0, ncol = 3)), 0)
  Xbar <- colMeans(X)
  expect_identical(ExEL1(X, mu = Xbar), 0)
})

test_that("The cut-off size depends on the sample size in 1D ExEL", {
  expect_identical(attr(ExEL1(1:2, 0), "xlim"), c(1.05, 1.95))
  expect_identical(attr(ExEL1(1:3, 0), "xlim"), c(1.1, 2.9))
  expect_identical(attr(ExEL1(c(1:3, 5), 0), "xlim"), c(1.15, 4.85))
  expect_identical(attr(ExEL1(1:5, 0), "xlim"), c(1.1, 4.9))
  expect_identical(attr(ExEL1(1:7, 0), "xlim"), c(1.5, 6.5))
  expect_identical(attr(ExEL1(1:10, 0), "xlim"), c(2L, 9L))
})

