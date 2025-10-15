test_that("Multi-variate ExEL2 extrapolation works", {
  set.seed(1)
  X <- cbind(rchisq(30, 3), rchisq(30, 3))
  ct <- runif(30)
  e1 <- ExEL1(X, mu = c(-1, -1),  ct = ct)
  e2 <- ExEL2(X, mu = c(-1, -1),  ct = ct)
  e3 <- ExEL2(X, mu = c(-1, -1))
  expect_lt(e1, 0)
  expect_lt(e2, 0)
  expect_lt(e3, e2)  # Because e2 weights are < 1 and e3 ct === 1

  # Comparing ExEL2 vs ExEL1 with bridges containing exp(x)
  z <- -4:4
  ct <- 9:1
  xseq <- seq(-7, 10.5, 0.1)
  xl <- range(xseq)
  a0 <- ExEL1(z, mu = xseq, ct = ct, exel.control = list(xlim = c(-11, 11)))
  a1 <- ExEL1(z, mu = xseq, ct = ct)
  a2 <- ExEL2(z, mu = xseq, ct = ct)
  expect_true(all(a1 >= a0))
  expect_true(all(a2 >= a0))

  expect_equal(ExEL2(z, 10, type = "EL0"), ExEL2(z, 10, type = "EL1"), tolerance = sqrt(.Machine$double.eps))
})

test_that("Errors are caught in ExEL2", {
  set.seed(1)
  X <- matrix(rnorm(15), ncol = 3)
  expect_error(ExEL2(X, exel.control = list(xlim = c(-1, 1))), "only xlim='auto'")
  expect_error(ExEL2(X, mu = matrix(0, ncol = 1)), "mu must have d")
  expect_lt(ExEL2(X, mu = matrix(0, ncol = 3)), 0)
  Xbar <- colMeans(X)
  expect_identical(ExEL2(X, mu = Xbar), 0)
})

test_that("Near-degenerate data sets are handled", {
  X <- cbind(1:5, 1:5)
  expect_lt(ExEL2(X, mu = c(0, 0)), 0)
})

test_that("Same-ray functionality works", {
  X <- matrix(c(-2:2, 1, 3, 7, 5, 4), ncol = 2)
  Xbar <- colMeans(X)
  v <- c(1, -1)
  tvec <- seq(0, 5, length.out = 11)
  mu.mat <- tvec * matrix(rep(v, 11), nrow = 11, byrow = TRUE)
  mu.mat <- sweep(mu.mat, 2, Xbar, "+")
  ex <- ExEL2(X, mu.mat)
  expect_identical(ex[1], 0)
  expect_true(all(ex[-1] < 0))
})


