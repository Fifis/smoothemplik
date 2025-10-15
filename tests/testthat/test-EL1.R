test_that("EL0 and EL1 are almost identical, but platform-dependent", {
  earth <- c(5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
             5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85)
  out0 <- EL0(earth, mu = 5.517, return.weights = TRUE)
  out1 <- EL1(earth, mu = 5.517, return.weights = TRUE)
  expect_true(setequal(intersect(names(out0), names(out1)), c("logelr", "lam", "wts", "iter", "deriv", "exitcode")))


  # This is the platform-dependent part
  diff.lam <- max(abs(out0$lam - out1$lam))
  diff.wts <- max(abs(out0$wts - out1$wts))
  # Actually, these differences could be 0
  # https://github.com/Fifis/smoothemplik/actions/runs/9457471399/job/26051339706
  expect_lt(diff.lam, 5e-10)
  expect_lt(diff.wts, 5e-10)
})

test_that("EL1 derivatives behave adequately", {
  # The derivative is zero at the mean
  set.seed(1)
  X <- matrix(rnorm(20), ncol = 2)
  Xbar <- colMeans(X)
  expect_equal(EL1(X, mu = Xbar, deriv = TRUE)$deriv[1], 0, tolerance = 1e-12)
})
