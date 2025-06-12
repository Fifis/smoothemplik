test_that("input validation for weightedEuL", {
  a <- -4:4
  expect_error(weightedEuL(c(a, NA)), "Non-finite observations")
  expect_error(weightedEuL(z = a, ct = 1/a), "Non-finite weights")
  expect_error(weightedEuL(z = a, vt = 1/a), "Non-finite variance weights")
  expect_warning(weightedEuL(z = 1:10, mu = pi, ct = c(9:1, 1e-12), verbose = TRUE),
                 "Counts closer to 0")
  expect_warning(weightedEuL(z = 1:10, mu = pi, vt = c(9:1, 1e-12), verbose = TRUE),
                 "Variance weights closer to 0")
  expect_error(weightedEuL(a, ct = rep(-1, 9)), "total sum of weights")
})

test_that("weighted EL works as expected with good inputs", {
  # Multi-variate case
  a <- cbind(seq(-9, -1, 1), c(3, 4, 5, 1, 2, 3, 6:8))
  expect_length(weightedEuL(a)$lam, 2)
})

test_that("negative weights are handled well", {
  a <- -4:4
  f <- weightedEuL(z = 1:6, ct = c(rep(1, 5), -0.9), mu = 3, return.weights = TRUE)
  expect_true(is.finite(f$lam))
  expect_error(weightedEuL(z = a, ct = a), "The total sum")
})

test_that("The FOC value at the optimum is near-zero.", {
  expect_lt(abs(weightedEuL(z = 1:10, ct = 10:1, mu = pi)$f.root), 8*.Machine$double.eps)
})

test_that("names are preserved in weights", {
  expect_named(weightedEuL(z = mtcars[, 1, drop = FALSE], mu = 20, return.weights = TRUE)$wts)
  expect_null(names(weightedEuL(z = -4:3, return.weights = TRUE)$wts))
})

test_that("very small counts are handled reasonably well", {
  z <- c(0.31, 0.15, -0.25, 0.14, -0.56, 0.71, 1.03, -0.19, -0.56, 0.31, -0.08,
         1.45, -0.02, 0.44, 0.02, -0.52, 0.13, -1.3, 1.06, 0.11, 1.62, 0.36,
         -0.53, 0.47, -0.76, -1.1, 0.29, -0.45, 0, 0.08, -0.62, -0.63, -0.16,
         1.4, -1.83, 0.73, 0.44, 1.44, -0.42, 0.51, 0.37, -0.79, 1.9, 1.87, 1.29, 2.99, 1.3, -3.42)
  ct <- c(4.2e-01, 3.7e-01, 1.1e-01, 7.9e-02, 4.5e-03, 4.1e-03, 1.9e-03, 1.6e-03,
          1.0e-03, 1.0e-03, 3.2e-04, 1.9e-04, 1.6e-04, 7.3e-05, 4.5e-05, 1.9e-05,
          1.7e-05, 1.1e-05, 1.0e-05, 6.8e-06, 6.6e-06, 6.4e-06, 5.8e-06, 4.3e-06,
          1.6e-06, 4.9e-07, 8.9e-08, 5.8e-08, 4.3e-08, 4.2e-08, 3.0e-08, 1.2e-08,
          5.0e-09, 3.9e-09, 3.1e-09, 2.1e-09, 7.6e-10, 4.3e-10, 3.0e-10, 2.8e-10,
          2.3e-10, 1.3e-10, 3.1e-11, 2.1e-11, 1.9e-12, 1.3e-12, 2.8e-14, 2.0e-15)
  vt <- ct
  EL0 <- weightedEuL(z, ct = ct, vt = vt, return.weights = TRUE, weight.tolerance = 0)
  EL1 <- weightedEuL(z, ct = ct, vt = vt, return.weights = TRUE)
  expect_equal(length(EL0$wts), length(EL1$wts))
  expect_equal(sum(EL1$wts == 0), 16) # If the defaults change, this will break
})

test_that("exit codes of weightedEuL", {
  expect_equal(weightedEuL(-4:3)$exitcode, 0)
  expect_equal(weightedEuL(1:5, chull.diag = FALSE)$exitcode, 0)
  expect_equal(weightedEuL(1:5, chull.diag = TRUE)$exitcode, 1)
  expect_equal(weightedEuL(matrix(1:2, nrow = 1))$exitcode, 2)
  expect_equal(weightedEuL(cbind(1:4, 2:5))$exitcode, 3)
})
