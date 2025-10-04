test_that("input validation for EL0", {
  a <- -4:4
  expect_error(EL0(c(a, NA)), "Non-finite observations")
  expect_error(EL0(z = a, ct = 1/a), "Non-finite weights")
  expect_error(EL0(z = cbind(-4:4, 1:9)), "Only one-dimensional")
  expect_warning(capture.output(EL0(z = 1:10, mu = pi, ct = c(9:1, 1e-12), verbose = TRUE)),
                 "Counts closer to 0")
  expect_error(EL0(a, ct = abs(a) * 1e-9), "Total weights")
})

test_that("EL works as expected with good inputs", {
  a <- seq(-9, -1, 1)
  ael <- suppressWarnings(EL0(a, log.control = list(order = 2)))
  expect_true(is.finite(ael$logelr))
  expect_true(ael$exitcode == 11)
  expect_false(ael$converged)
  expect_true(!is.finite(EL0(a)$logelr))
  expect_warning(EL0(a, verbose = TRUE), "in the convex hull")
})

test_that("negative weights are handled correctly", {
  a <- -4:4
  expect_error(EL0(z = 1:6, ct = c(rep(1, 5), -0.9), mu = 3))
  expect_error(EL0(z = a, ct = rep(0, 9)), "The total sum")
})

test_that("weight re-normalisation does not affect lambda", {
  expect_equal(EL0(z = 1:10, ct = 10:1, mu = pi)$lam,
               EL0(z = 1:10, ct = 10:1, mu = pi, renormalise = TRUE)$lam, tolerance = 1e-14)
})

test_that("SEL weight scaling does not affect lamdda", {
  expect_equal(EL0(z = 1:10, ct = 10:1, mu = pi)$lam,
               EL0(z = 1:10, ct = 10:1, mu = pi, renormalise = TRUE)$lam, tolerance = 1e-14)
})

test_that("names are preserved in weights", {
  expect_named(EL0(z = mtcars[, 1, drop = FALSE], mu = 20, return.weights = TRUE)$wts)
  expect_null(names(EL0(z = -4:3, return.weights = TRUE)$wts))
})

test_that("very small counts result in bad uniroot output", {
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
  # plot(ct, z, log = "x")
  EL0 <- EL0(z, ct = ct, return.weights = TRUE, weight.tolerance = 0)
  EL1 <- EL0(z, ct = ct, return.weights = TRUE)
  expect_identical(EL0$exitcode, 11L)  # A solution too close to the boundary
  expect_identical(EL1$exitcode, 0L)
  expect_identical(length(EL0$wts), length(EL1$wts))
  expect_identical(sum(EL1$wts == 0), 17L) # If the defaults change, this will break
})


test_that("exit codes of EL0", {
  expect_identical(EL0(-4:3)$exitcode, 0L)
  expect_identical(EL0(c(-1e-8, 1:9))$exitcode, 2L)
  expect_identical(EL0(-1:8, ct = c(1e-8, rep(1, 9)), weight.tolerance = 0)$exitcode, 3L)
  # expect_identical(EL0(z = -1:8, ct = c(1e-15, rep(1, 9)), weight.tolerance = 0)$exitcode, 4L)
  expect_identical(EL0(1:5)$exitcode, 5L)
  expect_identical(EL0(0:3)$exitcode, 7L)
  expect_identical(EL0(rep(pi, 10), mu = pi)$exitcode, 8L)
  # expect_identical(EL0(0:3, chull.fail = "taylor")$exitcode, 9L)
  # expect_identical(EL0(c(0.999, 1:9), chull.fail = "taylor")$exitcode, 10L)
  # expect_identical(EL0(1:5, chull.fail = "wald")$exitcode, 10L)
  expect_identical(EL0(z = -1:8, ct = c(1e-12, rep(1, 8), 1e-12), weight.tolerance = 0)$exitcode, 11L)
  # Come up with ideas for exit code = 1, 4, 9, 10!
})


