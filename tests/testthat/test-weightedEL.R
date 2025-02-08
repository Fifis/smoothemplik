test_that("weightedEL and cemplik are almost identical, but platform-dependent", {
  earth <- c(5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
             5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85)
  out1 <- cemplik(earth, mu = 5.517)[1:4]
  out2 <- weightedEL(earth, mu = 5.517, return.weights = TRUE)[1:4]
  expect_equal(names(out1), names(out2))

  # This is the platform-dependent part
  diff.lam <- max(abs(out1$lam - out2$lam))
  diff.wts <- max(abs(out1$wts - out2$wts))
  # Actually, these differences could be 0
  # https://github.com/Fifis/smoothemplik/actions/runs/9457471399/job/26051339706
  expect_lt(diff.lam, 5e-10)
  expect_lt(diff.wts, 5e-10)
})

test_that("input validation for weightedEL", {
  a <- -4:4
  expect_error(weightedEL(c(a, NA)), "Non-finite observations")
  expect_error(weightedEL(z = a, ct = 1/a), "Non-finite weights")
  expect_error(weightedEL(z = a, taylor.order = 10), "must be 2, 4")
  expect_error(weightedEL(z = a, taylor.order = 2, lower = 1, upper = -1), "must be less or equal")
  expect_error(weightedEL(z = cbind(-4:4, 1:9)), "Only one-dimensional")
  expect_warning(weightedEL(z = 1:10, mu = pi, ct = c(9:1, 1e-12), verbose = TRUE),
                 "Counts closer to 0")
  expect_error(weightedEL(a, ct = abs(a) * 1e-9), "Total weights")
})

test_that("weighted EL works as expected with good inputs", {
  a <- seq(-9, -1, 1)
  expect_true(is.finite(weightedEL(a, chull.fail = "taylor")$logelr))
  expect_true(!is.finite(weightedEL(a, chull.fail = "none")$logelr))
  expect_warning(weightedEL(a, chull.fail = "none", verbose = TRUE), "convex hull")
})

test_that("negative weights are handled correctly", {
  a <- -4:4
  f <- suppressWarnings(weightedEL(z = 1:6, ct = c(rep(1, 5), -0.9), mu = 3, return.weights = TRUE))
  expect_equal(f$logelr, -Inf)
  expect_true(is.finite(f$lam))
  expect_error(weightedEL(z = a, ct = a), "The total sum")
})

test_that("weight re-normalisation does not affect lambda", {
  expect_equal(weightedEL(z = 1:10, ct = 10:1, mu = pi)$lam,
               weightedEL(z = 1:10, ct = 10:1, mu = pi, SEL = TRUE)$lam, tolerance = 1e-14)
})

test_that("SEL weight scaling does not affect lamdda", {
  expect_equal(weightedEL(z = 1:10, ct = 10:1, mu = pi)$lam,
               weightedEL(z = 1:10, ct = 10:1, mu = pi, SEL = TRUE)$lam, tolerance = 1e-14)
})

test_that("names are preserved in weights", {
  expect_named(weightedEL(z = mtcars[, 1, drop = FALSE], mu = 20, return.weights = TRUE)$wts)
  expect_null(names(weightedEL(z = -4:3, return.weights = TRUE)$wts))
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
  EL0 <- weightedEL(z, ct = ct, return.weights = TRUE, weight.tolerance = 0)
  EL1 <- weightedEL(z, ct = ct, return.weights = TRUE)
  expect_equal(EL0$exitcode, 3)
  expect_equal(EL1$exitcode, 0)
  expect_equal(length(EL0$wts), length(EL1$wts))
  expect_equal(sum(EL1$wts == 0), 26) # If the defaults change, this will break
})

test_that("exit codes of weightedEL", {
  expect_equal(weightedEL(-4:3)$exitcode, 0)
  expect_equal(weightedEL(c(-1e-8, 1:9), chull.fail = "none")$exitcode, 2)
  expect_equal(weightedEL(-1:8, ct = c(1e-8, rep(1, 9)), weight.tolerance = 0)$exitcode, 3)
  expect_equal(weightedEL(1:5, chull.fail = "none")$exitcode, 5)
  expect_equal(weightedEL(0:3, chull.fail = "none")$exitcode, 7)
  expect_equal(weightedEL(rep(pi, 10), mu = pi, chull.fail = "none")$exitcode, 8)
  expect_equal(weightedEL(rep(pi, 10), mu = pi, chull.fail = "taylor")$exitcode, 8)
  expect_equal(weightedEL(0:3, chull.fail = "taylor")$exitcode, 9)
  expect_equal(weightedEL(c(0.999, 1:9), chull.fail = "taylor")$exitcode, 10)
  expect_equal(weightedEL(z = -1:8, lower = 0.1, upper = 0.1, taylor.order = 2)$exitcode, 13)
  # Come up with ideas for 1, &, 4 exit code!
})


