test_that("input validation for weightedEL0", {
  a <- -4:4
  expect_error(weightedEL0(c(a, NA)), "Non-finite observations")
  expect_error(weightedEL0(z = a, ct = 1/a), "Non-finite weights")
  expect_error(weightedEL0(z = cbind(-4:4, 1:9)), "Only one-dimensional")
  expect_warning(weightedEL0(z = 1:10, mu = pi, ct = c(9:1, 1e-12), verbose = TRUE),
                 "Counts closer to 0")
  expect_error(weightedEL0(a, ct = abs(a) * 1e-9), "Total weights")
})

test_that("weighted EL works as expected with good inputs", {
  a <- seq(-9, -1, 1)
  expect_true(is.finite(weightedEL0(a, chull.fail = "taylor")$logelr))
  expect_true(!is.finite(weightedEL0(a, chull.fail = "none")$logelr))
  expect_warning(weightedEL0(a, chull.fail = "none", verbose = TRUE), "convex hull")
})

test_that("negative weights are handled correctly", {
  a <- -4:4
  expect_error(weightedEL0(z = 1:6, ct = c(rep(1, 5), -0.9), mu = 3))
  expect_error(weightedEL0(z = a, ct = rep(0, 9)), "The total sum")
})

test_that("weight re-normalisation does not affect lambda", {
  expect_equal(weightedEL0(z = 1:10, ct = 10:1, mu = pi)$lam,
               weightedEL0(z = 1:10, ct = 10:1, mu = pi, SEL = TRUE)$lam, tolerance = 1e-14)
})

test_that("SEL weight scaling does not affect lamdda", {
  expect_equal(weightedEL0(z = 1:10, ct = 10:1, mu = pi)$lam,
               weightedEL0(z = 1:10, ct = 10:1, mu = pi, SEL = TRUE)$lam, tolerance = 1e-14)
})

test_that("names are preserved in weights", {
  expect_named(weightedEL0(z = mtcars[, 1, drop = FALSE], mu = 20, return.weights = TRUE)$wts)
  expect_null(names(weightedEL0(z = -4:3, return.weights = TRUE)$wts))
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
  EL0 <- weightedEL0(z, ct = ct, return.weights = TRUE, weight.tolerance = 0)
  EL1 <- weightedEL0(z, ct = ct, return.weights = TRUE)
  expect_identical(EL0$exitcode, 11L)  # A solution too close to the boundary
  expect_identical(EL1$exitcode, 0L)
  expect_identical(length(EL0$wts), length(EL1$wts))
  expect_identical(sum(EL1$wts == 0), 26L) # If the defaults change, this will break
})

test_that("weightedEL0 can handle ELR spanning condition failure", {
  x <- -4:5
  w <- 1:10
  mugrid <- seq(-6, 7, 0.25)
  ELR1 <- sapply(mugrid, function(m) -2*weightedEL0(z = x, ct = w, mu = m, chull.fail = "none")$logelr)
  ELR2 <- sapply(mugrid, function(m) -2*weightedEL0(z = x, ct = w, mu = m, chull.fail = "taylor")$logelr)
  ELR3 <- sapply(mugrid, function(m) -2*weightedEL0(z = x, ct = w, mu = m, chull.fail = "wald")$logelr)
  ELR4 <- sapply(mugrid, function(m) -2*weightedEL0(z = x, ct = w, mu = m, chull.fail = "adjusted")$logelr)
  ELR5 <- sapply(mugrid, function(m) -2*weightedEL0(z = x, ct = w, mu = m, chull.fail = "balanced")$logelr)
  expect_identical(ELR1[1], Inf)
  expect_true(all(is.finite(ELR2)))
  expect_true(all(is.finite(ELR3)))
  expect_true(all(is.finite(ELR4)))
  expect_true(all(is.finite(ELR5)))
  # plot(mugrid, ELR1, bty = "n", ylim = c(0, 600))
  # lines(mugrid, ELR2, col = 2)
  # lines(mugrid, ELR3, col = 4)
  # lines(mugrid, ELR4, col = 3)
  # lines(mugrid, ELR5, col = 5)
  # cbind(mugrid, ELR1, ELR2, ELR3, ELR4, ELR5)
})

test_that("exit codes of weightedEL0", {
  expect_identical(weightedEL0(-4:3)$exitcode, 0L)
  expect_identical(weightedEL0(c(-1e-8, 1:9), chull.fail = "none")$exitcode, 2L)
  expect_identical(weightedEL0(-1:8, ct = c(1e-8, rep(1, 9)), weight.tolerance = 0)$exitcode, 3L)
  expect_identical(weightedEL0(z = -1:8, ct = c(1e-15, rep(1, 9)), weight.tolerance = 0)$exitcode, 4L)
  expect_identical(weightedEL0(1:5, chull.fail = "none")$exitcode, 5L)
  expect_identical(weightedEL0(0:3, chull.fail = "none")$exitcode, 7L)
  expect_identical(weightedEL0(rep(pi, 10), mu = pi, chull.fail = "none")$exitcode, 8L)
  expect_identical(weightedEL0(rep(pi, 10), mu = pi, chull.fail = "taylor")$exitcode, 8L)
  expect_identical(weightedEL0(0:3, chull.fail = "taylor")$exitcode, 9L)
  expect_identical(weightedEL0(c(0.999, 1:9), chull.fail = "taylor")$exitcode, 10L)
  expect_identical(weightedEL0(1:5, chull.fail = "wald")$exitcode, 10L)
  expect_identical(weightedEL0(z = -1:8, ct = c(1e-12, rep(1, 8), 1e-12), weight.tolerance = 0)$exitcode, 11L)
  # Come up with ideas for exit code = 1!
})


