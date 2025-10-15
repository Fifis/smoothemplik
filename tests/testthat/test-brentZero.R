test_that("brentZero works well for good inputs", {
  f <- function (x) x - 1/3
  xroot <- brentZero(f, c(0, 1), tol = 0.0001)
  expect_lt(abs(xroot$root - 1/3), 16*.Machine$double.eps)
  expect_lt(abs(xroot$f.root), 16*.Machine$double.eps)
  sink(tempfile())
  xroot <- brentZero(f, c(-20, -19), extendInt = "upX", trace = 2)
  sink()
  expect_gt(xroot$init.it, 5)
  expect_lt(abs(xroot$f.root), 1e-13)
  expect_error(brentZero(f, c(-1, 1), extendInt = "something"), "Invalid 'extendInt' argument")

  f <- function (x) 1 - x
  sink(tempfile())
  xroot1 <- brentZero(f, c(9, 10), extendInt = "downX", trace = 2)
  xroot2 <- brentZero(f, c(9, 10), extendInt = "yes", trace = 2)
  xroot3 <- brentZero(f, c(9, 10), extendInt = "left", trace = 2)
  xroot4 <- brentZero(f, c(-10, -9), extendInt = "right", trace = 2)
  sink()
  expect_equal(xroot1$root, xroot2$root, tolerance = 1e-12)
  expect_equal(xroot1$root, xroot3$root, tolerance = 1e-12)
  expect_equal(xroot1$root, xroot4$root, tolerance = 1e-12)
  # Identical to 'left' and 'right'
  xroot5 <- brentZero(f, c(9, 10), extendInt = "lower")
  xroot6 <- brentZero(f, c(-10, -9), extendInt = "upper")
  expect_identical(xroot5, xroot3)
  expect_identical(xroot6, xroot4)

  expect_equal(unname(unlist(brentZero(function(x) x-1, c(0, 1)))), c(1, 0, 0, 0, 0, 0))
})

test_that("brentMin and brentZero correctly handle bad inputs", {
  expect_error(brentMin(sin, 1), "must be a vector of length 2")
  expect_error(brentZero(), "is missing, with no")
  expect_error(brentZero(sin), "must be strictly less")
  expect_error(brentZero(sin, interval = 1), "'interval' must be a vector of length 2")
  # TODO: there should be an error in expect_error(brentMin(sin, c(0, -1)), "must be strictly less")
  f <- function (x) x - 1
  expect_warning(brentZero(f, c(-20, -19), extendInt = "upX", maxiter = 2), "extension step limit")
  f <- function (x) log(x) - 5
  expect_warning(brentZero(f, c(1, 1000), maxiter = 3), "limit reached, no convergence")
  f <- function(x) c(x, 2*x)
  expect_error(brentZero(f, c(-1, 1), maxiter = 3), "Expecting a single value")
})

test_that("brentZero correctly handles bad functions", {
  # No root exists
  expect_warning(xroot <- brentZero(function(x) 1/x, c(-2, -1), extend = "left", maxiter = 100),
                 "Left-only extension step limit")
  expect_warning(xroot <- brentZero(function(x) 1/x, c(-2, -1), extend = "downX", maxiter = 100),
                 "Lower extension step limit")
  expect_identical(xroot$exitcode, 1L)

  # NA function
  f <- function(x) ifelse(x < 2, NA, x^2 - 1)
  sink(tempfile())
  expect_warning(xroot <- brentZero(f, c(3, 4), extend = "left", maxiter = 10, trace = 2),
                 "Left-only extension step")
  sink()
  f <- function(x) ifelse(x > -2, NA, x^2 - 1)
  expect_error(brentZero(f, c(-4, 3), extend = "right", maxiter = 10, trace = 2), "f\\(upper\\) is NA")
  sink(tempfile())
  expect_warning(xroot <- brentZero(f, c(-4, -3), extend = "right", maxiter = 10, trace = 2),
                 "Right-only extension step")
  sink()

  # Failure to extend in "yes" mode
  expect_warning(brentZero(function(x) x^2+1, c(-4, 5), extend = "yes", maxiter = 10),
                 "Extension step limit hit")

  # Shrink to finite
  f <- function(x) ifelse(abs(x) < 1, -1/(x-1)/(x+1) - 1000, NA)
  # curve(f, -2, 2)
  xroot <- brentZero(f, c(0.8, 0.98), maxiter = 100, extendInt = "right")
  expect_gt(xroot$init.it, 0)
  xroot <- brentZero(f, c(0.8, 0.95), maxiter = 100, extendInt = "yes")
  expect_gt(xroot$init.it, 0)


  # A zero function
  f <- function(x) return(0)
  xroot <- brentZero(f, c(-1, 1))
  expect_identical(xroot$estim.prec, 0)
})
