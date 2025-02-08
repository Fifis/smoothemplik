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

  # f <- function (x) 1 - x
  # xroot <- brentZero(f, c(9, 10), extendInt = "yes", trace = 2)
  # TODO: "downX" does not work here!

  expect_equal(unname(unlist(brentZero(function(x) x-1, c(0, 1)))), c(1, 0, 0, 0, 0))
})

test_that("brentMin correctly handles bad inputs", {
  expect_error(brentZero(), "is missing, with no")
  expect_error(brentZero(sin), "must be strictly less")
  # TODO: check if the interval is missing
  # TODO: there should be an error in expect_error(brentMin(sin, c(0, -1)), "must be strictly less")
  expect_error(brentMin(sin, 1), "must be a vector of length 2")
  f <- function (x) x - 1
  expect_error(brentZero(f, c(-20, -19), extendInt = "upX", maxiter = 2), "No sign change found")
  f <- function (x) log(x) - 5
  expect_warning(brentZero(f, c(1, 1000), maxiter = 3), "limit reached, no convergence")
})
