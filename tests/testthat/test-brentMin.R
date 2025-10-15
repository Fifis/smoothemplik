test_that("brentMin works well for good inputs", {
  f <- function (x) (x - 1/3)^2
  xmin <- brentMin(f, c(0, 1), tol = 0.0001)
  expect_lt(abs(xmin$root - 1/3), 16*.Machine$double.eps)
  expect_lt(abs(xmin$f.root), 16*.Machine$double.eps)
  sink(tempfile())
  expect_gt(brentMin(f, c(-20, -20+1e-4), trace = 2)$iter, 8)
  sink()
})

test_that("brentMin correctly handles bad inputs", {
  expect_error(brentMin(), "is missing, with no")
  expect_error(brentMin(sin), "must be strictly less")
  # TODO: check if the interval is missing
  # TODO: there should be an error in expect_error(brentMin(sin, c(0, -1)), "must be strictly less")
  expect_error(brentMin(sin, 1), "must be a vector of length 2")
  expect_error(brentMin(function(x) -cos(x), c(-0.5, 0.4), tol = -1), "must be > 0 and finite")
  suppressWarnings(expect_warning(brentMin(function(x) NULL, c(-0.5, 0.4), maxiter = 1), "Function returned NULL"))
  suppressWarnings(expect_warning(brentMin(function(x) NA, c(-0.5, 0.4), maxiter = 1), "Function returned NA"))
})
