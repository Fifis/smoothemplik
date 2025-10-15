test_that("Damped Newton search works as expected", {
  f1 <- function(x) list(fn = x - log(x), gradient = 1 - 1/x, Hessian = matrix(1/x^2, 1, 1))
  o1 <- optim(2, function(x) f1(x)[["fn"]], gr = function(x) f1(x)[["gradient"]], method = "BFGS")
  o2 <- dampedNewton(f1, 2)
  expect_equal(o1$par, o2$par, tolerance = 1e-6)

  # The minimum of f3 should be roughly at -0.57
  f3 <- function(x) list(fn = sum(exp(x) + 0.5 * x^2), gradient = exp(x) + x, Hessian =  diag(exp(x) + 1))
  o3 <- dampedNewton(f3, seq(0.1, 5, length.out = 11))
  expect_identical(o3$convergence, 0L)
})
