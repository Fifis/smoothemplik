test_that("LM via SVD works and recovers from failures", {
  X1 <- 1:10
  X2 <- X1
  X2[1] <- X2[1] + 0.000001
  X <- cbind(X1, X2)
  Y <- 1:10
  Y[10] <- Y[10] + 0.000001
  mod0 <- lm(Y ~ X)
  s <- svdlm(X, Y)
  expect_true(any(!is.finite(mod0$coefficients)))
  expect_true(all(is.finite(s)))

  # Setting up for failure
  X <- cbind(X1, X1)
  s <- svdlm(X, Y)
  expect_true(all(is.finite(s)))

  s <- svdlm(X, Y, rel.tol = 0)
  expect_true(all(is.finite(s)))

  # It takes a lot of effort to break SVD-LM
  X <- matrix(c(1e400, 0, 0, 0, 0, 1e-300), ncol = 2)
  Y <- c(1, 0, 1)
  expect_error(svdlm(X, Y), "pinv failed")

  # Vectors work with an internal conversion to matrix
  expect_identical(svdlm(1:10, 1:10), 1)
})
