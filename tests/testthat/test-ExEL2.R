test_that("Multi-variate ExEL2 extrapolation works", {
  set.seed(1)
  X <- cbind(rchisq(30, 3), rchisq(30, 3))
  ct <- runif(30)
  e1 <- ExEL1(X, mu = c(-1, -1),  ct = ct)
  e2 <- ExEL2(X, mu = c(-1, -1),  ct = ct)
  expect_lt(e1, 0)
  expect_lt(e2, 0)
})

