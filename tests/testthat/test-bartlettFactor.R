test_that("Bartlett factor works in dimension 1", {
  set.seed(1)
  x <- rchisq(50, df = 4)
  expect_gt(bartlettFactor(x), bartlettFactor(x, bias.adj = FALSE))
  expect_gt(bartlettFactor(x-4, centre = FALSE), 0)

  x <- matrix(c(-0.59, -0.35, 11.44, -0.45,  0.81,  2.41, -0.17, -0.05, -4.15, 0.08), ncol = 2)
  expect_warning(bartlettFactor(x), "cannot be negative; removing the adjustment")

  expect_warning(bartlettFactor(1:3), "requires n > 4")
  expect_error(bartlettFactor(cbind(1:5, 6:2)), "rank-deficient")
})


test_that("Bartlett factor works in multiple dimensions", {
  n <- 100
  g <- cbind(rchisq(n, 4)-4, rchisq(n, 3)-3, rchisq(n, 6)-6, rnorm(n), rexp(n)-1)
  expect_gt(bartlettFactor(g), 0)    # Bias-adjusted, centred
  expect_gt(bartlettFactor(g, centre = FALSE), 0)  # The true average was used in g
})

