test_that("Balanced EL preserves the sample mean", {
  set.seed(1)
  a <- EL(x, mu = c(0, 0), chull.fail = "balanced")
  x1 <- attr(a, "point1")
  x2 <- attr(a, "point2")
  expect_equal(colMeans(x), colMeans(rbind(x, x1, x2)), tolerance = 1e-12)
})
