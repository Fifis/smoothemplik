test_that("Balanced EL preserves the sample mean", {
  x <- cbind(1:5, c(4, 7, 3, 9, 1))
  a <- EL(x, mu = c(0, 0), chull.fail = "balanced")
  x1 <- attr(a, "point1")
  x2 <- attr(a, "point2")
  expect_equal(colMeans(x), colMeans(rbind(x, x1, x2)), tolerance = 1e-12)
})

test_that("Adjusted2 EL is less aggressive than Adjusted", {
  x <- cbind(1:5, c(4, 7, 3, 9, 1))
  el1 <- EL(x, mu = c(0, 0), chull.fail = "adjusted")
  el2 <- suppressWarnings(EL(x, mu = c(0, 0), chull.fail = "adjusted2"))
  expect_gt(sum(attr(el1, "point1")^2), sum(attr(el2, "point1")^2))
})

test_that("ExEL can be invoked correctly from EL", {
  x <- cbind(1:5, c(4, 7, 3, 9, 1))
  el1 <- EL(x, mu = c(0, 0), chull.fail = "taylor")
  el2 <- EL(x, mu = c(0, 0), chull.fail = "wald")
  expect_gt(el2$logelr, el1$logelr)  # Taylor penalises more strongly
})
