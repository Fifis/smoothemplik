test_that("tlog works", {
  x <- seq(1.5, 3.5, 0.1)
  expect_equal(tlog(x, a = 2, k = 1), x/2 - 1 + log(2), tolerance = 1e-15)

  # 2nd derivative = -1/x^2
  expect_equal(tlog(x, a = 3, k = 5, d = 2),
               -1/9 + 2*(x-3)/27 - 1/27*(x-3)^2 + 4/243*(x-3)^3, tolerance = 1e-15)
})

test_that("logTaylor works", {
  x <- seq(0.05, 0.95, length.out = 9) # Cut-off: 1/9 = 0.111
  l1 <- logTaylor(x)
  l2 <- log(x)
  expect_identical(l1[-1], l2[-1])
  expect_gt(abs(l1[1] - l2[1]), 0)
})
