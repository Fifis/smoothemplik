test_that("ctracelr constructs a correct path", {
  set.seed(1)
  xy <- matrix(rexp(200), ncol = 2)
  # ch <- chull(xy)
  # plot(xy)
  # points(xy[ch, ], pch = 16)
  # lines(c(0.5, 1.5), c(0.5, 1.5), lwd = 3, col = 2)
  trc <- ctracelr(xy, mu0 = c(0.5, 0.5), mu1 = c(1.5, 1.5), N = 10)
  expect_equal(nrow(trc), 11)
  expect_true(all(trc[, "exitcode"] == 0))
  expect_lt(max(trc[, "gradnorm"]), 1e-4)
})

