test_that("ctracelr constructs a correct path", {
  set.seed(1)
  xy <- matrix(rexp(200), ncol = 2)
  trc <- ctracelr(xy, mu0 = c(0.5, 0.5), mu1 = c(1.5, 1.5), N = 100)
  expect_equal(nrow(trc), 101)
  expect_true(all(trc[, "conv"] == 1))
<<<<<<< HEAD
  expect_true(all(trc[, "gnorm"] < 1e-7))
=======
  expect_lt(max(trc[, "gnorm"]), 1e-6)
>>>>>>> aa077fb (New unit tests for kernel and SEL functions)
})
