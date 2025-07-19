test_that("EL and EL0 are almost identical, but platform-dependent", {
  earth <- c(5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
             5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85)
  out1 <- EL(earth, mu = 5.517, return.weights = TRUE)
  out2 <- EL0(earth, mu = 5.517, return.weights = TRUE)
  expect_true(setequal(intersect(names(out1), names(out2)), c("logelr", "lam", "wts", "exitcode", "iter")))


  # This is the platform-dependent part
  diff.lam <- max(abs(out1$lam - out2$lam))
  diff.wts <- max(abs(out1$wts - out2$wts))
  # Actually, these differences could be 0
  # https://github.com/Fifis/smoothemplik/actions/runs/9457471399/job/26051339706
  expect_lt(diff.lam, 5e-10)
  expect_lt(diff.wts, 5e-10)
})
