test_that("Trimmed weighted mean works as expected", {
  set.seed(1)
  z <- rnorm(10, sd = 2); z[which.max(z)] <- z[which.max(z)] * 100
  w <- pmin(1, 1 / abs(z)^2)  # Far-away observations tails get lower weight
  expect_identical(mean(z, trim = 0.20), trimmed.weighted.mean(z, trim = 0.20))
  expect_identical(weighted.mean(z, w), trimmed.weighted.mean(z, w = w))
  expect_lt(abs(trimmed.weighted.mean(z, trim = 0.20, w = w)), abs(mean(z))) # Weighted trimmed mean

  expect_identical(trimmed.weighted.mean(z, w = w, trim = 0.8), trimmed.weighted.mean(z, w = w, trim = 0.5))
  expect_identical(trimmed.weighted.mean(c(z[1:9], NA), w = w, trim = 0.05), NA_real_)
  expect_identical(trimmed.weighted.mean(NA, w = 1, trim = 0.05, na.rm = TRUE), NA_real_)

  # Handling errors
  expect_error(trimmed.weighted.mean(z, trim = 0.20, w = w[-1]), "the same length")
  expect_error(trimmed.weighted.mean(z, trim = 0.20, w = -w), "non-negative numbers")
  expect_error(trimmed.weighted.mean(z, trim = -0.20, w = w), "must be in")
})
