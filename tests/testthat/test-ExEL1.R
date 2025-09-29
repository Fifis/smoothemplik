test_that("Extrapolation works for uni-variate data", {
  z <- c(1, 4, 5, 5, 6, 6)
  ct <- 1:6
  e1 <- ExEL1(z, 0.5, ct = ct)
  e2 <- ExEL2(z, 0.5, ct = ct)
  expect_lt(e1, 0)
  expect_lt(e2, 0)

  # The values in the middle must coincide
  xseq <- seq(0, 7, 0.2)
  ctrl0 <- list(xlim = c(-1, 8))  # No extrapolation
  ctrl1 <- list(xlim = c(2.5, 5.5))
  e0 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl0)
  e1 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl1)
  expect_identical(sum(e0 == e1), 15L)

  # Root searching
  ctrl2 <- list(fmax = qchisq(0.999, 1))
  e2 <- ExEL1(z, xseq, ct = ct, exel.control = ctrl2)
  xlim <- attr(e2, "xlim")
  cond <- xseq > xlim[1] & xseq < xlim[2]
  expect_identical(e0[cond], e2[cond])
})

test_that("Extrapolation works with EL0 and EL1 similarly", {
  z <- c(1, 4, 5, 5, 6, 6)
  ct <- 1:6
  xseq <- seq(0.8, 1.4, length.out = 101)
  # e0 <- ExEL1(z, xseq, ct = ct, exel.control = list(xlim = c(0, 7)))
  e1 <- ExEL1(z, xseq, ct = ct, type = "EL0")
  e2 <- ExEL1(z, xseq, ct = ct, type = "EL1")
  expect_lt(max(abs(e1 - e2)), sqrt(.Machine$double.eps))
})

