test_that("kernelDiscreteDensitySmooth works", {
  set.seed(1)
  x <- sort(rnorm(400))
  y <- as.numeric(runif(400) < 0.5*pnorm(x) + 0.25)
  ux <- -4:3+0.5
  x2 <- as.numeric(as.character(cut(x, -4:4, labels = ux)))
  yhat.comp <- kernelDiscreteDensitySmooth(x = x2, y = y, compact = TRUE)
  z <- prepareKernel(x = x2, y = y, bw = 1)
  wm <- sapply(ux[2:(length(ux)-1)], function(xx) {
    ii <- z$x == xx
    weighted.mean(z$y[ii], w = z$weights[ii])
  })
  den <- sapply(ux[2:(length(ux)-1)], function(xx) sum(z$weights[z$x == xx] / sum(z$weights)))
  expect_equal(as.numeric(yhat.comp$fhat), den)
  expect_equal(yhat.comp$y, wm)
})
