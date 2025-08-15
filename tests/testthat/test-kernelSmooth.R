test_that("kernelSmooth correctly handles bad inputs", {
  expect_error(kernelSmooth(c(1:9, NA), 1:10), "finite numeric values only")
  expect_warning(suppressMessages(kernelSmooth(c(1:9, NA), 1:10, bw = 2)),
                 "Some smoothed values are NA or Inf")
  expect_warning(suppressMessages(kernelSmooth(c(1:9, NA), 1:10, bw = 2, kernel = "triangular")),
                 "Some smoothed values are NaN")
})

test_that("kernelSmooth linear and higher-order smoothers address boundary bias", {
  set.seed(1)
  x <- sort(runif(300, -6, 6))
  g <- seq(-6, 6, 0.1)
  f <- function(x) 1 + x + sin(x)
  y <- f(x) + rt(300, df = 4)

  m2lc <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", no.dedup = TRUE)
  m4lc <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", order = 4, no.dedup = TRUE)
  m2ll <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", degree = 1, no.dedup = TRUE)
  expect_gt(mean(abs(f(g) - m2lc)), mean(abs(f(g) - m4lc)))
  expect_gt(mean(abs(f(g) - m2lc)), mean(abs(f(g) - m2ll)))
})

test_that("de-duplication works in kernelSmooth", {
  set.seed(1)
  n.uniq <- 100
  n <- 1000
  inds <- sort(ceiling(runif(n, 0, n.uniq)))
  x.uniq <- sort(rnorm(n.uniq))
  y.uniq <- 1 + x.uniq + sin(x.uniq*2) + rnorm(n.uniq)
  x <- x.uniq[inds]
  y <- y.uniq[inds]
  xout <- x.uniq[sort(ceiling(runif(n.uniq*3, 0, n.uniq)))]
  w <- runif(n)
  m1 <- kernelSmooth(x, y, xout, weights = w, kernel = "triangular", bw = 1)
  m2 <- kernelSmooth(x, y, xout, weights = w, kernel = "triangular", bw = 1, no.dedup = TRUE)
  expect_equal(as.numeric(m1), as.numeric(m2), tolerance = 1e-5)
})

test_that("prediction on arbitrary grids works", {
  set.seed(1)
  x <- rnorm(50)
  y <- 1 + rt(50, df = 2)
  xgrid <- seq(-1.2, 1.2, 0.05)
  yhat0 <- kernelSmooth(x, y, xgrid, bw = 1, degree = 1, robust.iterations = 0, no.dedup = TRUE)
  yhat1 <- kernelSmooth(x, y, xgrid, bw = 1, degree = 1, robust.iterations = 1, no.dedup = TRUE)
  # plot(x, y); lines(xgrid, yhat0); lines(xgrid, yhat1, col = 2)
  expect_gt(var(yhat0), var(yhat1))
})

test_that("chunking for smoothing works", {
  set.seed(1)
  x <- rnorm(1000)
  y <- 1 + rt(1000, df = 2)
  xgrid <- seq(-1.5, 1.5, 0.05)
  yhat0 <- kernelSmooth(x, y, xgrid, bw = 0.5, no.dedup = TRUE, chunks = 1)
  yhat1 <- kernelSmooth(x, y, xgrid, bw = 0.5, no.dedup = TRUE, chunks = 2)
  # expect_identical(as.numeric(yhat0), as.numeric(yhat1))
  expect_equal(as.numeric(yhat0), as.numeric(yhat1), tolerance = 1e-12)

  # Now a slightly more difficult case -- LOO
  yhat2 <- kernelSmooth(x, y, bw = 0.5, no.dedup = TRUE, LOO = TRUE, chunks = 1)
  yhat3 <- kernelSmooth(x, y, bw = 0.5, no.dedup = TRUE, LOO = TRUE, chunks = 2)
  expect_equal(as.numeric(yhat2), as.numeric(yhat3), tolerance = 1e-12)
  # expect_identical fails on Windows...
})

