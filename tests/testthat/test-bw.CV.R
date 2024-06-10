test_that("LSCV works for least-squares cross-validation", {
  set.seed(1)
  n.uniq <- 100
  n <- 400
  inds <- sort(ceiling(runif(n, 0, n.uniq)))
  x.uniq <- sort(rnorm(n.uniq))
  y.uniq <- 1 + 0.1*x.uniq + sin(x.uniq) + rnorm(n.uniq)
  x <- x.uniq[inds]
  y <- y.uniq[inds]
  w <- 1 + runif(n, 0, 2) # Relative importance
  bw.grid <- seq(0.1, 1.3, 0.2)
  CV2 <- LSCV(x, y = y, bw = bw.grid, weights = w)
  expect_true(all(is.finite(CV2)))
})
