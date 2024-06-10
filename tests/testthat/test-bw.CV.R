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

test_that("bw.CV de-duplicates correctly and minimises the CV criterion", {
  set.seed(1)
  n.uniq <- 100
  n <- 500
  inds <- sort(ceiling(runif(n, 0, n.uniq)))
  x.uniq <- sort(rnorm(n.uniq, sd = 2))
  y.uniq <- 1 + 0.1*x.uniq + sin(x.uniq) + rnorm(n.uniq)
  x <- x.uniq[inds]
  y <- y.uniq[inds]
  w <- 1 + runif(n, 0, 2)
  bw.grid <- seq(0.2, 1.3, 0.1)
  CV <- LSCV(x, y, bw.grid, weights = w)
  min.ind <- which.min(CV)
  bw.opt <- bw.CV(x, y, w)
  bw.opt2 <- bw.CV(x, y, w, no.dedup = TRUE)
  expect_equal(bw.opt, bw.opt2)

  # The optimal bandwidth must be close to the grid minimum
  expect_gt(bw.opt, bw.grid[min.ind-1])
  expect_lt(bw.opt, bw.grid[min.ind+1])
})
