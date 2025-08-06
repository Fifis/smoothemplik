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

test_that("LSCV handles multivariate inputs with unique and distinct bws, finding an improvement", {
  set.seed(1)
  x <- matrix(rlnorm(300), ncol = 3)
  y <- round(100 * x[, 1] + 200 * x[, 2] + 300 * x[, 3]) + rt(100, 1)
  x <- round(x, 1)

  system.time(expect_true(all(is.finite(suppressWarnings(
    b1 <- bw.CV(x = x, y = y, kernel = "epanechnikov", PIT = TRUE, tol = 1e-3, try.grid = FALSE, same = TRUE))))))
  system.time(expect_true(all(is.finite(suppressWarnings(
    b2 <- bw.CV(x = x, y = y, kernel = "epanechnikov", PIT = TRUE, tol = 1e-2, try.grid = FALSE, start.bw = rep(b1, 3)))))))
  expect_gt(LSCV(x = x, y = y, kernel = "epanechnikov", PIT = TRUE, bw = b1),
            LSCV(x = x, y = y, kernel = "epanechnikov", PIT = TRUE, bw = b2))
})

test_that("bw.CV de-duplicates correctly and minimises the CV criterion", {
  set.seed(1)
  n.uniq <- 50
  n <- 200
  inds <- sort(ceiling(runif(n, 0, n.uniq)))
  x.uniq <- sort(rnorm(n.uniq, sd = 2))
  y.uniq <- 1 + 0.1*x.uniq + sin(x.uniq) + rnorm(n.uniq)
  x <- x.uniq[inds]
  y <- y.uniq[inds]
  w <- 1 + runif(n, 0, 2)
  bw.grid <- seq(0.2, 1.3, 0.1)
  CV <- LSCV(x, y, bw.grid, weights = w)
  min.ind <- which.min(CV)
  bw.opt  <- bw.CV(x, y, w, attach.attributes = TRUE)
  bw.opt1 <- bw.CV(x, y, w, degree = 1, attach.attributes = TRUE)
  expect_gt(bw.opt1, bw.opt)

  # The optimal bandwidth must be close to the grid minimum
  expect_gt(bw.opt, bw.grid[min.ind-1])
  expect_lt(bw.opt, bw.grid[min.ind+1])
})
