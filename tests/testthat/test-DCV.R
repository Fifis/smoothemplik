test_that("DCV works for density cross-validation", {
  set.seed(1)
  n.uniq <- 100
  n <- 400
  inds <- sort(ceiling(runif(n, 0, n.uniq)))
  x.uniq <- sort(rnorm(n.uniq))
  x <- x.uniq[inds]
  w <- 1 + runif(n, 0, 2) # Relative importance
  bw.grid <- seq(0.1, 1.3, 0.2)
  # TODO: remove these messages!
  expect_true(all(is.finite(suppressMessages(DCV(x, bw = bw.grid, weights = w)))))
})
