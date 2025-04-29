test_that("de-duplication in kernel preparation works correctly", {
  set.seed(1)
  n.uniq <- 1000
  n <- 6000
  inds <- ceiling(runif(n, 0, n.uniq))
  x.uniq <- matrix(rnorm(n.uniq*10), ncol = 10)
  x <- x.uniq[inds, ]  # x and y have many duplicates now
  y <- runif(n.uniq)[inds]
  xout <- x.uniq[ceiling(runif(n.uniq*3, 0, n.uniq)), ]
  w <- runif(n)
  a1 <- prepareKernel(x, y, xout, w, bw = 0.5)
  a2 <- prepareKernel(x, y, xout, w, bw = 0.5,
                       deduplicate.x = FALSE, deduplicate.xout = FALSE)
  a3 <- prepareKernel(x, y, xout, w, bw = 0.5, no.dedup = TRUE)
  expect_identical(a2, a3)
  expect_equal(sum(a1$weights), sum(a2$weights))
  expect_equal(unname(a2$duplicate.stats), rep(NA, 4))
  expect_type(a1$duplicate.stats, "double")
  expect_length(a1, 12)
  expect_length(a2, 12)
  expect_lt(object.size(a1), object.size(a2))
})
