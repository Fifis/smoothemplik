test_that("chunking for density works", {
  set.seed(1)
  x <- rnorm(1000)
  xgrid <- seq(-1.5, 1.5, 0.05)
  fhat0 <- kernelDensity(x, xgrid, bw = 0.5, no.dedup = TRUE, chunks = 1)
  fhat1 <- kernelDensity(x, xgrid, bw = 0.5, no.dedup = TRUE, chunks = 2)
  # expect_identical(as.numeric(fhat0), as.numeric(fhat1))
  # expect_identical fails on Windows...
  expect_equal(as.numeric(fhat0), as.numeric(fhat1), tolerance = 1e-12)
})

test_that("grid is returned correctly", {
  set.seed(1)
  x <- rnorm(100)
  xgrid <- seq(-1.5, 1.5, 0.05)
  fhat0 <- kernelDensity(x, xgrid, bw = 0.5, no.dedup = TRUE, return.grid = FALSE)
  fhat1 <- kernelDensity(x, xgrid, bw = 0.5, no.dedup = TRUE, return.grid = TRUE)
  expect_identical(as.numeric(fhat0), fhat1[, 2])
  expect_identical(colnames(fhat1)[2], "density")
})

