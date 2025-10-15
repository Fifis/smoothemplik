test_that("Silverman rule of thumb handles bad inputs", {
  expect_error(bw.rot(1:9, kernel = "custom"), "Wrong kernel type")
  expect_error(bw.rot(c(1:9, NA)), "There are missing values")
  expect_warning(bw.rot(c(1:9, NA), na.rm = TRUE), "input data contain missing values")

  expect_error(bw.rot(matrix(1, nrow=1, ncol=2)), "Could not compute")
  expect_identical(bw.rot(matrix(1, nrow=2, ncol=2)), c(1, 1))

  x <- c(1:9, 100)
  expect_lt(bw.rot(x, robust = TRUE), bw.rot(x, robust = FALSE))
})

test_that("Silverman rule of thumb works", {
  x <- 1:32 / 4.69041575982343 # n = 32, n^(-1/5) = 1/2, sd(x) = 2
  expect_equal(bw.rot(x), 1.059223841048812, tolerance = 1e-15)
})
