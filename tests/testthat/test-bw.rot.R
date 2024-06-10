test_that("Silverman rule of thumb works", {
  x <- 1:32 / 4.69041575982343 # n = 32, n^(-1/5) = 1/2, sd(x) = 2
  expect_equal(bw.rot(x), 1.059223841)
})
