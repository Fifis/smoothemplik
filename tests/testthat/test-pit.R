 test_that("the probability integral transform works", {
   x1 <- c(1, 1, 3, 4, 6)
   expect_equal(pit(x1), c(0.3, 0.3, 0.5, 0.7, 0.9))
   expect_equal(pit(x1, c(0, 2, 2, 5.9, 7, 8)), c(0.15, 0.4, 0.4, 0.8, 28/30, 29/30))
 })
