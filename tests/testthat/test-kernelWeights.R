test_that("the rule-of-thumb bandwidth is invoked with a warning", {
  expect_warning(kernelWeights(1:10, PIT = TRUE), "No bandwidth supplied")
})
