 test_that("kernelSmooth linear and higher-order smoothers address boundary bias", {
   set.seed(1)
   x <- sort(runif(300, -6, 6))
   g <- seq(-6, 6, 0.1)
   f <- function(x) 1 + x + sin(x)
   y <- f(x) + rt(300, df = 4)

   m2lc <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", no.dedup = TRUE)
   m4lc <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", order = 4, no.dedup = TRUE)
   m2ll <- kernelSmooth(x, y, g, bw = 1.5, kernel = "triangular", degree = 1, no.dedup = TRUE)
   expect_gt(mean(abs(f(g) - m2lc)), mean(abs(f(g) - m4lc)))
   expect_gt(mean(abs(f(g) - m2lc)), mean(abs(f(g) - m2ll)))
 })

 test_that("de-duplication works in kernelSmooth", {
   set.seed(1)
   n.uniq <- 100
   n <- 1000
   inds <- sort(ceiling(runif(n, 0, n.uniq)))
   x.uniq <- sort(rnorm(n.uniq))
   y.uniq <- 1 + x.uniq + sin(x.uniq*2) + rnorm(n.uniq)
   x <- x.uniq[inds]
   y <- y.uniq[inds]
   xout <- x.uniq[sort(ceiling(runif(n.uniq*3, 0, n.uniq)))]
   w <- runif(n)
   m1 <- kernelSmooth(x, y, xout, w, kernel = "triangular", bw = 1)
   m2 <- kernelSmooth(x, y, xout, w, kernel = "triangular", bw = 1, no.dedup = TRUE)
   expect_equal(as.numeric(m1), as.numeric(m2), tolerance = 1e-5)
 })
