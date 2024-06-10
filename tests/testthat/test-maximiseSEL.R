 test_that("maximiseSEL works for a simple linear model", {
   set.seed(1)
   x <- sort(rlnorm(30))
   # Heteroskedastic DGP
   y <- 1 + 1*x + rnorm(30) * (1 + x + sin(x))
   mod.OLS <- lm(y ~ x)
   rho <- function(theta, ...) y - theta[1] - theta[2]*x  # Moment function
   w <- kernelWeights(x, PIT = TRUE, bw = 0.2, kernel = "triangular")
   w <- w / rowSums(w)

   mod.SEL  <- maximiseSEL(start.values = coef(mod.OLS), rho = rho, sel.weights = w)
   mod.SEL2 <- maximiseSEL(start.values = coef(mod.OLS), rho = rho, sel.weights = w, optmethod = "BFGS")
   expect_true(all(is.finite(mod.SEL$par)))
   expect_equal(mod.SEL$par, mod.SEL2$par, tolerance = 1e-5)
   expect_lt(mod.SEL$value, 0)
   expect_lte(mod.SEL$code, 2)
   expect_equal(mod.SEL2$code, 0)
 })
