 test_that("smoothEmplik works for a simple linear model", {
   set.seed(1)
   x <- sort(rlnorm(50))
   # Heteroskedastic DGP
   y <- 1 + 1*x + rnorm(50) * (1 + x + sin(x))
   mod.OLS <- lm(y ~ x)
   rho <- function(theta, ...) y - theta[1] - theta[2]*x  # Moment function
   w <- kernelWeights(x, PIT = TRUE, bw = 0.3)
   w <- w / rowSums(w)
   SEL <- function(b) smoothEmplik(rho = rho, theta = b, sel.weights = w)
   expect_type(SEL(coef(mod.OLS)), "double")
   b.SEL <- optim(coef(mod.OLS), SEL, method = "BFGS", control = list(fnscale = -1, reltol = 1e-6))
   expect_equal(b.SEL$convergence, 0)
   expect_lt(b.SEL$value, 0)

   SELopt <- smoothEmplik(rho = rho, theta = b.SEL$par, sel.weights = w, attach.attributes = "all")
   a <- attributes(SELopt)
   expect_equal(names(a), c("ELRs", "residuals", "lam", "nabla", "converged", "exitcode", "probabilities"))
   r <- rho(b.SEL$par)
   expect_equal(a$residuals, r)
   expect_equal(a$exitcode, rep(0, 50))
   expect_lt(max(abs(sapply(1:50, function(i) sum(w[i, ]*r / (1 + a$lam[i]*r))))), 1e-15)
 })
