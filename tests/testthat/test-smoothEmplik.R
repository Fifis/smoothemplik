test_that("smoothEmplik works for a simple linear model", {
  set.seed(1)
  x <- sort(rlnorm(50))
  # Heteroskedastic DGP
  y <- 1 + 1*x + rnorm(50) * (1 + x + sin(x))
  mod.OLS <- lm(y ~ x)
  rho <- function(theta, ...) y - theta[1] - theta[2]*x  # Moment function
  w <- kernelWeights(x, PIT = TRUE, bw = 0.25)
  w <- w / rowSums(w)
  SEL <- function(b) smoothEmplik(rho = rho, theta = b, sel.weights = w)
  expect_type(SEL(coef(mod.OLS)), "double")
  gradSEL <- function(b) c(SEL(b + c(1e-5, 0)) - SEL(b - c(1e-5, 0)),
                           SEL(b + c(0, 1e-5)) - SEL(b - c(0, 1e-5)))/2e-5
  b.SEL <- optim(coef(mod.OLS), SEL, gr = gradSEL,
                 method = "BFGS", control = list(fnscale = -1, reltol = 1e-6))
  expect_equal(b.SEL$convergence, 0)
  expect_lt(b.SEL$value, 0)

  SELopt <- smoothEmplik(rho = rho, theta = b.SEL$par, sel.weights = w, type = "EL0", attach.attributes = "all")
  a <- attributes(SELopt)
  expect_equal(names(a), c("ELRs", "residuals", "lam", "nabla", "converged", "exitcode", "probabilities"))
  r <- rho(b.SEL$par)
  expect_equal(a$residuals, r)
  expect_equal(a$exitcode, rep(0, 50))
  expect_lt(max(abs(sapply(1:50, function(i) sum(w[i, ]*r / (1 + a$lam[i]*r))))), 1e-14)
})

test_that("Chunking of weight matrices works", {
  set.seed(1)
  x <- sort(rlnorm(50))
  y <- rlnorm(50) - 1
  w <- kernelWeights(pit(x), bw = 0.2, kernel = "epanechnikov")
  w <- w / rowSums(w)
  d <- data.frame(x = x, y = y)
  wfun <- function(ii, data) {
    px <- pit(data$x)
    k <- kernelWeights(x = px, xout = px[ii], bw = 0.2, kernel = "epanechnikov")
  }
  rho <- function(theta, ...) y - theta  # For mean estimation

  a1 <- smoothEmplik(rho, theta = 0, data = NULL, sel.weights = w, chunks = 1)
  a2 <- smoothEmplik(rho, theta = 0, data = d, sel.weights = wfun, chunks = 5)
  expect_equal(a1, a2)
})

test_that("The parabola calculator is working", {
  expect_equal(getParabola3(c(-1, 0, 2), c(-5, -4, 10)), c(2, 3, -4))
  expect_equal(getParabola(3, 23, 15, 4), c(2, 3, -4))
})

test_that("Sparse matrix to list conversion is working", {
  a <- rbind(rep(0,5), c(1,0,1,0,0), c(0,0.0001,0,0,0), c(0,1,0,1,1), c(1, 1, 1, 1, 1))
  al <- sparseMatrixToList(a)
  al2 <- sparseMatrixToList(a, trim = function(x) rep(0.01, length(x)))
  expect_length(al, 5)
  expect_length(al[[3]]$ct, 1)
  expect_length(al2[[3]]$ct, 0)
})

