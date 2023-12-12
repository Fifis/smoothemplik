# Codes for the empirical part of the paper, designs 1 and 2

#' Design functions for simulation
#'
#' The propensity score functions from Wang and Rao with possibly one of three specifications, as well as two special cases (always 1 or sigmoid).
#' @param x A numeric vector for which the propensity function is evaluated.
#' @param spec An integer indicating the propensity specification type, an integer from 1 to 5.
#' @param r A numeric scalar of the location parameter (the shift of the propensity or skedastic function).
#' @param lower A numeric scalar, the lower bound for the sigmoid (specification 5).
#' @param upper A numeric scalar, the upper bound for the sigmoid (specification 5).
#' @param sd A numeric scalar, the standard deviation of the normal distribution that yields the sigmoid for specification 5 (i.e. sigmoid smoothness).
#' @param v The power to which the absolute value of the error is raised in the skedastic function.
#' @param s The values of the propensity score function at \eqn{X = 0} and \eqn{X = 1} (discrete case).
#' @return A numeric vector of the same length as \code{x} with the propensity score function values (between 0 and 1).
#' @export
propensityScore <- function(x, spec = 5, r = 0.1, lower = 0.95, upper = 0.25, sd = 0.5) {
  a <- switch(spec,
    ifelse(abs(x - r) <= 1, 0.8 + 0.2 * abs(x - r), 0.95),
    ifelse(abs(x - r) <= 4, 0.9 - 0.2 * abs(x - r), 0.1),
    rep(0.6, length(x)),
    rep(1, length(x)),
    stats::pnorm(x, mean = r, sd = sd) * (upper - lower) + lower
  )
  if (is.null(a)) stop("propensityScore: Incorrect missing function specification number! Should be 1, 2, 3 or 4.") else return(a)
}

#' @export
#' @describeIn propensityScore The propensity score for the discrete case.
propensityScore.discrete <- function(x, s) s[x + 1]

#' @export
#' @describeIn propensityScore The skedastic function for the continuous case.
sigma2cragg <- function(x, r = -1 / 3, v = 2) abs(x - r)^v + 1 / 15

#' @export
#' @describeIn propensityScore The skedastic function for the discrete case.
sigma2.discrete <- function(x, s) s[x + 1]


#' Generation of data for Design 1 (continuous case)
#'
#' @param n Number of observations.
#' @param distfun Density function of X
#' @param eq.struc A function computing the deterministic part of the structural equation
#' @param eq.reduc A function computing the deterministic part of the reduced form
#' @param propensity Propensity score function
#' @param sigma2X Skedastic function (conditional variance of U given X)
#' @param varUV Covariance matrix of U and V
#' @param throwaway Points on the right end to discard
#' @param seed Seed for Monte-Carlo simulation
#'
#' @return A data frame with variables necessary for Design 1 simulations.
#' @export
#'
#' @examples
generateData <- function(n = 1000,
                         distfun = stats::runif,
                         eq.struc = function(z) 1 + 1 * z,
                         eq.reduc = function(x) 1 + 1 * x,
                         propensity = function(x) propensityScore(x),
                         sigma2X = function(x) sigma2cragg(x),
                         varUV = matrix(c(1, 1, 1, 2), ncol = 2),
                         throwaway = 0,
                         seed = 1
) {
  set.seed(seed)
  n <- n + throwaway
  X <- sort(distfun(n))
  UV <- matrix(stats::rnorm(n * 2), ncol = 2, nrow = n)
  e <- eigen(varUV, symmetric = TRUE)
  UV <- t(e$vectors %*% diag(sqrt(e$values)) %*% t(UV))
  U <- UV[, 1]
  VarUX <- sigma2X(X)
  Usked <- U * sqrt(VarUX)
  V <- UV[, 2]
  Z <- eq.reduc(X) + V # Reduced form
  Ystar <- eq.struc(Z) + Usked # Structural equation
  prop.score <- propensity(X)
  D <- stats::runif(n) < prop.score
  Y <- Ystar
  Y[!D] <- NA
  DY <- Y
  DY[!D] <- 0
  alldata <- data.frame(Z, X, D = as.numeric(D), pi = prop.score, Ystar, Y, DY, U, Usked, V, VarUX)
  return(alldata[1:(n - throwaway), ])
}


#' Generation of data for Design 2 (discrete case)
#'
#' @param n Number of observations.
#' @param p Probability of X=1 (Bernoulli)
#' @param params.struc Structural model parameter vector
#' @param params.reduc Reduced-form parameter vector
#' @param sigma2x Structural model error variance vector: Var(U | X) = x where x is 0 or 1
#' @param varU Structural model error variance multiplier
#' @param varV Reduced-form model error variance
#' @param covUV Covariance of the two errors, determining the degree of endogeneity
#' @param pix Propensity score evaluated at 0 and 1
#' @param seed Seed for Monte-Carlo simulation
#'
#' @return A data frame with variables necessary for Design 2 simulations.
#' @export
#'
#' @examples
generateDataDiscrete <- function(n = 1000,
                                 p = 0.3,
                                 params.struc = c(1, 1),
                                 params.reduc = c(0, 1),
                                 sigma2x = c(4, 1),
                                 varU = 1,
                                 varV = 2,
                                 covUV = 1,
                                 pix = c(0.3, 0.6),
                                 seed = 1
) {
  set.seed(seed)
  X <- stats::rbinom(n, 1, p)
  UV <- matrix(stats::rnorm(n * 2), ncol = 2, nrow = n)
  e <- eigen(matrix(c(varU, covUV, covUV, varV), nrow = 2, ncol = 2), symmetric = TRUE)
  UV <- t(e$vectors %*% diag(sqrt(e$values)) %*% t(UV))
  U <- UV[, 1]
  VarUX <- sigma2.discrete(x = X, s = sigma2x)
  Usked <- U * sqrt(VarUX)
  V <- UV[, 2]
  Z <- as.numeric(cbind(1, X) %*% params.reduc + V > 0)
  Ystar <- cbind(1, Z) %*% params.struc + Usked # Structural equation
  propensity <- pix[X + 1] # 1nd or 2nd element because X is 0 or 1
  D <- stats::runif(n) < propensity
  Y <- Ystar
  Y[!D] <- NA
  DY <- Y
  DY[!D] <- 0
  alldata <- data.frame(Z, X, D = as.numeric(D), Ystar, Y, DY, U, Usked, V, VarUX, pi = propensity)
  return(alldata)
}


#' Calculation of efficiency bounds
#'
#' This function is to be used if the distribution of the conditioning variable (X in our paper) is continuous.
#'
#' @param distfun The distribution function of the exogenous variable.
#' @param support The support over which \code{distfun} and its modifications should be integrated (passed to \code{integrate}).
#' @param eq.reduc The equation of the reduced form for the endogenous variable.
#' @param propensity The propensity score function.
#' @param sigma2X The skedastic function (conditional variance of errors).
#' @param varUV The covariance matrix of the baseline errors (U and V in our paper) of the structural and reduced equations.
#' @param rel.tol Desired relative accuracy of the numerical integration. Passed to \code{integrate}.
#'
#' @return A vector with 3 elements: the asymptotic variances of the inefficient slope estimator and the efficient slope estimator, and their ratio.
#' @export
#'
#' @examples
lb.design1.theor <- function(distfun = stats::dunif,
                             support = c(0, 1),
                             eq.reduc = function(x) 1 + 1 * x,
                             propensity = function(x) propensityScore(x, lower = 0.75, upper = 0.25, r = 0.4, sd = 0.05, spec = 5),
                             sigma2X = function(x) sigma2cragg(x, r = -1 / 3, v = 2),
                             varUV = matrix(c(1, 1, 1, 2), ncol = 2),
                             rel.tol = 1e-6) {
  denom <- function(x) sigma2X(x) * (1 / propensity(x) * (varUV[1, 1] - varUV[1, 2]^2 / varUV[2, 2]) + varUV[1, 2]^2 / varUV[2, 2])
  lb <- lbVS <- matrix(NA_real_, 2, 2)
  lb[1, 1] <- stats::integrate(function(x) 1 / denom(x) * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  lb[1, 2] <- lb[2, 1] <- stats::integrate(function(x) eq.reduc(x) / denom(x) * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  lb[2, 2] <- stats::integrate(function(x) eq.reduc(x)^2 / denom(x) * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  # VS
  mult <- function(x) propensity(x) / sigma2X(x)
  lbVS[1, 1] <- stats::integrate(function(x) mult(x) * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  lbVS[1, 2] <- lbVS[2, 1] <- stats::integrate(function(x) mult(x) * eq.reduc(x) * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  lbVS[2, 2] <- stats::integrate(function(x) mult(x) * eq.reduc(x)^2 * distfun(x), lower = support[1], upper = support[2], rel.tol = rel.tol)$value
  ineff <- solve(lbVS) * varUV[1, 1]
  eff <- solve(lb)
  return(c(VS = ineff[2, 2], eff = eff[2, 2], gains = ineff[2, 2] / eff[2, 2]))
}


#' Calculation of efficiency bounds
#'
#' This function is to be used if the distribution of the conditioning variable (X in our paper) has a two-point support and no numerical integration is needed.
#'
#' @param p The expected value of X (the probability of being equal to 1 in the Bernoulli distribution).
#' @param eq.reduc The equation of the reduced form for the endogenous variable.
#' @param propensity The propensity score vector (\eqn{\pi(0, 1)}).
#' @param sigma2X The skedastic function vector (\eqn{\sigma^2(0, 1)}).
#' @param varUV The covariance matrix of the baseline errors (U and V in our paper) of the structural and reduced equations.
#'
#' @return A vector with 3 elements: the asymptotic variances of the inefficient slope estimator and the efficient slope estimator, and their ratio.
#' @export
#'
#' @examples
lb.design2.theor <- function(p = 0.6,
                             eq.reduc = function(x) 0 + 1 * x,
                             propensity = c(0.25, 0.9),
                             sigma2X = c(16, 1),
                             varUV = matrix(c(1, 1, 1, 2), ncol = 2)) {
  XpZ.V <- function(x) eq.reduc(x) / sqrt(varUV[2, 2])
  G <- function(x) stats::dnorm(x) / (stats::pnorm(x) * stats::pnorm(-x))
  num <- function(x) {
    bp <- stats::pnorm(XpZ.V(x))
    return(c(1, bp, bp^2))
  }
  denom <- function(x) {
    frac1 <- ((1 - x) * sigma2X[1] + x * sigma2X[2]) / ((1 - x) * propensity[1] + x * propensity[2])
    exp1 <- varUV[1, 1] - varUV[1, 2]^2 / varUV[2, 2] * G(XpZ.V(x)) * stats::dnorm(XpZ.V(x))
    exp2 <- ((1 - x) * sigma2X[1] + x * sigma2X[2]) * varUV[1, 2]^2 / varUV[2, 2] * G(XpZ.V(x)) * stats::dnorm(XpZ.V(x))
    return(frac1 * exp1 + exp2)
  }
  Efrac <- num(0) / denom(0) * (1 - p) + num(1) / denom(1) * p
  lb <- lbVS <- matrix(NA_real_, 2, 2)
  lb[1, 1] <- Efrac[1]
  lb[1, 2] <- lb[2, 1] <- Efrac[2]
  lb[2, 2] <- Efrac[3]
  # VS
  mult <- function(x) ((1 - x) * propensity[1] + x * propensity[2]) / ((1 - x) * sigma2X[1] + x * sigma2X[2])
  Efrac <- num(0) * mult(0) * (1 - p) + num(1) * mult(1) * p
  lbVS[1, 1] <- Efrac[1]
  lbVS[1, 2] <- lbVS[2, 1] <- Efrac[2]
  lbVS[2, 2] <- Efrac[3]
  #
  ineff <- solve(lbVS) * varUV[1, 1]
  eff <- solve(lb)
  return(c(VS = ineff[2, 2], eff = eff[2, 2], gains = ineff[2, 2] / eff[2, 2]))
}


#' Estimate a Design-2 model with Owen's EL
#'
#' @param seed The seed for data generation (passed to \code{generateDataDiscrete}).
#' @param design A list with design options (can have the following named elements: number of observations, structural parameters, reduced-form parameters, skedastic function values, propensity score values).
#' @param do.SE Logical: compute EL-based standard errors by inverting the numerically estimated Hessian and include results in the output?
#' @param do.restr Logical: Test the true null hypothesis that the coefficients are equal to their true values (one by one and jointly) via the ELR?
#' @param do.CI Logical: Compute the confidence regions for 90%, 95%, 97.5% and 99% levels?
#'
#' Since the slope confidence regions are one-dimensional, unimodality of the ELR function is assumed, so the regions become intervals, possibly unbounded, and are found numerically)
#'
#' @return A list with the unrestricted estimates, constrained estimates and ELR test statistics (if \code{do.restr}),
#' EL confidence intervals (if \code{do.CI}), and EL-based standard errors (if \code{do.SE}).
#' @export
#'
#' @examples
getCoefELDiscrete <- function(seed, design = list(n = 500, p = 0.6, params.struc = c(1, 1), params.reduc = c(0, 1), sigma2x = c(16, 1), pix = c(0.25, 0.9)),
                              do.SE = TRUE, do.restr = TRUE, do.CI = TRUE) {
  data <- generateDataDiscrete(n = design$n, p = design$p, params.reduc = design$params.reduc, sigma2x = design$sigma2x, pix = design$pix, seed = seed)
  data.VS <- data[as.logical(data$D), ]
  tic0 <- Sys.time()
  mod1sobs <- stats::lm(Z ~ X, data = data.VS)
  weak.VS <- summary(mod1sobs)$coefficients[2, 4] >= 0.05
  # We want the coefficient on X to be significant at least at 5% level
  # This gets used in the efficient estimation part
  Zhat <- mod1sobs$fitted.values
  m <- stats::lm(data.VS$Y ~ Zhat)
  mod2sobs <- m$coefficients

  # Likelihood-based variance
  if (do.SE) {
    margs <- list(eps = 1e-4, d = 1e-2)
    vcov.vs <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.vs(x, data = data))$logelr, x = mod2sobs, method.args = margs))
    se.lik <- unname(sqrt(diag(vcov.vs)))
    if (!all(is.finite(se.lik))) {
      margs <- list(eps = 1e-4, d = 1e-3)
      vcov.vs <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.vs(x, data = data))$logelr, x = mod2sobs, method.args = margs))
      se.lik <- unname(sqrt(diag(vcov.vs)))
      warning("Using a smaller initial step (d=0.001) for Hessian computation.")
    }
    if (!all(is.finite(se.lik))) {
      margs <- list(eps = 1e-4, d = 1e-4)
      vcov.vs <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.vs(x, data = data))$logelr, x = mod2sobs, method.args = margs))
      se.lik <- unname(sqrt(diag(vcov.vs)))
      warning("Using a MUCH smaller initial step (d=0.0001) for Hessian computation.")
    }
  } else {
    se.lik <- NULL
  }
  rm(mod1sobs, m)
  # Inefficient complete-case SEL with RoT bandwidth
  el.complete <- mod2sobs
  LR.complete <- cemplik(g.unconditional.vs(c(1, 1), data = data))$logelr

  if (do.restr) {
  frestr <- function(theta0 = 1, theta1 = 1) {
    r <- -cemplik(g.unconditional.vs(c(theta0, theta1), data = data))$logelr
    if (r > 100) return(NA) else return(r)
  }
  proflik1 <- function(theta1) {
    if (abs(theta1) > 1e6) return(NA)
    ret <- tryCatch(stats::nlm(f = function(x) frestr(theta0 = x, theta1 = theta1), p = mean(data.VS$Y) - theta1 * mean(data.VS$Z), gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik[1] / 5), error = function(e) return(list(code = 5)))
    if (!(ret$code %in% c(4, 5))) {
      return(ret$minimum)
    } else {
      ret <- tryCatch(stats::optim(fn = function(x) frestr(theta0 = x, theta1 = theta1), par = mean(data.VS$Y) - theta1 * mean(data.VS$Z), gr = NULL, method = "BFGS", control = list(reltol = 1e-8))$value, error = function(e) return(NA))
      return(ret)
    }
  }
  # Getting the semiparametrically efficient estimator under constraints
  sigmax <- sqrt(data.VS$VarUX)
  Yp <- kernelDiscreteDensitySmooth(x = data.VS$X, y = data.VS$Y)$y / sigmax
  Cp <- 1 / sigmax
  Zp <- kernelDiscreteDensitySmooth(x = data.VS$X, y = data.VS$Z)$y / sigmax
  theta0c <- (mean(Yp * Cp) - 1 * mean(Cp * Zp)) / mean(Cp^2)
  theta1c <- (mean(Yp * Zp) - 1 * mean(Cp * Zp)) / mean(Zp^2)
  restr0 <- tryCatch(stats::nlm(f = function(x) frestr(theta0 = x, theta1 = 1), p = theta0c, gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik[1] / 5), error = function(e) return(list(code = 5)))
  if (restr0$code %in% c(4, 5)) {
    restr0 <- stats::optim(fn = function(x) frestr(theta0 = x, theta1 = 1), par = theta0c, method = "BFGS", gr = NULL, control = list(reltol = 1e-8))
    restr0 <- list(minimum = restr0$value, estimate = restr0$par, gradient = NA, code = restr0$convergence, iterations = unname(restr0$counts[2]))
  }
  restr1 <- tryCatch(stats::nlm(f = function(x) frestr(theta0 = 1, theta1 = x), p = theta1c, gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik[2] / 5), error = function(e) return(list(code = 5)))
  if (restr1$code %in% c(4, 5)) {
    restr1 <- stats::optim(fn = function(x) frestr(theta0 = 1, theta1 = x), par = theta1c, method = "BFGS", gr = NULL, control = list(reltol = 1e-8))
    restr1 <- list(minimum = restr1$value, estimate = restr1$par, gradient = NA, code = restr1$convergence, iterations = unname(restr1$counts[2]))
  }
  } else {
    restr0 <- restr1 <- NULL
  }

  cat("Seed ", .lead0(seed, 4), ", ineff. estim. with I(Xi = Xj) (1/6), theta1 = ", sprintf("%1.3f", el.complete[2]), "\n")
  if (do.restr) cat("Seed ", .lead0(seed, 4), ", ineff. estim. with theta1 = 1 (2/6), theta0 = ", sprintf("%1.3f", restr0$estimate), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(restr0$minimum, 1)), "\n",
      "Seed ", .lead0(seed, 4), ", ineff. estim. with theta0 = 1 (3/6), theta1 = ", sprintf("%1.3f", restr1$estimate), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(restr1$minimum, 1)), "\n", sep = "")
  tic1 <- Sys.time()

  # Empirical confidence intervals
  if (do.CI) {
  ps <- c(0.1, 0.05, 0.025, 0.01)
  qs <- stats::qchisq(1 - ps, df = 1) / 2
  start.mult <- stats::qnorm(1 - ps/2)
  rCI <- lCI <- vector("list", length(ps))
  rCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1(p),
                                      lower = el.complete[2] + 0.95*start.mult[1]*se.lik[2], upper = el.complete[2] + 1.1*start.mult[1]*se.lik[2],
                                      extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
  if(is.na(rCI[[1]]$root)) rCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1(p),
                                                               lower = el.complete[2] + 0.1*start.mult[1]*se.lik[2], upper = el.complete[2] + 0.11*start.mult[1]*se.lik[2],
                                                               extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
  lCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1(p),
                                      lower = el.complete[2] - 1.1*start.mult[1]*se.lik[2], upper = el.complete[2] - 0.95*start.mult[1]*se.lik[2],
                                      extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  if(is.na(lCI[[1]]$root)) lCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1(p),
                                                               lower = el.complete[2] - 0.11*start.mult[1]*se.lik[2], upper = el.complete[2] - 0.10*start.mult[1]*se.lik[2],
                                                               extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  for (i in 2:4) {
    rCI[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - proflik1(p), lower = rCI[[i-1]]$root, f.lower = rCI[[i-1]]$f.root - qs[i-1] + qs[i],
                                        upper = el.complete[2] + (rCI[[i-1]]$root - el.complete[2]) / start.mult[i-1] * start.mult[i],
                                        extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
    lCI[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - proflik1(p), upper = lCI[[i-1]]$root, f.upper = lCI[[i-1]]$f.root - qs[i-1] + qs[i],
                                        lower = el.complete[2] + (lCI[[i-1]]$root - el.complete[2]) / start.mult[i-1] * start.mult[i],
                                        extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  }
  rCI <- unlist(lapply(rCI, "[[", "root"))
  lCI <- unlist(lapply(lCI, "[[", "root"))
  }

  tic2 <- Sys.time()

  # Efficient EL
  pihat <- kernelDiscreteDensitySmooth(x = data$X, y = data$D)$y
  pihat[pihat < 1 / nrow(data)] <- 1 / nrow(data)
  pihat[pihat > 1 - 1 / nrow(data)] <- 1 - 1 / nrow(data)
  data$XZ <- data$X * 2 + data$Z
  ystarhat <- kernelDiscreteDensitySmooth(x = data$XZ, y = data$DY)$y / pihat

  # If there are weak instruments in the VS, then using el.complete as starting values is dangerous!
  # Getting a better first-stage projection is key, then
  start.efficient <- if (weak.VS) unname(stats::lm(data$Y ~ stats::lm(Z ~ X, data = data)$fitted.values)$coefficients) else unname(el.complete)
  nlm.step.max <- sqrt(sum(diag(vcov.vs)))
  el.opt <- tryCatch(stats::nlm(function(theta) -cemplik(g.unconditional.eff(theta, data = data, pihat = pihat, ystarhat = ystarhat))$logelr,
    p = start.efficient, gradtol = 1e-8, steptol = 1e-8, stepmax = nlm.step.max), error = function(e) return(list(code = 5)))
  if (el.opt$code %in% c(4, 5)) {
    el.opt <- stats::optim(fn = function(theta) -cemplik(g.unconditional.eff(theta, data = data, pihat = pihat, ystarhat = ystarhat))$logelr,
      par = start.efficient, method = "BFGS", gr = NULL, control = list(reltol = 1e-8))
    el.opt <- list(minimum = el.opt$value, estimate = el.opt$par, gradient = NA, code = el.opt$convergence, iterations = unname(el.opt$counts[2]))
  }
  el.efficient <- el.opt$estimate
  LR.efficient <- cemplik(g.unconditional.eff(c(1, 1), data = data, pihat = pihat, ystarhat = ystarhat))$logelr
  if (do.SE) {
    margs <- list(eps = 1e-4, d = 1e-2)
    vcov.eff <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.eff(x, data = data, pihat = pihat, ystarhat = ystarhat))$logelr, el.efficient, method.args = margs))
    se.lik.eff <- sqrt(diag(vcov.eff))
    if (!all(is.finite(se.lik.eff))) {
      margs <- list(eps = 1e-4, d = 1e-3)
      vcov.eff <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.eff(x, data = data, pihat = pihat, ystarhat = ystarhat))$logelr, el.efficient, method.args = margs))
      se.lik.eff <- sqrt(diag(vcov.eff))
      warning("Using a smaller initial step (d=0.001) for Hessian computation.")
    }
    if (!all(is.finite(se.lik.eff))) {
      margs <- list(eps = 1e-4, d = 1e-4)
      vcov.eff <- solve(-numDeriv::hessian(function(x) cemplik(g.unconditional.eff(x, data = data, pihat = pihat, ystarhat = ystarhat))$logelr, el.efficient, method.args = margs))
      se.lik.eff <- sqrt(diag(vcov.eff))
      warning("Using a MUCH smaller initial step (d=0.0001) for Hessian computation.")
    }
  } else {
    se.lik.eff <- NULL
  }
  if (do.restr) {
    frestr.eff <- function(theta0 = 1, theta1 = 1) {
    r <- -cemplik(g.unconditional.eff(c(theta0, theta1), data = data, pihat = pihat, ystarhat = ystarhat))$logelr
    if (r > 100) return(NA) else return(r)
  }
  restr0.eff <- tryCatch(stats::nlm(f = function(x) frestr.eff(theta0 = x, theta1 = 1), p = theta0c, gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik.eff[1] / 5), error = function(e) return(list(code = 5)))
  if (restr0.eff$code %in% c(4, 5)) {
    restr0.eff <- stats::optim(fn = function(x) frestr.eff(theta0 = x, theta1 = 1), par = theta0c, method = "BFGS", gr = NULL, control = list(reltol = 1e-8))
    restr0.eff <- list(minimum = restr0.eff$value, estimate = restr0.eff$par, gradient = NA, code = restr0.eff$convergence, iterations = unname(restr0.eff$counts[2]))
  }
  restr1.eff <- tryCatch(stats::nlm(f = function(x) frestr.eff(theta0 = 1, theta1 = x), p = theta1c, gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik.eff[2] / 5), error = function(e) return(list(code = 5)))
  if (restr1.eff$code %in% c(4, 5)) {
    restr1.eff <- stats::optim(fn = function(x) frestr.eff(theta0 = 1, theta1 = x), par = theta1c, method = "BFGS", gr = NULL, control = list(reltol = 1e-8))
    restr1.eff <- list(minimum = restr1.eff$value, estimate = restr1.eff$par, gradient = NA, code = restr1.eff$convergence, iterations = unname(restr1.eff$counts[2]))
  }
  } else {
    restr0.eff <- restr1.eff <- NULL
  }
  cat("Seed ", .lead0(seed, 4), ", effic. estim. with I(Xi = Xj) (4/6), theta1 = ", sprintf("%1.3f", el.efficient[2]), "\n")
  if (do.restr) cat("Seed ", .lead0(seed, 4), ", effic. estim. with theta1 = 1 (5/6), theta0 = ", sprintf("%1.3f", restr0.eff$estimate), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(restr0.eff$minimum, 1)), "\n",
      "Seed ", .lead0(seed, 4), ", effic. estim. with theta0 = 1 (6/6), theta1 = ", sprintf("%1.3f", restr1.eff$estimate), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(restr1.eff$minimum, 1)), "\n", sep = "")
  tic3 <- Sys.time()

  if (do.CI) {
  proflik1.eff <- function(theta1) {
    if (abs(theta1) > 1e6) return(NA)
    p0 <- mean(data.VS$Y) - theta1 * mean(data.VS$Z)
    ret <- tryCatch(stats::nlm(f = function(x) frestr.eff(theta0 = x, theta1 = theta1), p = p0, gradtol = 1e-8, steptol = 1e-8, stepmax = se.lik.eff[1] / 5), error = function(e) return(list(code = 5)))
    if (!(ret$code %in% c(4, 5))) {
      return(ret$minimum)
    } else {
      ret <- tryCatch(stats::optim(fn = function(x) frestr.eff(theta0 = x, theta1 = theta1), par = p0, gr = NULL, method = "BFGS", control = list(reltol = 1e-8))$value, error = function(e) return(NA))
      return(ret)
    }
  }
  rCI.eff <- lCI.eff <- vector("list", length(ps))
  rCI.eff[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1.eff(p),
                                          lower = el.efficient[2] + 0.95*start.mult[1]*se.lik.eff[2], upper = el.efficient[2] + 1.1*start.mult[1]*se.lik.eff[2],
                                          extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
  if (is.na(rCI.eff[[1]]$root)) rCI.eff[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1.eff(p),
                                                                        lower = el.efficient[2] + 0.10*start.mult[1]*se.lik.eff[2], upper = el.efficient[2] + 0.11*start.mult[1]*se.lik.eff[2],
                                                                        extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
  lCI.eff[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1.eff(p),
                                          lower = el.efficient[2] - 1.1*start.mult[1]*se.lik.eff[2], upper = el.efficient[2] - 0.95*start.mult[1]*se.lik.eff[2],
                                          extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  if (is.na(rCI.eff[[1]]$root)) lCI.eff[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - proflik1.eff(p),
                                                                        lower = el.efficient[2] - 0.11*start.mult[1]*se.lik.eff[2], upper = el.efficient[2] - 0.10*start.mult[1]*se.lik.eff[2],
                                                                        extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  for (i in 2:4) {
    rCI.eff[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - proflik1.eff(p), lower = rCI.eff[[i-1]]$root, f.lower = rCI.eff[[i-1]]$f.root - qs[i-1] + qs[i],
                                            upper = el.efficient[2] + (rCI.eff[[i-1]]$root - el.efficient[2]) / start.mult[i-1] * start.mult[i],
                                            extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
    lCI.eff[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - proflik1.eff(p), upper = lCI.eff[[i-1]]$root, f.upper = lCI.eff[[i-1]]$f.root - qs[i-1] + qs[i],
                                            lower = el.efficient[2] + (lCI.eff[[i-1]]$root - el.efficient[2]) / start.mult[i-1] * start.mult[i],
                                            extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
  }
  rCI.eff <- unlist(lapply(rCI.eff, "[[", "root"))
  lCI.eff <- unlist(lapply(lCI.eff, "[[", "root"))
  }
  tic4 <- Sys.time()

  names(el.complete) <- names(el.efficient) <- c("const", "Z")
  if (do.CI) {
    CI.complete <- matrix(c(lCI, rCI), ncol = 2)
    CI.efficient <- matrix(c(lCI.eff, rCI.eff), ncol = 2)
    row.names(CI.complete) <- row.names(CI.efficient) <- paste0("q", c(90, 95, 97.5, 99))
  } else {
    CI.complete <- CI.efficient <- NULL
  }

  times <- as.numeric(c(difftime(tic1, tic0, units = "secs"), difftime(tic2, tic1, units = "secs"), difftime(tic3, tic2, units = "secs"), difftime(tic4, tic3, units = "secs")))

  return(list(
    el.complete = el.complete,
    el.complete.restr0 = c(const = restr0$estimate, Z = 1),
    el.complete.restr1 = c(const = 1, Z = restr1$estimate),
    LR.complete = c(bothEQ1 = LR.complete, slopeEQ1 = if (do.restr) -restr0$minimum else NULL, interceptEQ1 = if (do.restr) -restr1$minimum else NULL),
    CI.complete = CI.complete,
    se.complete = se.lik,
    el.efficient = el.efficient,
    el.efficient.restr0 = c(const = restr0.eff$estimate, Z = 1),
    el.efficient.restr1 = c(const = 1, Z = restr1.eff$estimate),
    LR.efficient = c(bothEQ1 = LR.efficient, slopeEQ1 = if (do.restr) -restr0.eff$minimum else NULL, interceptEQ1 = if (do.restr) -restr1.eff$minimum else NULL),
    CI.efficient = CI.efficient,
    se.efficient = se.lik.eff,
    seconds = times
  ))
}

#' Estimation of one model via SEL
#'
#' @param start.mod A list that must contain the following elements: \code{coefficients}
#'   (initial values), \code{vcov} (a plausible guess about the variance-covariance of the estimator;
#'   used when teremining the optimiser step), and \code{residuals} (true residuals
#'   for the validation sample, i.e. subset: data$D == 1)
#' @param data A data frame on which the linear model is estimated.
#' @param sel.bw SEL smoothing bandwidth.
#' @param eff.control NULL for the VS estimator, or a list passed to \code{rho.full.sample}
#'   that may contain any of the following: \code{pi.hat}, \code{helper},
#'   \code{helper.predicted}, \code{pi.bw}, \code{helper.bw}, \code{helper.degree}, \code{PIT}.
#' @param do.SE Logical: should the standard errors be computed numerically via
#'   SEL Hessian inversion? Will be overriden with \code{TRUE} if \code{CI.lev}
#'   is not \code{NULL} since the intervals rely on standard errors.
#' @param do.restr # Estimate models under constraints where one of the coefficients is equal to the true value?
#' @param CI.lev A numberic vector of probabilities for confidence intervale;
#'   a good value is c(0.9, 0.95, 0.99). If \code{NULL}, the confidence intervals are not returned.
#' @param try.ellipse Logical: should the optimiser safeguard against potential
#'   local optima (in very small samples) in some likely places based on the VCOV matrix?
#' @param ... Passed to maximiseSEL
#'
#' @return A list with model type, SEL estimates, SEL-based variance-covariance
#'   matrix (if \code{do.SE}), constrained estimates and ELR test statistics
#'   (if \code{do.restr}), ELR-based confidence intervals (if CI.lev is not NULL), and smoothing bandwidths.
#' @export
#'
#' @examples
estimateOneModelDesign1 <- function(start.mod,
                             data,
                             sel.bw,
                             eff.control = list(),
                             do.SE = FALSE,
                             do.restr = FALSE,
                             CI.lev = NULL, #
                             try.ellipse = TRUE, # Initial value search
                             ... #
) {
  if (!is.null(CI.lev)) do.SE <- TRUE
  data.VS <- data[as.logical(data$D), ]
  eff.controls <- list(pi.hat = NULL, pi.bw = bw.rot(data$X), helper.bw = bw.rot(data$X),
                       helper = "gstar", helper.predicted = NULL, helper.degree = 1, PIT = TRUE)
  if (is.list(eff.control) && length(eff.control) > 0) {
    eff.controls[names(eff.control)] <- eff.control
    rho <- function(thetarho) rho.full.sample(theta = thetarho, data = data, pi.hat = eff.controls$pi.hat,
                                              helper = eff.controls$helper, helper.predicted = eff.controls$helper.predicted,
                                              pi.bw = eff.controls$pi.bw, helper.bw = eff.controls$helper.bw,
                                              helper.degree = eff.controls$helper.degree, PIT = eff.controls$PIT)
    model.type <- "full-sample-efficient"
  } else {
    rho <- function(thetarho) rho.complete.case(theta = thetarho, data = data)
    model.type <- "validation-sample-inefficient"
  }
  sel.weights <- getSELWeights(data$X, bw = sel.bw)
  if (try.ellipse && !is.null(dim(start.mod$coefficients))) { # Adding candidate points for optimisation
    extra.points <- sampleEllipse(centre = start.mod$coefficients[1, ], vcov = start.mod$vcov, radius = sqrt(stats::qchisq(c(0.75, 0.9, 0.99), df = 2)))
    start.mod$coefficients <- rbind(start.mod$coefficients, do.call(rbind, extra.points))
  }
  start.point <- start.mod$coefficients
  if (!is.null(dim(start.mod$coefficients)[1])) { # If there are multiple candidate values
    vals <- apply(start.mod$coefficients, 1, function(p) tryCatch(maximiseSEL(rho = rho, restricted.params = p, sel.weights = sel.weights, ...),  error = .fail))
    vals <- unlist(lapply(vals, "[[", "value"))
    best.point <- which.max(vals)
    start.point <- start.mod$coefficients[best.point, ]
  }
  sel <- tryCatch(maximiseSEL(rho = rho, start.values = start.point, sel.weights = sel.weights, ...),  error = .fail)
  restr.both  <- tryCatch(maximiseSEL(rho = rho, restricted.params = c(1, 1), sel.weights = sel.weights, ...),  error = .fail)

  vcovar <- se <- NULL
  if (do.SE) {
    margs <- list(eps = 1e-4, d = 1e-2)
    vcovar <- solve(-numDeriv::hessian(function(th) smoothEmplik(rho, sel.weights, thetarho = th), x = sel$par, method.args = margs))
    se <- sqrt(diag(vcovar))
    if (!all(is.finite(se))) {
      margs <- list(eps = 1e-4, d = 1e-3)
      vcovar <- solve(-numDeriv::hessian(function(th) smoothEmplik(rho, sel.weights, thetarho = th), x = sel$par, method.args = margs))
      se <- sqrt(diag(vcovar))
      warning("Using a smaller initial step (d=0.001) for Hessian computation.")
    }
    if (!all(is.finite(se))) {
      margs <- list(eps = 1e-4, d = 1e-4)
      vcovar <- solve(-numDeriv::hessian(function(th) smoothEmplik(rho, sel.weights, thetarho = th), x = sel$par, method.args = margs))
      se <- sqrt(diag(vcovar))
      warning("Using a MUCH smaller initial step (d=0.0001) for Hessian computation.")
    }
    colnames(vcovar) <- rownames(vcovar) <- names(se) <- names(start.point)
  }

  stepmax0 <- if (do.SE) se[1] / 3 else sqrt(diag(start.mod$vcov)[1]) / 3
  stepmax1 <- if (do.SE) se[2] / 3 else sqrt(diag(start.mod$vcov)[2]) / 3

  if (do.restr || !is.null(CI.lev)) {
    # Getting the semiparametrically efficient estimator under constraints
    U2 <- start.mod$residuals^2
    VarUX <- kernelSmooth(data.VS$X, U2, bw = bw.CV(x = data.VS$X, y = U2))
    sigmax <- sqrt(VarUX)
    Yp <- kernelSmooth(x = data.VS$X, y = data.VS$Y, bw = bw.rot(data.VS$X), degree = 1) / sigmax
    Cp <- 1 / sigmax
    Zp <- kernelSmooth(x = data.VS$X, y = data.VS$Z, bw = bw.rot(data.VS$X), degree = 1) / sigmax
  }

  if (do.restr) {
    theta0c <- (mean(Yp * Cp) - 1 * mean(Cp * Zp)) / mean(Cp^2)
    theta1c <- (mean(Yp * Zp) - 1 * mean(Cp * Zp)) / mean(Zp^2)
    restr0 <- maximiseSEL(rho = rho, start.values = theta0c, restricted.params = c(NA, 1),
                          sel.weights = sel.weights, optmethod = "nlm", nlm.step.max = stepmax0, ...)
    restr1 <- maximiseSEL(rho = rho, start.values = theta1c, restricted.params = c(1, NA),
                          sel.weights = sel.weights, optmethod = "nlm", nlm.step.max = stepmax1, ...)
    LR0 <- 2 * (sel$value - restr0$value)
    LR1 <- 2 * (sel$value - restr1$value)
  } else {
    restr0 <- restr1 <- NULL
    LR0 <- LR1 <- NA
  }
  LR.both <- 2 * (sel$value - restr.both$value)

  if (!is.null(CI.lev)) {
    LR2 <- function(theta1) {
      theta0r <- (mean(Yp * Cp) - theta1 * mean(Cp * Zp)) / mean(Cp^2)
      2*(sel$value - maximiseSEL(rho = rho, start.values = theta0r, restricted.params = c(NA, theta1), sel.weights = sel.weights, optmethod = "nlm", nlm.step.max = stepmax0, ...)$value)
    }

    ps <- 1 - CI.lev
    qs <- stats::qchisq(1 - ps, df = 1)
    start.mult <- stats::qnorm(1 - ps/2)
    rCI <- lCI <- vector("list", length(ps))
    rCI[[1]] <-                          tryCatch(stats::uniroot(f = function(p) qs[1] - LR2(p),
                                                                 lower = sel$par[2] + 0.95*start.mult[1]*se[2], upper = sel$par[2] + 1.10*start.mult[1]*se[2],
                                                                 extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
    if(is.na(rCI[[1]]$root)) rCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - LR2(p),
                                                                 lower = sel$par[2] + 0.10*start.mult[1]*se[2], upper = sel$par[2] + 0.11*start.mult[1]*se[2],
                                                                 extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
    lCI[[1]] <-                          tryCatch(stats::uniroot(f = function(p) qs[1] - LR2(p),
                                                                 lower = sel$par[2] - 1.10*start.mult[1]*se[2], upper = sel$par[2] - 0.95*start.mult[1]*se[2],
                                                                 extendInt = "upX", tol = 1e-8),   error = function(e) return(list(root = NA)))
    if(is.na(lCI[[1]]$root)) lCI[[1]] <- tryCatch(stats::uniroot(f = function(p) qs[1] - LR2(p),
                                                                 lower = sel$par[2] - 0.11*start.mult[1]*se[2], upper = sel$par[2] - 0.10*start.mult[1]*se[2],
                                                                 extendInt = "upX", tol = 1e-8),   error = function(e) return(list(root = NA)))
    if (length(CI.lev) > 1) {
      for (i in 2:length(CI.lev)) {
        rCI[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - LR2(p),
                                            lower = rCI[[i-1]]$root, f.lower = rCI[[i-1]]$f.root - qs[i-1] + qs[i],
                                            upper = sel$par[2] + (rCI[[i-1]]$root - sel$par[2]) / start.mult[i-1] * start.mult[i],
                                            extendInt = "downX", tol = 1e-8), error = function(e) return(list(root = NA)))
        lCI[[i]] <- tryCatch(stats::uniroot(f = function(p) qs[i] - LR2(p),
                                            upper = lCI[[i-1]]$root, f.upper = lCI[[i-1]]$f.root - qs[i-1] + qs[i],
                                            lower = sel$par[2] + (lCI[[i-1]]$root - sel$par[2]) / start.mult[i-1] * start.mult[i],
                                            extendInt = "upX", tol = 1e-8), error = function(e) return(list(root = NA)))
      }
    }
    rCI <- unlist(lapply(rCI, "[[", "root"))
    lCI <- unlist(lapply(lCI, "[[", "root"))
    CI <- cbind(lCI, rCI)
    rownames(CI) <- as.character(CI.lev)
  } else {
    CI <- NULL
  }

  return(list(model.type = model.type, sel = sel, vcov = vcovar,
              restricted = list(`R:theta=(1,1)` = restr.both, `R:theta1=1` = restr0, `R:theta0=1` = restr1),
              LR = c(`R:theta=(1,1)` = LR.both, `R:theta1=1` = LR0, `R:theta0=1` = LR1),
              CI = CI, eff.control = if (length(eff.control) > 0) eff.controls else NULL, sel.bw = sel.bw))
}


#' Simulate and estimate one model with one seed and many hyper-parameters
#'
#' @param seed One numeric seed for Monte-Carlo simulations.
#' @param design A list containing the number of observations, skedastic function, and propensity score.
#' @param sel.bw A numeric vector of SEL smoothing bandwidths. Mandatory.
#' @param pi.bw NULL or a numeric vector of smoothing bandwidths for the propensity.
#' @param helper.bw NULL or a numeric vector of smoothing bandwidths for mu.
#' @param helper Type of intermediate object used in the prediction of mu.
#' @param degree Degree of non-parametric smoothing of the helper.
#' @param CV Logical: pick the smoothing bandwidth for pi(X) and the helper via cross-validation with PIT?
#' @param weak.instrument.f A scalar indicating the threshold F statistic: if the first-stage F is less than that, discard the simulation. Set to 0 to allow all simulations regardless of the instrument strength.
#' @param ... Passed to estimateOneModelDesign1.
#'
#' @return A list with seed used, design parameters to recreate the data set,
#'   2SLS estimates, weak instrument diagnostics, optimal IV estimates, VS SEL
#'   estimates and, if proper smoothing parameters were passed, efficient SEL
#'   estimates, and time taken to compute the latter two.
#' @export
#'
#' @examples
getCoefSELContinuous <- function(seed, design = list(),
                                 sel.bw,
                                 pi.bw = NULL, helper.bw = NULL, helper = NULL,
                                 degree = NULL,
                                 CV = FALSE,
                                 weak.instrument.f = 10,
                                 ...
) {
  this.design <- list(n = 500, sigma2X = function(x) sigma2cragg(x, r = -1/3, v = 2), propensity = function(x) propensityScore(x, lower = 0.95, upper = 0.25, r = 0.1, sd = 0.5, spec = 5))
  if (is.list(design) && length(design) > 0) this.design[names(design)] <- design
  data <- generateData(n = this.design$n, distfun = stats::runif, propensity = this.design$propensity, sigma2X = this.design$sigma2X, seed = seed)
  data.VS <- data[as.logical(data$D), ]
  tic0 <- Sys.time()
  mod1s.VS <- stats::lm(Z ~ X, data = data.VS)
  mod2s.VS <- unname(stats::lm(data.VS$Y ~ mod1s.VS$fitted.values)$coefficients)
  weak.pvalue <- summary(mod1s.VS)$coefficients[2, 4]
  first.stage.f <- stats::qf(1 - weak.pvalue, 1, this.design$n - 2)
  if (first.stage.f < weak.instrument.f) return(list(
    seed = seed,
    model.2SLS = mod2s.VS,
    weak.pvalue = weak.pvalue,
    model.Newey = NULL,
    models.VS = NULL,
    models.eff = NULL,
    seconds = NA
  ))
  # An initial model based on Newey's 1993 efficient instruments.
  mod.init <- lmEff(y = data.VS$Y, incl = NULL, endog = data.VS$Z, excl = data.VS$X, iterations = 2, coef.names = c("(Intercept)", "Z"))

  models.VS <- vector("list", length(sel.bw))
  for (j in seq_along(sel.bw)) {
    start.mod <- mod.init
    if (j > 1) start.point <- tryCatch(models.VS[[j-1]]$sel$par, error = function(e) return(NA)) else start.point <- mod.init$coefficients
    tryel <- FALSE
    if (!all(is.finite(start.point))) {
      start.point <- mod.init$coefficients # If there is no numeric value in the previous iteration, try the optimal 2SLS
      tryel <- TRUE # And a more thorough initial value search should be carried out
    }
    if (j == 1) tryel <- TRUE
    start.mod$coefficients <- start.point
    start.mod$coefficients <- rbind(start.mod$coefficients, mod2s.VS)
    models.VS[[j]] <- estimateOneModelDesign1(start.mod = start.mod, data = data, sel.bw = sel.bw[j], try.ellipse = tryel, ...)
    cat("Seed ", .lead0(seed, 4), ", SEL bw = ", .trail0(sel.bw[j], 3),  ", ^theta_VS = ", paste(.trail0(models.VS[[j]]$sel$par, 3), collapse = " "), " (", j, "/", length(sel.bw), ")", "\n", sep = "")
  }

  tic1 <- Sys.time()
  if (CV) { # Cross-validated bw for pi
    pi.bw <- tryCatch(bw.CV(x = data$X, y = data$D), error = function(e) {cat("Seed", seed, "CV error:\n"); print(e); return(bw.CV(data$X))})
    u <- data.VS$Y - as.numeric(cbind(1, data.VS$Z) %*% models.VS[[1]]$sel$par)
    helper.bw <- tryCatch(bw.CV(x = as.matrix(data.VS[, c("X", "Z")]), y = u, same = TRUE), error = function(e) {cat("Seed", seed, "CV error:\n"); print(e); return(bw.CV(as.matrix(data.VS[, c("X", "Z")]), same = TRUE))})
  }

  models.eff <- NULL
  if ((length(pi.bw) > 0) && (length(helper.bw) > 0) && (length(helper) > 0) && (length(degree) > 0)) {
    hyperpar.df <- expand.grid(sel.bw = sel.bw, pi.bw = pi.bw, helper.bw = helper.bw, helper = helper, degree = degree)
    hyperpar.num <- as.matrix(as.data.frame(lapply(hyperpar.df, as.numeric)))
    hyperpar.sd <- apply(hyperpar.num, 2, stats::sd)
    hyperpar.sd[hyperpar.sd < 1e-8 | is.na(hyperpar.sd)] <- 1
    hyperpar.df[, c("theta0", "theta1")] <- NA
    models.eff <- vector("list", nrow(hyperpar.df))
    tryel <- TRUE # Always try ellipse search at the first iteration
    for (j in seq_len(nrow(hyperpar.df))) {
      start.mod <- mod.init
      if (j == 2 && all(is.finite(models.eff[[1]]$sel$par))) {
        start.mod$coefficients <- models.eff[[1]]$sel$par
      } else if (j > 2) {
        j.where.theta.notNA <- which(!is.na(hyperpar.df$theta1))
        done.so.far <- matrix(hyperpar.num[j.where.theta.notNA, ], ncol = 5)
        done.so.far <- sweep(done.so.far, 2, hyperpar.sd, "/")
        model.distance <- sweep(done.so.far, 2, hyperpar.num[j, ] / hyperpar.sd, "-")
        closest.model <- which.min(rowSums(model.distance^2))
        if (length(closest.model) > 0) { # If something was computed at all and the first two simulations did not fail
          if (all(is.finite(models.eff[[closest.model]]$sel$par))) start.mod$coefficients <- models.eff[[j.where.theta.notNA[closest.model]]]$sel$par
          tryel <- !all(is.finite(models.eff[[closest.model]]$sel$par)) # If all went smooth, no advanced initial value search is required
        }
      } else {
        start.mod$coefficients <- rbind(start.mod$coefficients, mod2s.VS, models.VS[[1]]$sel$par)
      }
      ec <- list(pi.bw = hyperpar.df$pi.bw[j], helper.bw = hyperpar.df$helper.bw[j],
                 helper = as.character(hyperpar.df$helper[j]), helper.degree = hyperpar.df$degree[j], PIT = TRUE)
      models.eff[[j]] <- estimateOneModelDesign1(start.mod = start.mod, data = data, sel.bw = hyperpar.df$sel.bw[j], eff.control = ec, try.ellipse = tryel, ...)
      hyperpar.df[j, c("theta0", "theta1")] <- models.eff[[j]]$sel$par
      cat("Seed ", .lead0(seed, 4), ", SEL bw = ", .trail0(hyperpar.df$sel.bw[j], 3),
          ", ^theta_ef = ", paste(.trail0(models.eff[[j]]$sel$par, 3), collapse = " "),
          " [b_pi = ", .trail0(hyperpar.df$pi.bw[j], 3),
          ", b_mu = ", .trail0(hyperpar.df$helper.bw[j], 3),
          ", help = ", as.character(hyperpar.df$helper[j]),
          ", deg = ", hyperpar.df$degree[j], "]",
          " (", j, "/", nrow(hyperpar.df), ") {closest: ",
          if (j > 2) j.where.theta.notNA[closest.model] else "", "}\n", sep = "")
    }
  }

  tic2 <- Sys.time()

  times <- as.numeric(c(VS = difftime(tic1, tic0, units = "secs"), Efficient = difftime(tic2, tic1, units = "secs")))

  return(list(
    seed = seed,
    design = this.design,
    model.2SLS = mod2s.VS,
    weak.pvalue = weak.pvalue,
    model.Newey = mod.init$coefficients,
    models.VS = models.VS,
    models.eff = models.eff,
    seconds = times
  ))
}

doTheRestSELContinuous <- function(SELmodel, do.SE = TRUE, do.restr = TRUE, CI.lev = c(0.90, 0.95, 0.99),
                                   print.inference = TRUE,
                                   ... # Passed to estimateOneModelDesign1
                                   ) {
  num.VS  <- length(SELmodel$models.VS)
  num.eff <- length(SELmodel$models.eff)
  if ((num.VS > 0) || (num.eff > 0)) {
    data <- generateData(n = SELmodel$design$n, distfun = stats::runif, propensity = SELmodel$design$propensity, sigma2X = SELmodel$design$sigma2X, seed = SELmodel$seed)
    data.VS <- data[as.logical(data$D), ]
    mod1s.VS <- stats::lm(Z ~ X, data = data.VS)
    mod2s.VS <- unname(stats::lm(data.VS$Y ~ mod1s.VS$fitted.values)$coefficients)
    mod.init <- lmEff(y = data.VS$Y, incl = NULL, endog = data.VS$Z, excl = data.VS$X, iterations = 2, coef.names = c("(Intercept)", "Z"))
  }

  if (num.VS > 0) {
    for (i in 1:num.VS) {
      submodel <- SELmodel$models.VS[[i]]
      this.do.SE <- is.null(submodel$vcov) & do.SE # Do we need to do anything? If something has been computed already, skip
      this.do.restr <- any(sapply(submodel$restricted, is.null)) & do.restr
      this.CI.lev <- if (is.null(submodel$CI) && (!is.null(CI.lev))) CI.lev else NULL
      mod.init$coefficients <- submodel$sel$par

      if (isTRUE(all(is.finite(mod.init$coefficients)))) {
        this.mod <- estimateOneModelDesign1(start.mod = mod.init, data = data, sel.bw = submodel$sel.bw, do.SE = this.do.SE, do.restr = this.do.restr, CI.lev = this.CI.lev, try.ellipse = FALSE, ...)
        SELmodel$models.VS[[i]] <- this.mod
        if (print.inference) {
          cat("Seed ", .lead0(SELmodel$seed, 5), ", ^theta_VS = ", paste(.trail0(this.mod$sel$par, 4), collapse = " "), "\n", sep = "")
          if (this.do.SE) cat("Standard errors:       (", paste(.trail0(sqrt(diag(this.mod$vcov)), 4), collapse = " "), ")\n", sep = "")
          if (this.do.restr) cat("Restricted true slope: ^theta = ", paste(.trail0(this.mod$restricted$`R:theta1=1`$par, 4), collapse = " "), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(this.mod$LR[2], 1)), "\n",
                                 "Restricted true const: ^theta = ", paste(.trail0(this.mod$restricted$`R:theta0=1`$par, 4), collapse = " "), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(this.mod$LR[3], 1)), "\n", sep = "")
          if(!is.null(this.CI.lev)) cat("Confidence interval points:", .trail0(c(this.mod$CI[rev(seq_along(CI.lev)), 1], this.mod$CI[seq_along(CI.lev), 2]), 4), "\n")
        }
      }
    }
  }

  if (num.eff > 0) {
    for (i in 1:num.eff) {
      submodel <- SELmodel$models.eff[[i]]
      this.do.SE <- is.null(submodel$vcov) & do.SE # Do we need to do anything?
      this.do.restr <- any(sapply(submodel$restricted, is.null)) & do.restr
      this.CI.lev <- if (is.null(submodel$CI) && (!is.null(CI.lev))) CI.lev else NULL
      this.eff.control <- submodel$eff.control
      mod.init$coefficients <- submodel$sel$par

      if (isTRUE(all(is.finite(mod.init$coefficients)))) {
        this.mod <- estimateOneModelDesign1(start.mod = mod.init, data = data, sel.bw = submodel$sel.bw, eff.control = this.eff.control, do.SE = this.do.SE, do.restr = this.do.restr, CI.lev = this.CI.lev, try.ellipse = FALSE, ...)
        SELmodel$models.eff[[i]] <- this.mod
        if (print.inference) {
          cat("Seed ", .lead0(SELmodel$seed, 5), ", ^theta_ef = ", paste(.trail0(this.mod$sel$par, 4), collapse = " "), "\n", sep = "")
          if (this.do.SE) cat("Standard errors:       (", paste(.trail0(sqrt(diag(this.mod$vcov)), 4), collapse = " "), ")\n", sep = "")
          if (this.do.restr) cat("Restricted true slope: ^theta = ", paste(.trail0(this.mod$restricted$`R:theta1=1`$par, 4), collapse = " "), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(this.mod$LR[2], 1)), "\n",
                                 "Restricted true const: ^theta = ", paste(.trail0(this.mod$restricted$`R:theta0=1`$par, 4), collapse = " "), ", p(LR) = ", sprintf("%1.3f", stats::pchisq(this.mod$LR[3], 1)), "\n", sep = "")
          if(!is.null(this.CI.lev)) cat("Confidence interval points:", .trail0(c(this.mod$CI[rev(seq_along(CI.lev)), 1], this.mod$CI[seq_along(CI.lev), 2]), 4), "\n")
        }
      }
    }
  }

  return(SELmodel)
}

#' Moment functions from the discrete design
#'
#' These functions are similar to \eqn{g(\theta, x)} described in the documentation of the \code{gmm} package.
#'
#' @param theta A parameter vector at which the moment function is evaluated.
#' @param data A data frame on which the moment function is evaluated.
#'
#' @return A matrix with two columns of the moment function values.
#' @export
#'
#' @examples
g.unconditional.vs <- function(theta, data) {
  U <- data$Y - cbind(1, data$Z) %*% theta
  ret <- cbind(U, U * data$X)
  ret[!is.na(U), ]
}

#' @param pihat A numeric vector contatining the predicted propensity score. Must be the same length as \code{nrow(data)}.
#' @param ystarhat A numeric vector contatining the predicted Y* (notation from the paper). Must be the same length as \code{nrow(data)}.
#' @export
#' @describeIn g.unconditional.vs Moment function incorporating the propensity score.
g.unconditional.eff <- function(theta, data, pihat, ystarhat) {
  Ztheta <- as.numeric(cbind(1, data$Z) %*% theta)
  Dg <- data$Y - Ztheta
  Dg[is.na(data$Y)] <- 0
  mu <- ystarhat - Ztheta
  rho <- (Dg - data$D * mu) / pihat + mu
  cbind(rho, rho * data$X)
}

#' Moment functions with conditioning
#'
#' @param theta A parameter vector at which the moment function is evaluated.
#' @param data A data frame on which the moment function is evaluated.
#'
#' Unlike \code{g.unconditional.vs} and \code{g.unconditional.eff}, these functions are designed for models with conditional moment restrictions.
#'
#' @return A matrix with two columns of the moment function values.
#' @export
#'
#' @examples
rho.complete.case <- function(theta, data) {
  Ztheta <- as.numeric(cbind(1, data$Z) %*% theta)
  Dg <- data$Y - Ztheta
  Dg[is.na(data$Y)] <- 0
  return(Dg)
}

#' @param pi.hat A numeric vector contatining the predicted propensity score. Must be the same length as \code{nrow(data)}.
#' @export
#' @describeIn rho.complete.case Moment function depending on the propensity score.
rho.complete.case.pi <- function(theta, data, pi.hat) {
  Ztheta <- as.numeric(cbind(1, data$Z) %*% theta)
  Dg <- data$Y - Ztheta
  Dg[is.na(data$Y)] <- 0
  return(Dg / pi.hat)
}

#' Moment function for efficient full-sample estimation
#'
#' @param theta A parameter vector at which the moment function is evaluated.
#' @param data A data frame on which the moment function is evaluated.
#' @param pi.hat A numeric vector contatining the predicted propensity score. Must be the same length as \code{nrow(data)}.
#' @param helper A string indicating which helper to use. See the vignette for a better explanation.
#' @param helper.predicted A numeric vector contatining the predicted helper value (notation from the paper). Must be the same length as \code{nrow(data)}.
#' @param pi.bw Smoothing bandwidth for propensity score estimation.
#' @param helper.bw Smoothing bandwidth(s) for propensity score estimation (the defult variable order is Z, X).
#' @param helper.degree Integer, passed to \code{kernelSmooth}.
#' @param PIT Logical: transform the values in \code{data} column-wise prior to smoothing?
#'
#' @return A numeric vector with moment function values.
#' @export
#'
#' @examples
#' data <- generateData()
#' data.VS <- data[as.logical(data$D), ]
#' t0 <- c(1, 1) # True parameters
#' Ztheta <- t0[1] + t0[2] * data$Z
#' b.pi <- bw.CV(x = pit(data$X), y = data$D)
#' b.muY  <- bw.CV(x = apply(as.matrix(data.VS[, c("Z", "X"), ]), 2, pit), y = data.VS$Y)
#' b.mug <- bw.CV(x = apply(as.matrix(data.VS[, c("Z", "X"), ]), 2, pit), y = data.VS$U)
#' r.Ystar <- rho.full.sample(t0, data, helper = "Ystar", pi.bw = b.pi, helper.bw = b.muY)
#' r.DY    <- rho.full.sample(t0, data, helper = "DY"   , pi.bw = b.pi, helper.bw = b.muY)
#' r.gstar <- rho.full.sample(t0, data, helper = "gstar", pi.bw = b.pi, helper.bw = b.mug)
#' r.Dg    <- rho.full.sample(t0, data, helper = "Dg"   , pi.bw = b.pi, helper.bw = b.mug)
#'
#' cols <- rainbow(4, end = 0.8, v = 0.6, alpha = 0.6)
#' plot(NULL, NULL, xlim = c(0, 1), ylim = c(-10, 10), bty = "n", xlab = "X", ylab = "rho")
#' points(data$X, r.Ystar, col = cols[1], pch = 16, cex = 0.75)
#' points(data$X, r.DY,    col = cols[2], pch = 16, cex = 0.75)
#' points(data$X, r.gstar, col = cols[3], pch = 16, cex = 0.75)
#' points(data$X, r.Dg,    col = cols[4], pch = 16, cex = 0.75)
#' plot(r.Ystar, r.gstar, bty = "n", col = data$D + 1)
#'
#' # Smoothing linear functions with locally constant kernels is a bad idea
#' Zthetahat0 <- kernelSmooth(data.VS[, c("Z", "X")], Ztheta[as.logical(data$D)],
#'   xout = data[, c("Z", "X")], bw = c(0.5, 0.15), degree = 0)
#' Zthetahat1 <- kernelSmooth(data.VS[, c("Z", "X")], Ztheta[as.logical(data$D)],
#'   xout = data[, c("Z", "X")], bw = c(1.25, 0.25), degree = 1)
#' plot(data$Z, Zthetahat1, bty = "n")
#' points(data$Z, Zthetahat0, col = 2)
rho.full.sample <- function(theta, data,
                            pi.hat = NULL,
                            helper = c("gstar", "Ystar", "Dg", "DY", "mu"),
                            helper.predicted = NULL,
                            pi.bw = NULL,
                            helper.bw = NULL,
                            helper.degree = 0,
                            PIT = TRUE) {
  helper <- helper[1]
  Ztheta <- as.numeric(cbind(1, data$Z) %*% theta)
  Dg <- data$Y - Ztheta
  Dg[is.na(data$Y)] <- 0
  if (is.null(pi.hat)) {
    if (PIT) X0 <- pit(data$X) else X0 <- data$X
    pi.hat <- kernelSmooth(x = X0, y = data$D, bw = pi.bw, degree = 0)
    pi.hat[pi.hat < 1 / nrow(data)] <- 1 / nrow(data)
    pi.hat[pi.hat > 1 - 1 / nrow(data)] <- 1 - 1 / nrow(data)
  }
  if (is.null(helper.predicted)) {
    ZX <- as.matrix(data[, c("Z", "X")])
    # For estimation of mu while using only VS for training, we need to use only (Z, X) from the VS as the grid
    if (helper %in% c("Ystar", "gstar")) {
      ZX.VS <- ZX[as.logical(data$D), ]
      if (PIT) {
        ZX.support <- do.call(cbind, lapply(seq_len(ncol(ZX)), function(i) pit(x = ZX.VS[, i], xout = ZX[, i])))
      } else {
        ZX.support <- ZX
      }
      ZX.VS <- ZX.support[as.logical(data$D), ]
    } else if (helper %in% c("DY", "Dg")) {
      if (PIT) ZX <- apply(ZX, 2, pit)
    }
    if (helper == "Ystar") {
      Y.VS <- data$Y[as.logical(data$D)]
      Ystar.hat <- kernelSmooth(x = ZX.VS, y = Y.VS,    xout = ZX.support, bw = helper.bw, degree = helper.degree)
      mu.hat <- Ystar.hat - Ztheta
    } else if (helper == "gstar") {
      g.VS <- Dg[as.logical(data$D)]
      mu.hat    <- kernelSmooth(x = ZX.VS, y = g.VS,    xout = ZX.support, bw = helper.bw, degree = helper.degree)
    } else if (helper == "DY") {
      Ystar.hat <- kernelSmooth(x = ZX,    y = data$DY, xout = NULL,       bw = helper.bw, degree = helper.degree) / pi.hat
      mu.hat <- Ystar.hat - Ztheta
    } else if (helper == "Dg") {
      mu.hat    <- kernelSmooth(x = ZX,    y = Dg,      xout = NULL,       bw = helper.bw, degree = helper.degree) / pi.hat
    }
  } else { # A helper was passed
    if (helper == "mu") mu.hat <- helper.predicted else stop("The predicted helper should be 'mu'.")
  }
  rho <- (Dg - data$D * mu.hat) / pi.hat + mu.hat
  return(rho)
}

sampleEllipse <- function(centre, vcov, radius, angles = rep(seq(-45, 45, 15), 2) + rep(c(0, 180), each = 7)) {
  # Based on car::ellipse
  angles.rad <- angles / 180 * pi
  if (!isTRUE(is.vector(centre) && length(centre) == 2))
    stop("centre must be a vector of length 2")
  if (!isTRUE(is.matrix(vcov) && all(dim(vcov) == 2)))
    stop("vcov must be a 2 by 2 matrix")
  if (max(abs(vcov - t(vcov)))/max(abs(vcov)) > 1e-10)
    stop("vcov must be a symmetric matrix")
  circle.coords <- cbind(cos(angles.rad), sin(angles.rad))
  Q <- chol(vcov, pivot = TRUE)
  order <- order(attr(Q, "pivot"))
  if (length(radius) == 1) {
    ellipse <- t(centre + radius * t(circle.coords %*% Q[, order]))
  } else {
    ellipse <- lapply(radius, function(r) t(centre + r * t(circle.coords %*% Q[, order])))
  }
  return(ellipse)
}

generateGeneric <- function(n = 10000,
                            ndiscrXin = 6, ncatXin = 2, ncontXin = 0, corXin = 0.3,
                            ndiscrXex = 2, ncatXex = 0, ncontXex = 0, corXex = 0.3,
                            corXinXex = 0.2,
                            ncat = 10, ndiscr = 2,
                            ndiscrZ = 1, ncatZ = 0, ncontZ = 0,
                            corUV = 0.5,
                            seed = 1
) {
  set.seed(seed)
  # Generating weakly correlated excluded and inlcuded instruments
  nXin <- ncontXin + ncatXin + ndiscrXin
  nXex <- ncontXex + ncatXex + ndiscrXex
  nZ <- ncontZ + ncatZ + ndiscrZ
  nexog <- nXin + nXex
  varXin <- if (nXin > 0) matrix(corXin, ncol = nXin, nrow = nXin) + diag(nXin) * (1 - corXin) else NULL
  varXex <- if (nXex > 0) matrix(corXex, ncol = nXex, nrow = nXex) + diag(nXex) * (1 - corXex) else NULL
  covXinXex <- if (nXin > 0 && nXex > 0) matrix(corXinXex, ncol = nXex, nrow = nXin) else NULL
  tnull <- function(x) if (is.null(x)) NULL else t(x)
  varX <- rbind(cbind(varXin, covXinXex), cbind(tnull(covXinXex), varXex))
  X <- matrix(stats::rnorm(n * nexog), ncol = nexog, nrow = n)
  if (ncol(X) > 1) {
    eX <- eigen(varX, symmetric = TRUE)
    X <- tnull(eX$vectors %*% diag(sqrt(eX$values)) %*% tnull(X))
  }
  colnames(X) <- c(if (ncontXin > 0) paste0("Xin.cont", 1:ncontXin) else NULL,
                    if (ncatXin > 0) paste0("Xin.cat", 1:ncatXin) else NULL,
                    if (ndiscrXin > 0) paste0("Xin.discr", 1:ndiscrXin) else NULL,
                    if (ncontXex > 0) paste0("Xex.cont", 1:ncontXex) else NULL,
                    if (ncatXex > 0) paste0("Xex.cat", 1:ncatXex) else NULL,
                    if (ndiscrXex > 0) paste0("Xex.discr", 1:ndiscrXex) else NULL)
  # Making some of the instruments categorical, with 20 values per category
  cat.ind <- grep("cat", colnames(X))
  if (length(cat.ind) > 0) X[, cat.ind] <- sapply(cat.ind, function(i) scale(as.numeric(cut(X[, i], breaks = c(-Inf, stats::quantile(X[, i], (1:(ncat-1))/ncat), Inf), labels = 1:ncat))))
  # Discretising some of the instruments, 3 categories for each
  discr.ind <- grep("discr", colnames(X))
  if (length(discr.ind) > 0) X[, discr.ind] <- sapply(discr.ind, function(i) scale(as.numeric(cut(X[, i], breaks = c(-Inf, stats::quantile(X[, i], (1:(ndiscr-1))/ndiscr), Inf), labels = 1:ndiscr))))
  Xin <- if (nXin > 0) X[, 1:nXin, drop = FALSE] else NULL
  Xex <- if (nXin > 0) X[, -(1:nXin), drop = FALSE] else X

  UV <- matrix(stats::rnorm(n * (nZ + 1)), ncol = nZ + 1, nrow = n)
  covUV <- matrix(corUV, ncol = 1, nrow = nZ)
  varUV <- rbind(cbind(diag(nZ), covUV), c(covUV, 1))
  eUV <- eigen(varUV, symmetric = TRUE)
  UV <- tnull(eUV$vectors %*% diag(sqrt(eUV$values)) %*% tnull(UV))
  U <- UV[, ncol(UV), drop = FALSE]
  V <- UV[, -ncol(UV), drop = FALSE]
  colnames(V) <- paste0("V", 1:nZ)

  sigma2X <- function(x) {
    x <- x + 3
    return(0.1 + 0.2*sum(abs(x)) + 0.3*sum(x^2))
  }
  VarU.X <- apply(X, 1,  sigma2X)
  Usked <- U * sqrt(VarU.X) * 1.5
  # Grouping instruments by endogenous variable; some are going to be more relevant than others
  corresp.Xex.Z <- if (nZ > 1) as.integer(cut(1:nXex, breaks = nZ, labels = 1:nZ)) else rep(1, nXex)
  Z <- do.call(cbind, lapply(1:nZ, function(i) (if (nXin > 0) rowMeans(Xin) else 0) + 1 * rowSums(Xex[, corresp.Xex.Z == i, drop = FALSE]) + 0.1 * rowSums(Xex[, corresp.Xex.Z != i, drop = FALSE]) + V[, i] * sqrt(2 + sum(corresp.Xex.Z == i))
                             )) # The variance of the reuced-form error is proportional to the sum of the absolute values of the coefficients; the variance of the regressors is more or less the same
  # Now making the variance of Z close to unity
  Z <- sweep(Z, 2, apply(Z, 2, stats::sd), "/")
  colnames(Z) <- paste0("Z", 1:nZ)

  Ystar <- 1 + (if (nXin > 0) rowSums(Xin) else 0) + rowSums(Z) + Usked # Structural equation
  # The propensity score in this design will be simple: since all regressors are symmetrically distributed around zero,
  # it be proportional to the share of positive observations (between 0.95 if everything is negative to 0.70 if everything is positive)
  # The skedastic function is higher on the upper end, so the missingness will appear more frequently there
  prop.score <- function(x) 0.7 + 0.25 * mean(x >= 0)
  piXZ <- apply(cbind(Xin, Xex, Z), 1, prop.score)
  D <- stats::runif(n) < piXZ
  Y <- Ystar
  Y[!D] <- NA
  DY <- Y
  DY[!D] <- 0
  alldata <- if (nXin > 0) data.frame(Xin, Xex, Z, D = as.numeric(D), pi = piXZ, Ystar, Y, DY, U, Usked, V, VarU.X) else data.frame(Xex, Z, D = as.numeric(D), pi = piXZ, Ystar, Y, DY, U, Usked, V, VarU.X)
  return(alldata)
}


generateCMRData <- function(n = 10000,
                            ndiscrX = 6, ncatX = 2, ncontX = 0, corX = 0.2,
                            ncat = 10, ndiscr = 2,
                            seed = 1
) {
  set.seed(seed)
  # Generating weakly correlated excluded and included instruments
  nX <- ncontX + ncatX + ndiscrX
  varX <- matrix(corX, ncol = nX, nrow = nX) + diag(nX) * (1 - corX)
  X <- matrix(stats::rnorm(n * nX), ncol = nX, nrow = n)
  if (ncol(X) > 1) {
    eX <- eigen(varX, symmetric = TRUE)
    X <- t(eX$vectors %*% diag(sqrt(eX$values)) %*% t(X))
  }
  colnames(X) <- c(if (ncontX > 0) paste0("X.cont", 1:ncontX) else NULL,
                   if (ncatX > 0) paste0("X.cat", 1:ncatX) else NULL,
                   if (ndiscrX > 0) paste0("X.discr", 1:ndiscrX) else NULL
                   )
  # Making some of the instruments categorical, with 20 values per category
  cat.ind <- grep("cat", colnames(X))
  if (length(cat.ind) > 0) X[, cat.ind] <- sapply(cat.ind, function(i) scale(as.numeric(cut(X[, i], breaks = c(-Inf, stats::quantile(X[, i], (1:(ncat-1))/ncat), Inf), labels = 1:ncat))))
  # Discretising some of the instruments, 3 categories for each
  discr.ind <- grep("discr", colnames(X))
  if (length(discr.ind) > 0) X[, discr.ind] <- sapply(discr.ind, function(i) scale(as.numeric(cut(X[, i], breaks = c(-Inf, stats::quantile(X[, i], (1:(ndiscr-1))/ndiscr), Inf), labels = 1:ndiscr))))

  U <- stats::rnorm(n)
  sigma2X <- function(x) {
    x <- x + 3
    return(0.1 + 0.2*sum(abs(x)) + 0.3*sum(x^2))
  }
  VarU.X <- apply(X, 1,  sigma2X)
  Usked <- U * sqrt(VarU.X)

  Ystar <- 1 + rowSums(X) + Usked # Structural equation
  # The propensity score in this design will be simple: since all regressors are symmetrically distributed around zero,
  # it be proportional to the share of positive observations (between 0.95 if everything is negative to 0.70 if everything is positive)
  # The skedastic function is higher on the upper end, so the missingness will appear more frequently there

  # The propensity score in this design will be simple: since all regressors are symmetrically distributed around zero,
  # it be proportional to the share of positive observations (between 0.95 if everything is negative to 0.70 if everything is positive)
  # The skedastic function is higher on the upper end, so the missingness will appear more frequently there
  prop.score <- function(x) 0.7 + 0.25 * mean(x >= 0)
  piX <- apply(X, 1, prop.score)
  D <- stats::runif(n) < piX
  Y <- Ystar
  Y[!D] <- NA
  DY <- Y
  DY[!D] <- 0
  alldata <- data.frame(X, Ystar, D = as.numeric(D), pi = piX, Y, DY, U, Usked, VarU.X)
  return(alldata)
}
