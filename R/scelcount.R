#' Self-concordant empirical likelihood for a vector mean with counts
#'
#' Implements the empirical-likelihood-ratio test for the mean of the coordinates of \code{z}
#' (with the hypothesised value \code{mu}). The counts need not be integer;
#' in the context of local likelihoods, they can be kernel observation weights.
#'
#' @source This original code was written for \insertCite{owen2013self}{smoothemplik}
#' and [published online](https://artowen.su.domains/empirical/) by Art B. Owen
#' (February 2014) and slightly reworked to contain fewer inner functions and loops.
#'
#' @details
#' Negative weights are not allowed. They could be useful in some applications, but they can destroy
#' convexity or even boundedness. They also make the Newton step fail to be of least squares type.
#'
#' This function relies on the improved computational strategy for the empirical likelihood.
#' The search of the lambda multipliers is carried out via a dampened Newton method with guaranteed
#' convergence owing to the fact that the log-likelihood is replaced by its Taylor approximation
#' of any desired order (default: 4, the minimum value that ensures self-concordance).
#'
#' Tweak \code{alpha} and \code{beta} with extreme caution. See \insertCite{boyd2004convex}{smoothemplik},
#' pp. 464--466 for details. It is necessary that \code{0 < alpha < 1/2} and \code{0 < beta < 1}.
#' \code{alpha = 0.03} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).
#'
#' The argument names, except for \code{lambda.init}, are matching the original names in Art B. Owen's implementation.
#' The highly optimised one-dimensional counterpart, \code{weightedEL}
#'
#' @param z A numeric vector or a matrix with one data vector per column.
#' @param ct A numeric vector of non-negative counts.
#' @param mu Hypothesised mean, default (0 ... 0) in R^ncol(z)
#' @param lambda.init Starting lambda, default (0 ... 0)
#' @param return.weights Logical: if \code{TRUE}, returns the empirical probabilities. Default is memory-saving (\code{FALSE}).
#' @param lower Lower cut-off for [logTaylor()], default 1/nrow(z)
#' @param upper Upper cutoff for [logTaylor()], default Inf
#' @param order Positive integer such that the Taylor approximation of this order to log(x) is self-concordant; usually 4 or higher. Passed to [logTaylor()].
#' @param wttol Weight tolerance for counts to improve numerical stability
#' @param thresh Convergence threshold for log-likelihood (the default is aggressive)
#' @param itermax Upper bound on number of Newton steps (seems ample)
#' @param alpha Backtracking line search parameter: acceptance of a decrease in function value by ALPHA*f of the prediction
#'   based on the linear extrapolation.
#' @param beta Backtracking line search reduction factor. 0.1 corresponds to a very crude search, 0.8 corresponds
#'   to a less crude search.
#' @param backeps Backtrack threshold: the search can miss by this much. Consider setting it to 1e-10
#'   if backtracking seems to be failing due to round-off.
#' @param verbose Logical: print output diagnostics?
#'
#' @return A list with the following values:
#' \describe{
#'     \item{logelr}{Log of empirical likelihood ratio (equal to 0 if the hypothesised mean is equal to the sample mean)}
#'     \item{lam}{Vector of Lagrange multipliers}
#'     \item{wts}{Observation weights/probabilities (vector of length n)}
#'     \item{converged}{\code{TRUE} if algorithm converged. \code{FALSE} usually means that mu is not in the convex hull of the data. Then, a very small likelihood is returned (instead of zero).}
#'     \item{iter}{Number of iterations taken.}
#'     \item{ndec}{Newton decrement (see Boyd & Vandenberghe).}
#'     \item{gradnorm}{Norm of the gradient of log empirical likelihood.}
#' }
#'
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [logTaylor()], [weightedEL()]
#'
#' @examples
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' cemplik(earth, mu = 5.517, verbose = TRUE) # 5.517 is the modern accepted value
#'
#' # Linear regression through empirical likelihood
#' coef.lm <- coef(lm(mpg ~ hp + am, data = mtcars))
#' xmat <- cbind(1, as.matrix(mtcars[, c("hp", "am")]))
#' yvec <- mtcars$mpg
#' foc.lm <- function(par, x, y) {  # The sample average of this
#'   resid <- y - drop(x %*% par)   # must be 0
#'   resid * x
#' }
#' minusEL <- function(par) -cemplik(foc.lm(par, xmat, yvec), itermax = 10)$logelr
#' coef.el <- optim(c(mean(yvec), 0, 0), minusEL)$par
#' abs(coef.el - coef.lm) / coef.lm  # Relative difference
#'
#' # Likelihood ratio testing without any variance estimation
#' # Define the profile empirical likelihood for the coefficient on am
#' minusPEL <- function(par.free, par.am)
#'   -cemplik(foc.lm(c(par.free, par.am), xmat, yvec), itermax = 20)$logelr
#' # Constrained maximisation assuming that the coef on par.am is 3.14
#' coef.el.constr <- optim(coef.el[1:2], minusPEL, par.am = 3.14)$par
#' print(-2 * cemplik(foc.lm(c(coef.el.constr, 3.14), xmat, yvec))$logelr)
#' # Exceeds the critical value qchisq(0.95, df = 1)
#' @export
cemplik <- function(z, ct = NULL, mu = NULL, lambda.init = NULL,
                    return.weights = FALSE, lower = NULL, upper = NULL,
                    order = 4L, wttol = 0.001,
                    thresh = 1e-16, itermax = 100L, verbose = FALSE,
                    alpha = 0.3, beta = 0.8, backeps = 0) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- nrow(z)
  d <- ncol(z)
  if (length(mu) == 0) mu <- rep(0, d)
  if (length(mu) != d) stop("The length of mu must be the same as the dimension of z.")
  if (length(lambda.init) == 0) lambda.init <- rep(0, d)
  if (length(lambda.init) != d) stop("The length of mu must be the same as the dimension of z.")

  if (is.null(lower)) lower <- rep(1/n, n)
  if (is.null(upper)) upper <- rep(Inf, n)

  if (is.null(ct)) ct <- rep(1, n)
  if (min(ct) < 0) stop("Negative weights are not allowed.")
  if (any(0 < ct & ct < wttol)) {
    warning(paste("Positive counts below", wttol, "have been replaced by zero."))
    ct[ct < wttol] <- 0
  }
  nonz <- (ct > 0) # These are observations with nonzero weight used in the Newton step below
  if (sum(ct) <= 0) stop("Total weight must be positive.")

  weightedELCPP(z = z, ct = ct, mu = mu, lambda_init = lambda.init,
                return_weights = return.weights, lower = lower, upper = upper,
                order = order, wttol = wttol, thresh = thresh,
                itermax = itermax, verbose = verbose,
                alpha = alpha, beta = beta, backeps = backeps)
}

#' Compute empirical likelihood on a trajectory
#'
#' @param z Passed to \code{cemplik}.
#' @param ct Passed to \code{cemplik}.
#' @param mu0 Starting point of trajectory
#' @param mu1 End point of trajectory
#' @param N Number of segments into which the path is split (i. e. \code{N+1} steps are used).
#' @param verbose Logical: report iteration data?
#' @param verbose.solver Logical: report internal iteration data from the optimiser? Very verbose.
#' @param ... Passed to \code{cemplik}.
#'
#' This function does not accept the starting lambda because it is much faster (3--5 times)
#' to reuse the lambda from the previous iteration.
#'
#' @return A matrix with one row at each mean from mu0 to mu1 and a column for each EL return value (except EL weights).
#' @export
#'
#' @examples
#' # Plot 2.5 from Owen (2001)
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' cemplik(earth, mu = 5.1,  verbose = TRUE)
#' logELR <- ctracelr(earth, mu0 = 5.1, mu1 = 5.65, N = 55, verbose = TRUE)
#' hist(earth, breaks = seq(4.75, 6, 1/8))
#' plot(logELR[, 1], exp(logELR[, 2]), bty = "n", type = "l",
#'      xlab = "Earth density", ylab = "ELR")
#' # TODO: why is there non-convergence in row 0?
#'
#' # Two-dimensional trajectory
#' set.seed(1)
#' xy <- matrix(rexp(200), ncol = 2)
#' logELR2 <- ctracelr(xy, mu0 = c(0.5, 0.5), mu1 = c(1.5, 1.5), N = 100)
ctracelr <- function(z, ct = NULL, mu0, mu1, N = 5, verbose = FALSE,
                     verbose.solver = FALSE, ...) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  d <- ncol(z)

  ans <- matrix(0, nrow = N + 1, ncol = d + 1 + d + 4)
  colnames(ans) <- c(paste("z", 1:d, sep = "."),
                     "logelr",
                     paste("lambda", 1:d, sep = "."),
                     c("conv", "iter", "decr", "gnorm"))
  lam0 <- NULL
  for (i in 0:N) {
    mui <- (i*mu1 + (N-i)*mu0) / N
    if (i > 0) lam0 <- if (all(is.finite(ansi$lam))) drop(ansi$lam) else NULL
    ansi <- cemplik(z = z, ct = ct, mu = mui, lambda.init = lam0, verbose = verbose.solver, ...)
    ans[i + 1, ] <- c(mui, ansi$logelr, ansi$lam, ansi$converged, ansi$iter, ansi$ndec, ansi$gradnorm)
    if (verbose) cat("Point ", i, "/", N, ", ", if (ansi$converged == 0) "NOT " else "", "converged",
                     ", log(ELR) = ", ansi$logelr, "\n", sep = "")
  }
  return(ans)
}


#' Modified minus logarithm with derivatives
#'
#' @param x Numeric vector for which approximated logarithm is to be computed.
#' @param order Positive integer: Taylor approximation order. If \code{NA}, returns \code{log(x)} or its derivative.
#' @param lower Lower threshold below which approximation starts; can be a scalar of a vector of the same length as \code{x}.
#' @param upper Upper threshold above which approximation starts; can be a scalar of a vector of the same length as \code{x}.
#' @param der Non-negative integer: 0 yields the function, 1 and higher yields derivatives
#'
#' @details
#' Provides a family of alternatives to -log() and derivative thereof in order to attain self-concordance and
#' computes the modified negative logarithm and its first derivatives.
#' For lower <= x <= upper, returns just the logarithm. For x < lower and x > upper, returns the Taylor approximation of the given \code{order}.
#' 4th order is the lowest that gives self concordance.
#'
#' @return A numeric matrix with \code{(order+1)} columns containing the values of the modified log and its derivatives.
#' @export
#'
#' @examples
#' x <- seq(0.01^0.25, 2^0.25, length.out = 51)^4 - 0.11 # Denser where |f'| is higher
#' plot(x,  log(x)); abline(v = 0, lty = 2) # Observe the warning
#' lines(x, logTaylor(x, lower = 0.2), col = 2)
#' lines(x, logTaylor(x, lower = 0.5), col = 3)
#' lines(x, logTaylor(x, lower = 1, upper = 1.2, order = 6), col = 4)
#'
#' # Substitute log with its Taylor approx. around 1
#' x <- seq(0.1, 2, 0.05)
#' ae <- abs(sapply(2:6, function(o) log(x) - logTaylor(x, lower=1, upper=1, order=o)))
#' matplot(x[x!=1], ae[x!=1,], type = "l", log = "y", lwd = 2,
#'   main = "Abs. trunc. err. of Taylor expansion at 1", ylab = "")
logTaylor <- function(x, lower = NULL, upper = NULL, der = 0, order = 4) {
  n <- length(x)
  if (is.null(lower)) lower <- rep(1/n, n)
  if (is.null(upper)) upper <- rep(Inf, n)
  logTaylorCPP(x, lower, upper, der, order)
}


#' Least-squares regression via SVD
#'
#' @param x Model matrix.
#' @param y Response vector.
#' @param rel.tol Relative zero tolerance for generalised inverse via SVD.
#' @param abs.tol Absolute zero tolerance for generalised inverse via SVD.
#'
#' Newton steps for many empirical likelihoods are of least-squares type.
#' Denote \eqn{x^+} to be the generalised inverse of \code{x}.
#' If SVD algorithm failures are encountered, it sometimes helps to try
#' \code{svd(t(x))} and translate back. First check to ensure that
#' \code{x} does not contain \code{NaN}, or \code{Inf}, or \code{-Inf}.
#'
#' The tolerances are used to check the closeness of singular values to zero. The values of the
#' singular-value vector \code{d} that are less than
#' \code{max(rel.tol * max(d), abs.tol)} are set to zero.
#'
#' @return A vector of coefficients.
#' @export
#'
#' @examples
#' b.svd <- svdlm(x = cbind(1, as.matrix(mtcars[, -1])), y = mtcars[, 1])
#' b.lm  <- coef(lm(mpg ~ ., data = mtcars))
#' b.lm - b.svd  # Negligible differences
svdlm <- function(x, y, rel.tol = 1e-9, abs.tol = 1e-100) {
  if (is.null(dim(x))) x <- as.matrix(x)
  svdlmCPP(x = x, y = y, rel_tol = rel.tol, abs_tol = abs.tol)
}


#' d-th derivative of the k-th-order Taylor expansion of log(x)
#'
#' @param x Numeric: a vector of points for which the logarithm is to be evaluated
#' @param a Scalar: the point at which the polynomial approximation is computed
#' @param k Non-negative integer: maximum polynomial order in the Taylor expansion
#'   of the original function. \code{k = 0} returns a constant.
#' @param d Non-negative integer: derivative order
#'
#' Note that this function returns the d-th derivative of the k-th-order Taylor expansion, not the
#' k-th-order approximation of the d-th derivative. Therefore, the degree of the resulting polynomial
#' is \eqn{d-k}{d-k}.
#'
#' @return The approximating Taylor polynomial around \code{a} of the order \code{d-k} evaluated at \code{x}.
#' @export
#'
#' @examples
#' cl <- rainbow(9, end = 0.8, v = 0.8, alpha = 0.8)
#' a <- 1.5
#' x <- seq(a*2, a/2, length.out = 101)
#' f <- function(x, d = 0)  if (d == 0) log(x) else ((d%%2 == 1)*2-1) * 1/x^d * gamma(d)
#' par(mfrow = c(2, 3), mar = c(2, 2, 2.5, 0.2))
#' for (d in 0:5) {
#' y <- f(x, d = d)
#' plot(x, y, type = "l", lwd = 7, bty = "n", ylim = range(0, y),
#'        main = paste0("d^", d, "/dx^", d, " Taylor(Log(x))"))
#'   for (k in 0:8) lines(x, tlog(x, a = a, k = k, d = d), col = cl[k+1], lwd = 1.5)
#'   points(a, f(a, d = d), pch = 16, cex = 1.5, col = "white")
#' }
#' legend("topright", as.character(0:8), title = "Order", col = cl, lwd = 1)
tlog <- function(x, a = as.numeric(c(1.0)), k = 4L, d = 0L) tlogCPP(x, a, k, d)
