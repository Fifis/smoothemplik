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
#' Tweak \code{ALPHA} and \code{BETA} with extreme caution. See \insertCite{boyd2004convex}{smoothemplik},
#' pp. 464--466 for details. It is necessary that \code{0 < ALPHA < 1/2} and \code{0 < BETA < 1}.
#' \code{ALPHA = 0.03} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).
#'
#' The argument names are matching the original names in Art B. Owen's implementation.
#' The highly optimised one-dimensional counterpart, \code{weightedEL}
#'
#' @param z A numeric vector or a matrix with one data vector per column.
#' @param ct A numeric vector of non-negative counts.
#' @param mu Hypothesised mean, default (0 ... 0) in R^ncol(z)
#' @param lam Starting lambda, default (0 ... 0)
#' @param eps Lower cut-off for \code{logTaylor()}, default 1/nrow(z)
#' @param M Upper cutoff for \code{logTaylor()}, default Inf
#' @param order Positive integer such that the Taylor approximation of this order to log(x) is self-concordant; usually 4 or higher. Passed to [logTaylor()].
#' @param wttol Weight tolerance for counts to improve numerical stability
#' @param thresh Convergence threshold for log-likelihood (the default is aggressive)
#' @param itermax Upper bound on number of Newton steps (seems ample)
#' @param ALPHA Backtracking line search parameter: acceptance of a decrease in function value by ALPHA*f of the prediction
#'   based on the linear extrapolation.
#' @param BETA Backtracking line search reduction factor. 0.1 corresponds to a very crude search, 0.8 corresponds
#'   to a less crude search.
#' @param BACKEPS Backtrack threshold: the search can miss by this much. Consider setting it to 1e-10
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
#' @references
#' \insertAllCited{}
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
cemplik <- function(z, ct = NULL, mu = NULL,
                    lam = NULL, eps = NULL, M = Inf,
                    order = 4, wttol = 0.001,
                    thresh = 1e-30, itermax = 100, ALPHA = 0.3, BETA = 0.8, BACKEPS = 0,
                    verbose = FALSE) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  n <- nrow(z)
  d <- ncol(z)

  if (is.null(mu)) mu <- rep(0, d)
  if (length(mu) != d) stop("The length of mu must be the same as the dimension of z.")
  z <- t(t(z) - mu) # Faster than z <- sweep(z, 2, mu, "-")

  if (is.null(eps)) eps <- 1/n

  if (is.null(ct)) ct <- rep(1, n)
  if (min(ct) < 0) stop("Negative weights are not allowed.")

  if (any(0 < ct & ct < wttol)) {
    warning(paste("Positive counts below", wttol, "have been replaced by zero."))
    ct[ct < wttol] <- 0
  }
  nonz <- (ct > 0) # These are observations with nonzero weight used in the Newton step below

  if (sum(ct) <= 0) stop("Total weight must be positive.")

  # Initialising at lambda = 0
  init <- -logTaylor(rep(1, n), eps = eps, M = M, order = order, der = 2) # i.e. fn and derivatives at lam = 0
  init <- init * ct # Multiply every column: weights appear OUTSIDE of logTaylor
  # If the initial lambda is supplied, compare lam = 0 with it and pick the best
  if (!is.null(lam)) {
    init1 <- -logTaylor(1 + z %*% lam, eps = eps, M = M, order = order, der = 2)
    init1 <- init1 * ct
    if (sum(init1[, 1]) < sum(init[, 1])) { # Minimising --> lower is better
      init <- init1
    } else {
      lam <- rep(0, d) # Ignore the supplied lambda
    }
  }
  if (is.null(lam)) lam <- rep(0, d)

  # Initial f, g
  fold <- sum(init[, 1])
  gold <- apply(z * init[, 2], 2, sum)

  converged <- FALSE
  iter <- 0
  oldvals <- init

  while (!converged) {
    iter <- iter + 1

    # Get Newton step
    rootllpp <- sqrt(oldvals[, 3]) # sqrt 2nd deriv of -llog lik
    zt <- z * rootllpp # Multiplying each column is faster
    yt <- oldvals[, 2] / rootllpp
    step <- -svdlm(zt[nonz, ], yt[nonz]) #  SVD is more reliable than step <- -lm( yt~zt-1 )$coef
    # Also, we kick out the ones with ct=0 to avoid unnecessary NaNs
    backtrack <- FALSE
    tt <- 1 # usually called t, but R uses t for transpose
    while (!backtrack) {
      newlam <- lam + tt*step
      newvals <- -logTaylor(1 + z %*% newlam, eps = eps, M = M, order = order, der = 2)
      newvals <- newvals * ct
      fnew <- sum(newvals[, 1])
      targ <- fold + ALPHA * tt * sum(gold*step) + BACKEPS # (BACKEPS for round-off, should not be needed)
      if (fnew <= targ) { # Backtracking has converged
        backtrack <- TRUE
        oldvals <- newvals
        fold <- fnew
        gold <- apply(z*oldvals[, 2], 2, sum)
        # take the step
        lam <- newlam
      } else {
        tt <- tt * BETA
      }
    }

    # Newton decrement and gradient norm
    ndec <- sqrt(sum((step*gold)^2))
    gradnorm <- sqrt(sum(gold^2))

    if (verbose) cat("Iter. ", iter, ", lambda = ", paste0(formatC(lam, digits = 8), collapse = " "),
                     ", f = ", fold, ", |grad| = ", gradnorm, ", decr = ", ndec, "\n", sep = "")
    converged <- (ndec^2 <= thresh)
    if (iter > itermax) break
  }

  zlam <- drop(z %*% lam)
  wts <- (ct / sum(ct)) / (1 + zlam) # Without the weights, the numerator was 1/n
  logelr <- -sum(ct * logTaylor(1 + zlam, eps = eps, M = M, order = order, der = 0))

  return(list(logelr = logelr, lam = lam, wts = wts, converged = converged,
              iter = iter, ndec = ndec, gradnorm = gradnorm))
}

#' Compute empirical likelihood on a trajectory
#'
#' @param z Passed to \code{cemplik}.
#' @param ct Passed to \code{cemplik}.
#' @param mu0 Starting point of trajectory
#' @param mu1 End point of trajectory
#' @param N Number of segments into which the path is split (i. e. \code{N+1} steps are used).
#' @param verbose Logical: report iteration data?
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
#' logELR <- ctracelr(earth, mu0 = 5.1, mu1 = 5.65, N = 100, verbose = TRUE)
#' hist(earth, breaks = seq(4.75, 6, 1/8))
#' plot(logELR[, 1], exp(logELR[, 2]), bty = "n", type = "l",
#'      xlab = "Earth density", ylab = "ELR")
#' # TODO: why is there non-convergence in row 59?
#'
#' # Two-dimensional trajectory
#' set.seed(1)
#' xy <- matrix(rexp(200), ncol = 2)
#' logELR2 <- ctracelr(xy, mu0 = c(0.5, 0.5), mu1 = c(1.5, 1.5), N = 100)
ctracelr <- function(z, ct = NULL, mu0, mu1, N = 5, verbose = FALSE, ...) {
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
    ansi <- cemplik(z = z, ct = ct, mu = mui, lam = lam0, ...)
    ans[i + 1, ] <- c(mui, ansi$logelr, ansi$lam, ansi$converged, ansi$iter, ansi$ndec, ansi$gradnorm)
    if (verbose) cat("Point ", i, "/", N, ", ", if (ansi$converged == 0) "NOT " else "", "converged",
                     ", log(ELR) = ", ansi$logelr, "\n", sep = "")
  }
  return(ans)
}

# Derivatives of the logarithm
dlog <- function(x, d = 0)
  if (d == 0) log(x) else ((d%%2 == 1)*2-1) * 1/x^d * gamma(d)

#' Modified minus logarithm with derivatives
#'
#' @param x Numeric vector for which approximated logarithm is to be computed.
#' @param order Positive integer: Taylor approximation order. If \code{NA}, returns \code{log(x)} or its derivative.
#' @param eps Lower threshold below which approximation starts; can be a scalar of a vector of the same length as \code{x}.
#' @param M Upper threshold above which approximation starts; can be a scalar of a vector of the same length as \code{x}.
#' @param der Non-negative integer: 0 yields the function, 1 and higher yields derivatives
#' @param drop If \code{TRUE} and \code{x} is a vector, the output is a vector and not a 1-column matrix.
#'
#' @details
#' Provides a family of alternatives to -log() and derivative thereof in order to attain self-concordance and
#' computes the modified negative logarithm and its first derivatives.
#' For eps <= x <= M, returns just the logarithm. For x < eps and x > M, returns the Taylor approximation of the given \code{order}.
#' 4th order is the lowest that gives self concordance.
#'
#' @return A numeric matrix with \code{(order+1)} columns containing the values of the modified log and its derivatives.
#' @export
#'
#' @examples
#' x <- seq(0.01^0.25, 2^0.25, length.out = 51)^4 - 0.11 # Denser where |f'| is higher
#' plot(x,  log(x)); abline(v = 0, lty = 2) # Observe the warning
#' lines(x, logTaylor(x, eps = 0.2), col = 2)
#' lines(x, logTaylor(x, eps = 0.5), col = 3)
#' lines(x, logTaylor(x, eps = 1, M = 1.2, order = 6), col = 4)
#'
#' # Substitute log with its Taylor approx. around 1
#' x <- seq(0.1, 2, 0.05)
#' ae <- abs(sapply(2:6, function(o) log(x) - logTaylor(x, eps=1, M=1, order=o)))
#' matplot(x[x!=1], ae[x!=1,], type = "l", log = "y", lwd = 2,
#'   main = "Abs. trunc. err. of Taylor expansion at 1", ylab = "")
logTaylor <- function(x, eps = NULL, M = NULL, der = 0, order = 4, drop = TRUE) {
  if (is.null(eps)) eps <- rep(1/length(x), length(x))
  if (is.null(M)) M <- rep(Inf, length(x))
  if (length(eps) == 1) eps <- rep(eps, length(x))
  if (length(M) == 1) M <- rep(M, length(x))
  if (length(eps) != length(M))
    stop(paste0("eps has length ", length(eps), ", but M has ", length(M), ". Make them equal!"))
  if (any(eps > M)) stop("Thresholds not ordered. eps must be less than M.")

  lo <- x < eps
  hi <- x > M
  md <- (!lo) & (!hi)
  if (is.na(order)) {
    lo <- hi <- rep(FALSE, length(x))
    md <- rep(TRUE, length(x))
  }

  out <- vapply(0:der, function(d) {
    f <- numeric(length(x))
    f[md] <- dlog(x[md], d = d)
    if (any(lo)) f[lo] <- tlog(x[lo], a = eps[lo], k = order, d = d)
    if (any(hi)) f[hi] <- tlog(x[hi], a = M[hi],   k = order, d = d)
    return(f)
  }, FUN.VALUE = numeric(length(x)))
  colnames(out) <- paste0("deriv", 0:der)
  if ((der == 0) && drop) out <- drop(out)

  return(out)
}


#' Least-squares regression via SVD
#'
#' @param X Model matrix.
#' @param y Response vector.
#' @param rel.tol Relative zero tolerance for generalised inverse via SVD.
#' @param abs.tol Relative zero tolerance for generalised inverse via SVD.
#'
#' Empirical likelihood's Newton steps are of least-squares type.
#' Denote \eqn{X^+} to be the generalised inverse of X. If SVD algorithm failures are encountered,
#' it sometimes helps to try \code{svd(t(X))} and translate back. First check to ensure that
#' \code{X} does not contain \code{NaN}, or \code{Inf}, or \code{-Inf}.
#'
#' The tolerances are used to check the closeness of singular values to zero. The values of the
#' singular-value vector \code{d} that are
#' less than \code{rel.tol * max(d) + abs.tol} are set to zero.
#'
#' @return A vector of coefficients.
#' @export
#'
#' @examples
#' b.svd <- svdlm(X = cbind(1, as.matrix(mtcars[, -1])), y = mtcars[, 1])
#' b.lm  <- coef(lm(mpg ~ ., data = mtcars))
#' b.lm - b.svd  # Negligible differences
svdlm <- function(X, y, rel.tol = 1e-9, abs.tol = 1e-100) {
  svdX <- svd(X)
  d <- svdX$d
  lo <- d < (rel.tol * max(d) + abs.tol)
  dinv <- 1 / d
  dinv[lo] <- 0
  # nrow is necessary if dinv is 1x1 to avoid getting the identity matrix of size floor(dinv)
  X.plus <- svdX$v %*% diag(dinv, nrow = length(dinv)) %*% t(svdX$u)
  ret <- X.plus %*% matrix(y, ncol = 1)
  return(drop(ret))
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
tlog <- function(x, a = as.numeric(c(1.0)), k = 4L, d = 0L)
  .Call(`_smoothemplik_tlogCPP`, x, a, k, d)
# Taken from RcppExports

