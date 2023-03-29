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
#' @param eps Lower cut-off for \code{mllog()}, default 1/nrow(z)
#' @param M Upper cutoff for \code{mllog()}, default Inf
#' @param order Positive integer such that the Taylor approximation of this order to log(x) is self-concordant; usually 4 or higher. Passed to \code{mllog}.
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
#'     }
#'
#' @examples
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' cemplik(earth, mu = 5.517, verbose = TRUE) # 5.517 is the modern accepted value
#'
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
  init <- mllog(rep(1, n), eps = eps, M = M, order = order, der = 2) # i.e. fn and derivatives at lam = 0
  init <- init * ct # sweep(init, 1, ct, "*") # Weights appear OUTSIDE of mllog
  # If the initial lambda is supplied, compare lam = 0 with it and pick the best
  if (!is.null(lam)) {
    init1 <- mllog(1 + z %*% lam, eps = eps, M = M, order = order, der = 2)
    init1 <- init1 * ct # sweep(init1, 1, ct, "*")
    if (sum(init1[, 1]) < sum(init[, 1])) { # Minimising --> lower is better
      init <- init1
    } else lam <- rep(0, d) # Ignore the supplied lambda
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
    zt <- z * rootllpp # sweep(z, 1, rootllpp, "*")
    yt <- oldvals[, 2] / rootllpp
    step <- -svdlm(zt[nonz, ], yt[nonz]) #  SVD is more reliable than step <- -lm( yt~zt-1 )$coef
    # Also, we kick out the ones with ct=0 to avoid unnecessary NaNs
    backtrack <- FALSE
    tt <- 1 # usually called t, but R uses t for transpose
    while (!backtrack) {
      newlam <- lam + tt*step
      newvals <- mllog(1 + z %*% newlam, eps = eps, M = M, order = order, der = 2)
      newvals <- newvals * ct # sweep(newvals, 1, ct, "*")
      fnew <- sum(newvals[, 1])
      targ <- fold + ALPHA * tt * sum(gold*step) + BACKEPS # (BACKEPS for round-off, should not be needed)
      if (fnew <= targ) { # Backtracking has converged
        backtrack <- TRUE
        oldvals <- newvals
        fold <- fnew
        gold <- apply(z*oldvals[, 2], 2, sum)
        # take the step
        lam <- newlam
      } else tt <- tt * BETA
    }

    # Newton decrement and gradient norm
    ndec <- sqrt(sum((step*gold)^2))
    gradnorm <- sqrt(sum(gold^2))

    if (verbose) cat("Iter. ", iter, ", lambda = ", paste0(formatC(lam, digits = 8), collapse = " "),
                     ", f = ", fold, ", |grad| = ", gradnorm, ", decr = ", ndec, "\n", sep = "")
    converged <- (ndec^2 <= thresh)
    if (iter > itermax) break
  }

  wts <- (ct / sum(ct)) / (1 + z %*% lam) # Without the weights, the numerator was 1/n
  logelr <- sum(ct * mllog(1 + z %*% lam, eps = eps, M = M, order = order, der = 0))

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
#' @return A matrix with one row at each mean from mu0 to mu1 and a column for each EL return value (except EL weights).
#' @export
#'
#' @examples
ctracelr <- function(z, ct, mu0, mu1, N = 5, verbose = FALSE, ...) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  d <- ncol(z)

  ans <- matrix(0, N + 1, d + 1 + d + 4)
  colnames(ans)[1:d] = paste("z", 1:d, sep = ".")
  colnames(ans)[d+1] = "logelr"
  colnames(ans)[d+1 + 1:d] = paste("lambda", 1:d, sep = ".")
  colnames(ans)[2*d+ 2:5] <- c("conv", "iter", "decr", "gnorm")

  for (i in 0:N) {
    mui <- (i*mu1 + (N-i)*mu0) / N
    ansi <- cemplik(z = z, ct = ct, mu = mui, ...)
    ans[i + 1, ] = c(mui, ansi$logelr, ansi$lam, ansi$converged, ansi$iter, ansi$ndec, ansi$gradnorm)
    if (verbose) cat("Point ", i, "/", N, ", ", if (ansi$converged == 0) "NOT " else "", "converged",
                     ", log(ELR) = ", ansi$logelr, "\n", sep = "")
  }
  return(ans)
}

#' Modified minus logarithm with derivatives
#'
#' @param x Numeric vector for which approximated logarithm is to be computed.
#' @param order Positive integer: Taylor approximation order.
#' @param eps Lower threshold below which approximation starts.
#' @param M Upper threshold above which approximation starts.
#' @param der Non-negative integer: 0 yields the function, 1 and higher yields derivatives
#' @param flatten1d If TRUE and \code{x} is a vector, the output is a vector and not a 1-column matrix.
#'
#' @details
#' Provides a family of alternatives to -log() and derivative thereof in order to attain self-concordance and
#' computes the modified negative logarithm and its first derivatives.
#' For eps <= x <= M, returns just the logarithm. For x < eps and x > M, returns the Taylor approximation of the given \code{order}.
#' 4th order is the lowest that gives self concordance.

#' @return A numeric matrix with \code{(order+1)} columns containing the values of the modified log and its derivatives.
#' @export
#'
#' @examples
#' x <- seq(0.01^0.25, 2^0.25, length.out = 51)^4 - 0.11 # Denser where |f'| is higher
#' plot(x, -log(x))
#' lines(x, mllog(x, eps = 0.2), col = 2)
#' lines(x, mllog(x, eps = 0.5), col = 3)
#' lines(x, mllog(x, eps = 1, M = 1.2, order = 6), col = 4)
#' # Suppose that we substitute -log with its Taylor approx. around 1
#' x <- seq(0.1, 1.5, 0.2)
#' err.mat <- cbind(x, ErrOrder2 = -log(x) - mllog(x, eps = 1, M = 1, order = 2),
#'                     ErrOrder3 = -log(x) - mllog(x, eps = 1, M = 1, order = 3),
#'                     ErrOrder4 = -log(x) - mllog(x, eps = 1, M = 1, order = 4),
#'                     ErrOrder5 = -log(x) - mllog(x, eps = 1, M = 1, order = 5),
#'                     ErrOrder6 = -log(x) - mllog(x, eps = 1, M = 1, order = 6))
#' print(err.mat, digits = 2)
mllog <- function(x, eps = NULL, M = Inf, der = 0, order = 4, flatten1d = TRUE) {
  if (is.null(eps)) eps <- 1/length(x)
  if (eps > M) stop("Thresholds not ordered. eps must be less than M.")

  lo <- x < eps
  hi <- x > M
  md <- (!lo) & (!hi)

  dlog <- function(x, r = 0)  if (r == 0) log(x) else ((r%%2 == 1)*2-1) * 1/x^r * gamma(r)
  out <- sapply(0:der, function(r) {
    f <- x*0
    f[md] <- dlog(x[md], r = r)
    if (any(lo)) f[lo] <- logTaylor(x[lo], a = eps, k = order, r = r)
    if (any(hi)) f[hi] <- logTaylor(x[hi], a = M,   k = order, r = r)
    return(f)
  })
  if ((der == 0) & flatten1d) out <- as.vector(out)

  return(-out)
}

#' Least-squares regression via SVD
#'
#' @param X Model matrix.
#' @param y Response vector.
#'
#' Empirical likelihood's Newton steps are of least-squares type.
#' Denote \eqn{X^+} to be the generalised inverse of X. If SVD algorithm failures are encountered,
#' it sometimes helps to try \code{svd(t(X))} and translate back. First check to ensure that \code{X} does not contain \code{NaN}, or \code{Inf}, or \code{-Inf}.
#'
#' @return A matrix of coefficients.
#' @export
#'
#' @examples
svdlm <- function(X, y) {
  # Tolerances for generalised inverse via SVD
  RELTOL <- 1e-9
  ABSTOL <- 1e-100

  svdX <- svd(X)
  d <- svdX$d
  lo <- d < (RELTOL * max(d) + ABSTOL)
  dinv <- 1 / d
  dinv[lo] <- 0
  # nrow is necessary if dinv is 1x1 to avoid getting the identity matrix of size floor(dinv)
  Xplus <- svdX$v %*% diag(dinv, nrow = length(dinv)) %*% t(svdX$u)

  return(Xplus %*% matrix(y, ncol = 1))
}

#' r-th derivative of the kth-order Taylor expansion of log(x)
#'
#' @param x Numeric: a vector of points for which the logarithm is to be evaluated
#' @param a Scalar: the point at which the polynomial approximation is computed
#' @param k Positive integer: maximum polynomial order in the Taylor expansion of the original function.
#' @param r Non-negative integer: derivative order
#'
#' Note that this function returns the r-th derivative of the k-th-order Taylor expansion, not the
#' k-th-order approximation of the r-th derivative. Therefore, the degree of the resulting polynomial
#' is \eqn{r-k}.
#'
#' @return The approximating Taylor polynomial around \code{a} of the order \code{r-k} evaluated at \code{x}.
#' @export
#'
#' @examples
#' cl <- c("#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
#'         "#911eb4", "#42d4f4", "#f032e6", "#bfef45")
#' a <- 1.5
#' x <- seq(a*2, a/2, length.out = 101)
#' f <- function(x, r = 0)  if (r == 0) log(x) else ((r%%2 == 1)*2-1) * 1/x^r * gamma(r)
#' par(mfrow = c(2, 3), mar = c(2, 2, 2.5, 0.2))
#' for (r in 0:5) {
#' y <- f(x, r = r)
#' plot(x, y, type = "l", lwd = 7, bty = "n", ylim = range(0, y),
#'        main = paste0("d^", r, "/dx^", r, " Taylor(Log(x))"))
#'   for (k in 0:8) lines(x, logTaylor(x, a = a, k = k, r = r), col = cl[k+1], lwd = 1)
#'   points(a, f(a, r = r), pch = 16, cex = 1.5, col = "white")
#' }
#' legend("topright", as.character(0:8), title = "Order", col = cl, lwd = 1)
logTaylor <- function(x, a = 1, k = 4, r = 0) {
  if (r > k) return(x*0) # Polynomial derivatives of order > k are zero
  l <- length(x)
  xc <- (x-a)/a
  taylor <- sapply(r:k, function(n) { # Terms of the Taylor expansion
    if (n == r) { # Lowest order: constant
      if (r == 0) return(rep(log(a), l)) # Original function = the only special term
      return(rep(-1/(-a)^r * gamma(n), l))
    }
    mult <- if (r == 0) 1 else prod(n:(n-r+1)) # Multiplier from the polynomial power
    sgn <- if (n%%2 == 1) 1 else -1
    mult * sgn * xc^(n-r) / n / a^r
  })
  if (l == 1) taylor <- matrix(taylor, nrow = 1) # If sapply is not a matrix
  rowSums(taylor)
}

