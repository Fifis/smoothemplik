#' Self-concordant empirical likelihood for a vector mean
#'
#' Implements the empirical likelihood test for whether the mean of coordinates of \code{z} is \code{mu}.
#'
#' @source This code was written and published on his personal web site by Art B. Owen (February 2014),
#'   and slightly reworked by Andreï V. Kostyrka.
#'
#' @details
#' Tweak \code{ALPHA} and \code{BETA} with extreme caution. See See Boyd and Vandenberghe, pp. 464--466 for details.
#' It is necessary that \code{0 < ALPHA < 1/2} and \code{0 < BETA < 1}.
#' \code{ALPHA = 0.03} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).
#'
#'
#' @param z A numeric vector or a matrix with one data vector per column.
#' @param mu Hypothesised mean, default (0 ... 0) in R^ncol(z)
#' @param lam Starting lambda, default (0 ... 0)
#' @param eps Lower cut-off for \code{mllog()}, default 1/nrow(z)
#' @param M Upper cutoff for \code{mllog()}, default Inf
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
#' earth <- c(5.5,5.61,4.88,5.07,5.26,5.55,5.36,5.29,5.58,5.65,5.57,5.53,5.62,5.29,
#'   5.44,5.34,5.79,5.1,5.27,5.39,5.42,5.47,5.63,5.34,5.46,5.3,5.75,5.68,5.85)
#' emplik(earth, mu = 5.517) # 5.517 is the modern accepted value
#'
#' @export
emplik <- function(z, mu, lam, eps, M, thresh = 1e-30, itermax = 100, ALPHA = 0.3, BETA = 0.8, BACKEPS = 0, verbose = FALSE) {

  if (is.vector(z)) z <- matrix(z, ncol = 1)

  n <- nrow(z)
  d <- ncol(z)

  if (missing(mu)) mu <- rep(0, d)
  z = t(t(z) - mu) # subtract mu from each z[i,]

  if (missing(eps)) eps <- 1/n
  if (missing(M)) M <- Inf

  init0 <- mllog(rep(1, n), eps = eps, M = M, der = 2) # i.e. lam = 0

  if (missing(lam)) {
    init = init0
    lam = rep(0, d)
  } else {
    init = mllog(1 + z %*% lam, eps = eps, M = M, der = 2)
    if (sum(init0[, 1]) < sum(init[, 1])) {
      lam = rep(0, d)
      init = init0
    }
  }

  # Initial f, g
  fold = sum(init[, 1])
  gold = apply(z * init[, 2], 2, sum)

  converged = FALSE
  iter = 0
  oldvals = init
  while (!converged) {
    iter = iter + 1

    # Get Newton Step
    rootllpp = sqrt(oldvals[, 3]) # sqrt 2nd deriv of -llog lik
    zt = z
    for (j in 1:d)
      zt[, j] = zt[, j] * rootllpp
    yt = oldvals[, 2] / rootllpp
    step = -svdlm(zt, yt) #  more reliable than -lm(yt ~ zt - 1)$coefficients

    backtrack = FALSE
    tt = 1 # usually called t, but R uses t for transpose
    while (!backtrack) {
      newvals = mllog(1 + z %*% (lam + tt * step), eps = eps, M = M, der = 2)
      fnew = sum(newvals[, 1])
      targ = fold + ALPHA * tt * sum(gold * step) + BACKEPS # (BACKEPS for round-off, should not be needed)
      if (fnew <= targ) {
        # backtracking has converged
        backtrack = TRUE
        oldvals = newvals
        fold = fnew
        gold = apply(z * oldvals[, 2], 2, sum)
        # take the step
        lam = lam + tt * step
      } else {
        tt = tt * BETA
      }
    }

    # Newton decrement and gradient norm
    ndec = sqrt(sum((step * gold)^2))
    gradnorm = sqrt(sum(gold^2))

    if (verbose) print(c(fold, gradnorm, ndec, lam))

    converged = (ndec^2 <= thresh)
    if (iter > itermax) break
  }

  wts <- (1 / n) / (1 + z %*% lam)
  logelr <- sum(mllog(1 + z %*% lam, eps = eps, M = M, der = 0))

  return(list(logelr = logelr, lam = lam, wts = wts, converged = converged, iter = iter, ndec = ndec, gradnorm = gradnorm))
}


#' Compute empirical likelihood on a trajectory
#'
#' @param z Passed to \code{emplik}.
#' @param mu0 Starting point of trajectory
#' @param mu1 End point of trajectory
#' @param N Number of segments into which the path is split (i. e. \code{N+1} steps are used).
#' @param verbose Logical: report iteration data?
#' @param ... Passed to \code{emplik}.
#'
#' @return A matrix with one row at each mean from mu0 to mu1 and a column for each EL return value (except EL weights).
#' @export
#'
#' @examples
tracelr <- function(z, mu0, mu1, N, verbose = FALSE, ...) {
  if (is.vector(z)) z = matrix(z, ncol = 1)
  d = ncol(z)

  ans = matrix(0, N + 1, d + 1 + d + 4)
  colnames(ans) = rep("", 2 * d + 5)
  for (j in 1:d)
    colnames(ans)[j] = paste("z", j, sep = ".")
  colnames(ans)[d + 1] = "logelr"
  for (j in 1:d)
    colnames(ans)[d + 1 + j] = paste("lambda", j, sep = ".")
  colnames(ans)[2*d + 2:5] <- c("conv", "iter", "decr", "gradnorm")

  lam = rep(0, d)
  for (i in 0:N) {
    mui = (i * mu1 + (N - i) * mu0) / N
    ansi = emplik(z = z, mu = mui, lam = lam, ...)
    ans[i + 1, ] = c(mui, ansi$logelr, ansi$lam, ansi$converged, ansi$iter, ansi$ndec, ansi$gradnorm)
    if (verbose) print(c(i, ansi$iter))
  }

  return(ans)
}

#' Modified minus logarithm with derivatives
#'
#' @param x Numeric vector for which approximated logarithm is to be computed.
#' @param eps Lower threshold below which approximation starts.
#' @param M Upper threshold above which approximation starts.
#' @param der Integer: 0, 1, or 2 derivatives
#'
#' @details
#' Provides a family of alternatives to -log() and derivative thereof in order to attain self-concordance and
#' computes the modified negative logarithm and its first derivatives.
#' For eps <= x <= M, returns just the logarithm. For x < eps and x > M, returns the 4th-order Taylor approximation.
#' 4th order is the lowest that gives self concordance.

#' @return A numeric matrix with 1, 2, or 3 columns containing the modified log and possibly its two first derivatives.
#' @export
#'
#' @examples
mllog <- function(x, eps, M = Inf, der = 0) {
  if (eps > M) stop("Thresholds not ordered. eps must be less than M.")

  lo = x < eps
  hi = x > M
  md = (!lo) & (!hi)

  # Coefficients for 4th order Taylor approx below eps
  coefs = rep(0, 5)
  coefs[1] = -log(eps)
  coefs[2:5] = (-eps)^-(1:4) / (1:4)

  # Coefficients for 4th order Taylor approx above M
  Coefs = rep(0, 5)
  Coefs[1] = -log(M)
  Coefs[2:5] = (-M)^-(1:4) / (1:4)

  # degree 4 polynomial approx to log
  h = function(y, cvals) { # cvals are coefs at eps, Coefs at M
    # sum c[t+1] y^t
    tee = 1:4
    ans = y * 0
    ans = ans + cvals[1]
    for (j in tee)
      ans = ans + y^j * cvals[j + 1]
    ans
  }

  # first derivative of h at y, from approx at pt
  hp = function(y, pt) {
    tee = 0:3
    ans = y * 0
    for (j in tee)
      ans = ans + (-y / pt)^j
    ans = ans * (-pt)^-1
    ans
  }

  # second derivative of h at y, from approx at pt
  hpp = function(y, pt) {
    tee = 0:2
    ans = y * 0
    for (j in tee)
      ans = ans + (j + 1) * (-y / pt)^j
    ans = ans * (-pt)^-2
    ans
  }

  # function value
  f = x * 0
  f[lo] = h(x[lo] - eps, coefs)
  f[hi] = h(x[hi] - M, Coefs)
  f[md] = -log(x[md])

  if (der < 1) return(cbind(f))

  # first derivative
  fp = x * 0
  fp[lo] = hp(x[lo] - eps, eps)
  fp[hi] = hp(x[hi] - M, M)
  fp[md] = -1 / x[md]

  if (der < 2) return(cbind(f, fp))

  # second derivative
  fpp = x * 0
  fpp[lo] = hpp(x[lo] - eps, eps)
  fpp[hi] = hpp(x[hi] - M, M)
  fpp[md] = 1 / x[md]^2

  return(cbind(f, fp, fpp))
}

#' Least-squares regression via SVD
#'
#' @param X Model matrix.
#' @param y Response vector.
#'
#' Empirical likelihood's Newton steps are of least-squares type.
#' Denote \eqn{X^+} to be the generalised inverse of X. If SVD algorithm failures are encountered,
#' it sometimes helps to try \code{svd(t(X))} and translate back. First check to ensure that X does not contain \code{NaN}, or \code{Inf}, or \code{-Inf}.
#'
#' @return A matrix of coefficients.
#' @export
#'
#' @examples
svdlm <- function(X, y) {
  # Tolerances for generalised inverse via SVD
  RELTOL = 1e-9
  ABSTOL = 1e-100

  svdX = svd(X)
  d = svdX$d
  lo = d < (RELTOL * max(d) + ABSTOL)
  dinv = 1 / d
  dinv[lo] = 0
  Xplus = svdX$v %*% diag(dinv, nrow = length(dinv)) %*% t(svdX$u)
  # Taking care with diag when dinv is 1x1 to avoid getting the identity matrix of size floor(dinv)

  return(Xplus %*% matrix(y, ncol = 1))
}
