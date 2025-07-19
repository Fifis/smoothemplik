# d-th derivative of the k-th-order Taylor expansion of log(x)
dlog <- function(x, d = 0L) dlogCPP(x, d)

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
#' oldpar <- par(mfrow = c(2, 3), mar = c(2, 2, 2.5, 0.2))
#' for (d in 0:5) {
#' y <- f(x, d = d)
#' plot(x, y, type = "l", lwd = 7, bty = "n", ylim = range(0, y),
#'        main = paste0("d^", d, "/dx^", d, " Taylor(Log(x))"))
#'   for (k in 0:8) lines(x, tlog(x, a = a, k = k, d = d), col = cl[k+1], lwd = 1.5)
#'   points(a, f(a, d = d), pch = 16, cex = 1.5, col = "white")
#' }
#' legend("topright", as.character(0:8), title = "Order", col = cl, lwd = 1)
#' par(oldpar)
tlog <- function(x, a = as.numeric(c(1.0)), k = 4L, d = 0L) tlogCPP(x, a, k, d)


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


#' Modified logarithm with derivatives
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


# Solve for coefficients a, b, c of the parabola y = a*x^2 + b*x + c
# using the three points (x0, y0), (x1, y1), (x2, y2) using Cramer's rule
# smoothemplik:::getParabola3(c(-1, 0, 1), c(2, 0, 2))
# x0 <- c(-1, 0, 2); y0 <- c(0, 1, 1)
# abc <- smoothemplik:::getParabola3(x0, y0)
# f <- function(x) abc[1]*x^2 + abc[2]*x + abc[3]
# curve(f, -1.5, 2.5); points(x0, y0)
getParabola3 <- function(x, y) getParabola3CPP(x, y)


# Get the coefficients of a parabola with f(x), f'(x), f''
# f(x) = 2x^2 + 3x - 5
# f' = 4x + 3
# f'' = 4
# If x = -1, (f, f', f'') = (-6, -1, 4)
# smoothemplik:::getParabola(-1, -6, -1, 4)
getParabola <- function(x, f, fp, fpp) getParabolaCPP(x, f, fp, fpp)


#' Monotone interpolation between a function and a reference parabola
#'
#' Create *piece-wise monotone* splines that smoothly join an arbitrary function `f` to the
#' quadratic reference curve \eqn{(x-\mathrm{mean})^{2}/\mathrm{var}} at
#' a user–chosen abscissa \code{at}. The join occurs over a finite interval
#' of length \code{gap}, guaranteeing a C1-continuous transition (function and first
#' derivative are continuous) without violating monotonicity.
#'
#' @param x A numeric vector of evaluation points.
#' @param f Function: the original curve to be spliced into the parabola.
#'   It must be vectorised (i.e.\ accept a numeric vector and return a numeric
#'   vector of the same length).
#' @param mean Numeric scalar defining the shift of the reference parabola.
#' @param var Numeric scalar defining the vertical scaling of the reference parabola.
#' @param at Numeric scalar: beginning of the transition zone, i.e.\ the
#'   boundary where `f` stops being evaluated and merging into the parabola begins.
#' @param gap Positive numeric scalar.  Width of the transition window; the spline
#'   is constructed on `[at, at+gap]` (or `[at-gap, at]` when `at < mean`) when the
#'   reference parabola is higher. If the reference parabola is lower, it is the
#'   distance from the point `z` at which `f(z) = parabola(z)` to allow some growth and
#'   ensure monotonicity.
#'
#' @details
#' This function calls `interpToHigher()` when the reference parabola is *above* `f(at)`;
#' the spline climbs from `f` up to the parabola, and `interpToLower()` when the parabola is *below*
#' `f(at)`, and the transition interval has to be extended to ensure that the spline does not descend.
#'
#' Internally, the helpers build a **monotone Hermite cubic spline** via
#' Fritsch--Carlson tangents.  Anchor points on each side of the
#' transition window are chosen so that the spline’s one edge matches `f`
#' while the other edge matches the reference parabola, ensuring strict
#' monotonicity between the two curves.
#'
#' @returns A numeric vector of length \code{length(x)} containing the smoothly
#'   interpolated values.
#'
#' @seealso [splinefun()]
#' @export
#'
#' @examples
#' xx <- -4:5  # Global data for EL evaluation
#' w <- 10:1
#' w <- w / sum(w)
#'
#' f <- Vectorize(function(m) -2*EL0(xx, mu = m, ct = w, chull.fail = "none")$logelr)
#' museq <- seq(-6, 6, 0.1)
#' LRseq <- f(museq)
#' plot(museq, LRseq, bty = "n")
#' rug(xx, lwd = 4)
#'
#' wm <- weighted.mean(xx, w)
#' wv <- weighted.mean((xx-wm)^2, w) / sum(w)
#' lines(museq, (museq - wm)^2 / wv, col = 2, lty = 2)
#'
#' xr <- seq(4, 6, 0.1)
#' xl <- seq(-6, -3, 0.1)
#' lines(xl, interpTwo(xl, f, mean = wm, var = wv, at = -3.5, gap = 0.5), lwd = 2, col = 4)
#' lines(xr, interpTwo(xr, f, mean = wm, var = wv, at = 4.5, gap = 0.5), lwd = 2, col = 3)
#' abline(v = c(-3.5, -4, 4.5, 5), lty = 3)
interpTwo <- function(x, f, mean, var, at, gap) {
  fa <- function(x) (x-mean)^2/var
  if (f(at) <= fa(at))
    interpToHigherCPP(x, f, mean, var, at, gap)
  else
    interpToLowerCPP(x, f, mean, var, at, gap)
}


#' Damped Newton optimiser
#'
#' @param fn A function that returns a list: f, f', f''.
#' If the function takes vector arguments, the dimensions of the list components must
#' be 1, \code{dim X}, \code{(dim X) x (dim X)}. The function must be (must be twice continuously
#' differentiable at x)
#' @param par Numeric vector: starting point.
#' @param thresh A small scalar: stop when Newton decrement squared falls belowe \code{thresh}.
#' @param itermax Maximum iterations. Consider optimisation failed if the maximum is reached.
#' @param alpha Back-tracking parameter strictly between 0 and 0.5: acceptance of a decrease in function value by alpha*f of the prediction.
#' @param beta Back-tracking parameter strictly between 0 and 1: reduction of the step size until
#'   the stopping criterion is met. 0.1 corresponds to a very crude search, 0.8 corresponds
#'   to a less crude search.
#' @param backeps Back-tracking threshold: the search can miss by this much. Consider setting it to 1e-10
#'   if backtracking seems to be failing due to round-off.
#' @param verbose Logical: if true, prints the tracing infornation (iteration log).
#'
#' This is a translation of Algorithm 9.5 from \insertCite{boyd2004convex}{smoothemplik} into C++.
#'
#' @references
#' \insertAllCited{}
#'
#' @return A list:
#' @export
#'
#' @examples
#' f1 <- function(x)
#'   list(fn = x - log(x), gradient = 1 - 1/x, Hessian = matrix(1/x^2, 1, 1))
#' optim(2, function(x) f1(x)[["fn"]], gr = function(x) f1(x)[["gradient"]], method = "BFGS")
#' dampedNewton(f1, 2, verbose = TRUE)
#'
#' # The minimum of f3 should be roughly at -0.57
#' f3 <- function(x)
#'   list(fn = sum(exp(x) + 0.5 * x^2), gradient = exp(x) + x, Hessian =  diag(exp(x) + 1))
#' dampedNewton(f3, seq(0.1, 5, length.out = 11), verbose = TRUE)
dampedNewton <- function(fn, par, thresh  = 1e-30, itermax = 100,
                         verbose = FALSE, alpha = 0.3, beta = 0.8, backeps = 0.0) {
  dampedNewtonCPP(fn = fn, par = par, thresh = thresh, itermax = itermax, verbose = verbose,
                  alpha = alpha, beta = beta, backeps = backeps)
}

# Bartlett correction factor for the 1-dimensional Adjusted EL
computeBartlett <- function(x, ct = NULL) {
  # Expressions from Liu & Chen (2011, Annals of Statistic, p.1347)
  n <- length(x)
  if (is.null(ct)) ct <- rep(1, n)
  alpha1 <- stats::weighted.mean(x, ct)
  zc <- x - alpha1
  alpha2hat <- stats::weighted.mean(zc^2, ct)
  alpha3hat <- stats::weighted.mean(zc^3, ct)
  alpha4hat <- stats::weighted.mean(zc^4, ct)
  alpha6hat <- stats::weighted.mean(zc^6, ct)
  alpha2 <- n/(n-1) * alpha2hat
  alpha4 <- n/(n-4)*alpha4hat - 6/(n-4)*alpha2^2
  alpha22 <- alpha2^2 - alpha4/n
  alpha3 <- n/(n-3)*alpha3hat
  alpha33 <- alpha3^2 - (alpha6hat - alpha3^2)/n
  alpha222 <- alpha2^3
  b <- 0.5*alpha4/alpha22 - alpha33/alpha222/3
  return(b)
}


#' Weighted trimmed mean
#'
#' Compute a weighted trimmed mean, i.e. a mean that assigns non-negative weights
#' to the observations and (2)  discards an equal share of total weight from
#' each tail of the distribution before averaging.
#'
#' For example, `trim = 0.10` removes 10% of the weight from the left tail and 10%
#' from the right (20% in total), then takes the weighted mean of what is left.
#' Setting `trim = 0.5` returns the weighted median.
#'
#' @param x Numeric vector of data values.
#' @param trim Single number in \eqn{[0,\,0.5]}. Fraction of the total weight to cut
#'   from each tail.
#' @param w Numeric vector of non-negative weights of the same length as `x`.
#'   If `NULL` (default), equal weights are used.
#' @param na.rm Logical: should `NA` values in `x` or `w` be removed?
#' @param ... Further arguments passed to [`weighted.mean()`] (for compatibility).
#'
#' @details
#' The algorithm follows these steps:
#' \enumerate{
#'   \item Sort the data by `x` and accumulate the corresponding weights.
#'   \item Identify the lower and upper cut-points that mark the central
#'         share of the total weight.
#'   \item Drop observations whose cumulative weight lies entirely
#'         outside the cut-points and proportionally down-weight the two (at most)
#'         remaining outermost observations.
#'   \item Return the weighted mean of the retained mass.  If `trim == 0.5`,
#'         only the 50% quantile remains, so the function returns the weighted median.
#' }
#'
#' @return A single numeric value: the trimmed weighted mean of `x`. Returns `NA_real_`
#' if no non-`NA` observations remain after optional `na.rm` handling.
#' @export
#'
#' @seealso
#' [`mean()`] for the unweighted trimmed mean, [`weighted.mean()`] for the untrimmed weighted mean.
#'
#' @examples
#' set.seed(1)
#' z <- rt(100, df = 3)
#' w <- pmin(1, 1 / abs(z)^2)  # Far-away observations tails get lower weight
#'
#' mean(z, trim = 0.20)  # Ordinary trimmed mean
#' trimmed.weighted.mean(z, trim = 0.20)  # Same
#'
#' weighted.mean(z, w)   # Ordinary weighted mean (no trimming)
#' trimmed.weighted.mean(z, w = w)  # Same
#'
#' trimmed.weighted.mean(z, trim = 0.20, w = w)  # Weighted trimmed mean
#' trimmed.weighted.mean(z, trim = 0.5,  w = w)  # Weighted median
trimmed.weighted.mean <- function(x, trim = 0, w = NULL, na.rm = FALSE, ...) {
  if (trim == 0)  return(stats::weighted.mean(x, w, na.rm = na.rm, ...))
  if (is.null(w)) return(mean(x, trim = trim, na.rm = na.rm, ...))
  if (length(w) != length(x)) stop("'w' must have the same length as 'x'")
  if (any(w < 0, na.rm = TRUE)) stop("'w' must contain non-negative numbers only")
  if (trim < 0) stop("'trim' must be in [0, .5]")
  if (trim > 0.5) trim <- 0.5

  if (na.rm) {
    ok <- !(is.na(x) | is.na(w))
    x  <- x[ok]
    w <- w[ok]
  }
  if (length(x) < 1) return(NA_real_)

  # Sort all observations by value
  ord <- order(x)
  x   <- x[ord]
  w   <- w[ord]

  totw <- sum(w)
  cumwfrac <- cumsum(w) / totw  # Cumulative weight share
  prevfrac <- c(0, utils::head(cumwfrac, -1))  # Share up to before obs
  # Lower and upper cut points (fractions of weight)
  L <- trim
  U <- 1 - trim

  # Initialise kept weights after trimming
  keepw <- w
  keepw[cumwfrac <= L | prevfrac >= U] <- 0  # Fully trimmed observations
  ## Partially trimmed observations (one on each side at most)
  idxL <- prevfrac < L & cumwfrac > L
  if (any(idxL)) keepw[idxL] <- w[idxL] * (cumwfrac[idxL] - L) / (w[idxL]/totw)
  idxU <- prevfrac < U & cumwfrac > U
  if (any(idxU)) keepw[idxU] <- w[idxU] * (U - prevfrac[idxU]) / (w[idxU]/totw)

  stats::weighted.mean(x, keepw, na.rm = na.rm, ...)
}

