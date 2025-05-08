# Get the coefficients of a parabola passing through (x, y)
# getParabola3(c(-1, 0, 1), c(2, 0, 2))
# x0 <- c(-1, 0, 2); y0 <- c(0, 1, 1)
# abc <- getParabola3(x0, y0)
# f <- function(x) abc[1]*x^2 + abc[2]*x + abc[3]
# curve(f, -1.5, 2.5); points(x0, y0)
getParabola3 <- function(x, y) {
  if (length(x) != 3 || length(y) != 3) stop("Exactly 3 points are required to fit a parabola.")
  xprod <- (x[1] - x[2])*(x[1] - x[3])*(x[2] - x[3])
  a <- (x[1] * (y[3] - y[2]) + x[2] * (y[1] - y[3]) + x[3] * (y[2] - y[1])) / xprod
  b <- (x[1]^2 * (y[2] - y[3]) + x[2]^2 * (y[3] - y[1]) + x[3]^2 * (y[1] - y[2])) / xprod
  c <- (x[2] * (x[1]^2 * y[3] - x[3]^2 * y[1]) + x[2]^2 * (x[3] * y[1] - x[1] * y[3]) + x[1] * x[3] * y[2] * (x[3] - x[1])) / xprod
  c(a = a, b = b, c = c)
}

# Get the coefficients of a parabola with f(x), f'(x), f''
# f(x) = 2x^2 + 3x - 5
# f' = 4x + 3
# f'' = 4
# If x = -1, (f, f', f'') = (-6, -1, 4)
# getParabola(-1, -6, -1, 4)
getParabola <- function(x, f, fp, fpp) {
  a <- fpp/2
  b <- fp - fpp * x
  c <- f + 0.5 * x * (fpp * x - 2*fp)
  return(c(a = a, b = b, c = c))
}


#' Empirical likelihood for one-dimensional vectors
#'
#' Empirical likelihood with counts to solve one-dimensional problems efficiently with Brent's root search algorithm.
#' Conducts an empirical likelihood ratio test of the hypothesis that the mean of \code{z} is \code{mu}.
#' The names of the elements in the returned list are consistent with the original R code in \insertCite{owen2017weighted}{smoothemplik}.
#'
#' @param z Numeric data vector.
#' @param ct Numeric count variable with positive values that indicates the multiplicity of observations.
#'   Can be fractional. Very small counts below the threshold \code{weight.tolerance} are zeroed.
#' @param shift The value to add in the denominator (useful in case there are extra Lagrange multipliers): 1 + lambda'Z + shift.
#' @param taylor.order If \code{NA}, then the sum of log-weights is maximised. If 2, 4, 6, or 8,
#'   replaces the left branch of the log with its Taylor expansion, allowing for \code{mu}
#'   to be outside the convex hull of \code{z}.
#' @param lower A positive scalar or numeric vector denoting the left cut-off before which log(weight)
#'   in the maximisation problem is replaced by its Taylor expansion if \code{taylor.order} is 2, 4, 6, or 8.
#'   If \code{NULL} or \code{-Inf}, no such replacement is done. Recommended at a small number ct/sum(ct)
#'   to ensure self-concordance of the log-likelihood.
#' @param upper A positive scalar or numeric vector greater or equal to \code{lower} denoting the right cut-off
#'   after which log(weight) in the maximisation problem is replaced by its Taylor expansion if \code{taylor.order}
#'   is 2, 4, 6, or 8. If \code{NULL} or \code{Inf}, no such replacement is done. Recommended at 1
#'   if \code{mu} is not in the convex hull of \code{z} to ensure that a finite solution exists.
#' @param mu Hypothesized mean of \code{z} in the moment condition.
#' @param SEL If \code{FALSE}, then the boundaries for the lambda search are based on the total sum of counts, like in vanilla empirical likelihood,
#' due to formula (2.9) in \insertCite{owen2001empirical}{smoothemplik}, otherwise according to Cosma et al. (2019, p. 170, the topmost formula).
#' @param n.orig An optional scalar to denote the original sample size (useful in the rare cases re-normalisation is needed).
#' @param weight.tolerance Weight tolerance for counts to improve numerical stability
#'   (similar to the ones in Art B. Owen's 2017 code, but adapting to the sample size).
#' @param chull.fail A character: what to do if the convex hull of \code{z} does not contain \code{mu}
#'   (spanning condition does not hold). \code{"taylor"} requests a Taylor approximation
#'   of the logarithm outside \code{[lower, upper]} or, if any of those is not provided,
#'   \code{[ct/sum(ct), 1]}. \code{"wald"} eliminates the numerical lambda search and replaces the LR statistic
#'   with the Wald statistic (testing if the weighted mean of \code{z} is zero).
#' @param boundary.tolerance Relative tolerance for determining when the lambda is not an interior
#'   solution because it is too close to the boundary. Corresponds to a fraction of the
#'   interval range length.
#' @param trunc.to Counts under \code{weight.tolerance} will be set to this value.
#'   In most cases, setting this to \code{0} or \code{weight.tolerance} is a viable solution of the zero-denominator problem.
#' @param uniroot.control A list passed to the \code{uniroot}.
#' @param return.weights Logical: if TRUE, individual EL weights are computed and returned.
#'   Setting this to FALSE gives huge memory savings in large data sets, especially when smoothing is used.
#' @param verbose Logical: if \code{TRUE}, prints warnings.
#'
#' @details
#' This function provides the core functionality for univariate empirical likelihood.
#' The technical details is given in \insertCite{cosma2019inference}{smoothemplik},
#' although the algorithm used in that paper is slower than the one provided by this function.
#'
#' Since we know that the EL probabilities belong to (0, 1), the interval (bracket) for \eqn{\lambda}{l} search
#' can be determined in the spirit of formula (2.9) from \insertCite{owen2001empirical}{smoothemplik}. Let
#' \eqn{z^*_i := z_i - \mu}{Z[i] := z[i] - mu} be the recentred observations.
#' \deqn{p_i = c_i/N \cdot (1 + \lambda z^*_i + s)^{-1}}{p[i] = c[i]/N * 1/(1 + l*Z[i] + s)}
#' The probabilities are bounded from above: \eqn{p_i<1}{p[i] < 1} for all \emph{i}, therefore,
#' \deqn{c_i/N \cdot (1 + \lambda z^*_i + s)^{-1} < 1}{c[i]/N * 1/(1 + l*Z[i] + s) < 1}
#' \deqn{c_i/N - 1 - s < \lambda z^*_i}{c[i]/N - 1 - s < l*Z[i]}
#' Two cases: either \eqn{z^*_i<0}{Z[i] < 0}, or \eqn{z^*_i>0}{Z[i] > 0}
#' (cases with \eqn{z^*_i=0}{Z[i] = 0} are trivially excluded because they do not affect the EL). Then,
#' \deqn{(c_i/N - 1 - s)/z^*_i > \lambda,\ \forall i: z^*_i<0}{(c[i]/N - 1 - s)/Z[i] > l,  such i that Z[i]<0}
#' \deqn{(c_i/N - 1 - s)/z^*_i < \lambda,\ \forall i: z^*_i>0}{(c[i]/N - 1 - s)/Z[i] < l,  such i that Z[i]>0}
#' which defines the search bracket:
#' \deqn{\lambda_{\min} := \max_{i: z^*_i>0} (c_i/N - 1 - s)/z^*_i}{l > max_{i: Z[i]>0} (c_i/N - 1 - s)/Z[i]}
#' \deqn{\lambda_{\max} := \min_{i: z^*_i<0} (c_i/N - 1 - s)/z^*_i}{l < min_{i: Z[i]<0} (c_i/N - 1 - s)/Z[i]}
#' \deqn{\lambda_{\min} < \lambda < \lambda_{\max}}
#'
#' (This derivation contains \emph{s}, which is the extra shift that extends the
#' function to allow mixed conditional and unconditional estimation;
#' Owen's textbook formula corresponds to \eqn{s = 0}{s = 0}.)
#'
#' The actual tolerance of the lambda search in \code{uniroot} is
#' \eqn{2 |\lambda_{\max}| \epsilon_m + \mathrm{tol}/2}{2 * MachEps * l_max + tol/2},
#' where \code{tol} can be set in \code{uniroot.control} and
#' \eqn{\epsilon_m}{MachEps} is \code{.Machine$double.eps}.
#'
#' TODO: write a section of the concondant minus-log expansion
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{logelr}{Logarithm of the empirical likelihood ratio.}
#'   \item{lam}{The Lagrange multiplier.}
#'   \item{wts}{Observation weights/probabilities (of the same length as \code{z}).}
#'   \item{converged}{\code{TRUE} if the algorithm converged, \code{FALSE} otherwise (usually means that \code{mu} is not within the range of \code{z}, i.e. the one-dimensional convex hull of \code{z}).}
#'   \item{iter}{The number of iterations used (from \code{uniroot}).}
#'   \item{bracket}{The admissible interval for lambda (that is, yielding weights between 0 and 1).}
#'   \item{estim.prec}{The approximate estimated precision of lambda (from \code{uniroot}).}
#'   \item{f.root}{The value of the derivative of the objective function w.r.t. lambda at the root (from \code{uniroot}). Values \code{> sqrt(.Machine$double.eps)} indicate convergence problems.}
#'   \item{exitcode}{An integer indicating the reason of termination.}
#'   \item{message}{Character string describing the optimisation termination status.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Figure 2.4 from Owen (2001) -- with a slightly different data point
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' set.seed(1)
#' system.time(r1 <- replicate(40, cemplik(sample(earth, replace = TRUE), mu = 5.517)))
#' set.seed(1)
#' system.time(r2 <- replicate(40, weightedEL(sample(earth, replace = TRUE), mu = 5.517)))
#' plot(apply(r1, 2, "[[", "logelr"), apply(r1, 2, "[[", "logelr") - apply(r2, 2, "[[", "logelr"),
#'      bty = "n", xlab = "log(ELR) computed via dampened Newthon method",
#'      main = "Discrepancy between cemplik and weightedEL", ylab = "")
#' abline(h = 0, lty = 2)
#'
#' # Handling the convex hull violation differently
#' weightedEL(1:9, chull.fail = "none")
#' weightedEL(1:9, chull.fail = "taylor")
#' weightedEL(1:9, chull.fail = "wald")
#'
#' # Interpolation to well-defined branches outside the convex hull
#' mu.seq <- seq(-1, 7, 0.1)
#' wEL1 <- -2*sapply(mu.seq, function(m) weightedEL(1:9, mu = m, chull.fail = "none")$logelr)
#' wEL2 <- -2*sapply(mu.seq, function(m) weightedEL(1:9, mu = m, chull.fail = "taylor")$logelr)
#' wEL3 <- -2*sapply(mu.seq, function(m) weightedEL(1:9, mu = m, chull.fail = "wald")$logelr)
#' plot(mu.seq, wEL1)
#' lines(mu.seq, wEL2, col = 2)
#' lines(mu.seq, wEL3, col = 4)
#'
#' # Warning: depending on the compiler, the discrepancy between cemplik and weightedEL
#' # can be one million (1) times larger than the machine epsilon despite both of them
#' # being written in pure R
#' # The results from Apple clang-1400.0.29.202 and Fortran GCC 12.2.0 are different from
#' # those obtained under Ubuntu 22.04.4 + GCC 11.4.0-1ubuntu1~22.04,
#' # Arch Linux 6.6.21 + GCC 14.1.1, and Windows Server 2022 + GCC 13.2.0
#' out1 <- cemplik(earth, mu = 5.517)[1:4]
#' out2 <- weightedEL(earth, mu = 5.517, return.weights = TRUE)[1:4]
#' print(c(out1$lam, out2$lam), 16)
#'
#' # Value of lambda                            cemplik           weightedEL
#' # aarch64-apple-darwin20         -1.5631313955??????   -1.5631313957?????
#' # Windows, Ubuntu, Arch           -1.563131395492627   -1.563131395492627
#' @export
weightedEL <- function(z, mu = 0, ct = NULL, shift = NULL,
                       taylor.order = NA, lower = NULL, upper = NULL,
                       SEL = FALSE,
                       n.orig = NULL,
                       weight.tolerance = NULL, boundary.tolerance = 1e-9,
                       trunc.to = 0, chull.fail = c("taylor", "wald", "none"),
                       uniroot.control = list(),
                       return.weights = FALSE,
                       verbose = FALSE
) {
  me <- .Machine$double.eps
  if (is.data.frame(z)) z <- as.matrix(z, rownames.force = TRUE)
  if (any(!is.finite(z))) stop("Non-finite observations (NA, NaN, Inf) are not welcome.")
  if (is.null(ct)) ct <- rep(1, length(z))
  if (any(!is.finite(ct))) stop("Non-finite weights (NA, NaN, Inf) are not welcome.")
  if (min(ct) < 0 && verbose) warning("Negative weights are present.")
  if (sum(ct) <= 0) stop("The total sum of EL weights must be positive.")
  if (is.null(n.orig)) n.orig <- length(z)
  if (is.null(shift)) shift <- rep(0, length(z))
  if (is.finite(taylor.order) && (!(taylor.order %in% (1:4)*2)))
    stop("'taylor.order' must be 2, 4, 6, or 8.")
  if (is.null(weight.tolerance)) weight.tolerance <- if (!SEL) me^(1/3) else max(ct) * sqrt(.Machine$double.eps)
  chull.fail <- match.arg(chull.fail)

  # If originally the weights were too small, too many points would be truncated
  # Warn if any non-zero weights are smaller than weight.tolerance
  closeto0 <- (abs(ct) < weight.tolerance)
  if (any(closeto0)) {
    if (verbose)
      warning(paste0("Counts closer to 0 than ", sprintf("%1.2e", weight.tolerance),
                     " have been replaced with ",
                     if (trunc.to == 0) "0." else  paste0("+-", trunc.to, " of appropriate sign.")))
    ct[closeto0 & ct > 0] <- trunc.to
    ct[closeto0 & ct < 0] <- -trunc.to
  }

  if (is.matrix(z)) {
    if (ncol(z) > 1) stop("Only one-dimensional vectors or matrices are supported.") else z <- drop(z)
  }

  if (return.weights) {
    wts <- numeric(n.orig)
    names(wts) <- names(z)
  } else {
    wts <- NULL
  }

  # Not all observations contribute meaningfully to the likelihood; tiny weights
  # push lambda to the boundary
  nonz <- which(ct != 0)
  n.init <- length(z)
  n.final <- length(nonz)
  if (n.final < n.init) {
    z <- z[nonz]
    ct <- ct[nonz]
    shift <- shift[nonz]
    if (length(lower) == n.init) lower <- lower[nonz]
    if (length(upper) == n.init) upper <- upper[nonz]
  }
  if (SEL) ct <- ct / sum(ct) # We might have truncated some weights, so re-normalisation is needed!
  # The denominator for EL with counts is the sum of total counts, and for SEL, it is the number of observations!
  N <- if (SEL) n.orig else sum(ct)
  if (N <= 0) stop(paste0("Total weights after tolerance checks (", N, ") must be positive (check the counts and maybe decrease 'weight.tolerance', which is now ", sprintf("%1.1e", weight.tolerance), ".\n"))

  z <- z - mu

  z1 <- min(z)
  zn <- max(z)

  # The default output is for the worst case when the spanning condition is not met
  # It is overwritten if optimisation succeeds
  lam <- Inf
  logelr <- -Inf
  converged <- FALSE
  iter <- NA
  int <- c(-Inf, Inf)
  estim.prec <- NA
  f.root <- NA
  exitcode <- 5

  if (is.null(lower)) lower <- rep(-Inf, n.final)
  if (is.null(upper)) upper <- rep(Inf, n.final)
  if (any(lower > upper)) stop("'lower' must be less or equal to 'upper', e.g. lower = 0.01, upper = Inf.")

  # Weighted mean and variance for convenience
  wm <- stats::weighted.mean(z, ct)
  wv <- stats::weighted.mean((z-wm)^2, ct) / sum(ct)  # Variance of the mean; dividing by sum(ct) is not a mistake

  # Checking the spanning condition; 0 must be in the convex hull of z, that is, min(z) < 0 < max(z),
  # or some form of Taylor expansion must be used to avoid this check an simply return a strongly negative
  # value (e.g. Euclidean LL, quartic LL etc.) without w_i > 0 or a Wald statistic instead of LR
  zu <- sort(unique(z))
  # The values cannot be too close to each other because it may break numerical differentiation
  # Keeping only sufficiently distinct ones
  zu <- zu[c(TRUE, abs(diff(zu)) > 16 * .Machine$double.eps)]
  l <- length(zu)
  if (length(zu) >= 2) {
    z12 <- zu[1:2]
    znn <- zu[l-1:0]
  } else {
    z12 <- znn <- rep(zu[1], 2)
  }
  mu.llimit <- if (l > 2) mean(z12) else sum(z12*c(0.9, 0.1))
  mu.rlimit <- if (l > 2) mean(znn) else sum(znn*c(0.1, 0.9))

  spanning <- ((z1 < 0 && zn > 0) && chull.fail == "none") || ((mu.llimit <= 0 && mu.rlimit >= 0) && chull.fail == "taylor")

  if (n.init < 2) { # Codes > 5
    exitcode <- 6
  } else if (z1 == zn) { # The sample is degenerate without variance, no extrapolation possible
    exitcode <- 8
  } else if (!spanning) {  ######### TODO: test with 2 observations
    if ((z1 == 0 || zn == 0) && chull.fail == "none") {
      # mu is on the boundary
      logelr <- -Inf
      lam <- iter <- estim.prec <- f.root <- 0
      converged <- TRUE
      int <- c(0, 0)
      exitcode <- 7
    } else if (chull.fail == "taylor") {
      if (length(zu) < 2) stop("For Taylor extrapolation, at least two unique observations are required.")
      iqr <- stats::IQR(z)
      if (iqr == 0) iqr <- stats::sd(z)
      # With f(x + (0, 1, 2, 3)h), one can estimate the following derivatives
      # 2nd: weights 2 -5 4 -1, error O(h^2)
      # 1st: weights -11/6 3 -1.5 1/3, error O(h^3)
      # With f(x + (-1, 0, 1)h), it is the standard 1 -2 1 for f'' and -0.5 0.5 for f'
      # Step size for noisy functions: for 3 wrong digits, multiply the optimum by 10
      if (mean(sign(zu)) >= 0) {  # All non-negative values, sample average > 0
        # If the spanning condition holds but barely, we need to ignore one point
        # Larger step size for numerically more stable parabola estimates
        mu.limit <- mu.llimit
        # If the points are close, use IQR, but do not go outside the convex hull
        stepsize <- max(diff(z12)*0.01, min(iqr*0.001, diff(z12)*0.4))
      } else if (mean(sign(zu)) <= 0) {
        mu.limit <- mu.rlimit
        stepsize <- max(diff(znn)*0.01, min(iqr*0.001, diff(znn)*0.4))
      } else if (mean(sign(zu)) == 0) {  # There are two unique points?
        stop("Error checking the signs of the data for correct extrapolation. Please report this bug.")
      }
      zgrid <- mu.limit + stepsize*(-1:1)
      w.fp <- c(-0.5, 0, 0.5)  # These numbers are obtained from the pnd package
      w.fpp <- c(1, -2, 1)  # pnd::fdCoef(deriv.order = 1, stencil = -1:1)

      llgrid <- vapply(zgrid, function(m) weightedEL(z = z, mu = m, ct = ct, shift = shift, SEL = SEL, chull.fail = "none")$logelr, numeric(1))
      fp <- sum(llgrid * w.fp) / stepsize
      fpp <- sum(llgrid * w.fpp) / stepsize^2
      # Check with
      # pnd::Grad(function(m) weightedEL(z = z, mu = m, ct = ct)$logelr, zgrid[1],
      #   elementwise = FALSE, vectorised = FALSE, multivalued = FALSE, h = 1e-5)
      f <- llgrid[2]
      abc <- unname(getParabola(x = mu.limit, f, fp, fpp))
      parab  <- function(x) abc[1]*x^2  + abc[2]*x  + abc[3]
      logelr <- parab(0)
      # For z = 1:9
      # xgrid <- seq((z[1]+z[2])/2, (znn[1]+znn[2])/2, length.out = 51)
      # ygrid <- sapply(xgrid, function(m) weightedEL(z = z, mu = m, ct = ct, shift = shift, SEL = SEL)$logelr)
      # plot(xgrid, ygrid, xlim = range(z, 0) + c(-0.25, 0.25), ylim = c(min(ygrid)*1.5, 0), bty = "n")
      # points(zgrid, llgrid, col = 2, pch = 16)
      # xgrid2 <- seq(-0.1, 1.5, length.out = 31)
      # ygrid2 <- parab(xgrid2)
      # lines(xgrid2, ygrid2, col = 2)
      if (z1 == 0 || zn == 0) {
        converged <- TRUE
        exitcode <- 9
      } else {
      exitcode <- 10
      }
    } else if (chull.fail == "wald") {  # Expansion of EL at 0
      if (length(zu) < 2) stop("For Wald extrapolation, at least two unique observations are required.")

      if (mean(sign(zu)) >= 0) {  # All non-negative values, sample average > 0
        # Zero is to the left of the data --> left branch
        mu.limit <- mu.llimit
        gap <- if (l > 2) abs(diff(z12)) * 0.5 else abs(diff(z12)) * 0.05
      } else if (mean(sign(zu)) <= 0) {
        mu.limit <- mu.rlimit
        gap <- if (l > 2) abs(diff(znn)) * 0.5 else abs(diff(znn)) * 0.05
      }
      # TODO: speed up the evaluation; extrapolate only where necessary; check the gap location
      # Extract info from the interpTwo function
      f <- function(mm) vapply(mm, function(m) -2*weightedEL(z = z, mu = m, ct = ct, shift = shift, SEL = SEL, chull.fail = "none")$logelr, numeric(1))
      logelr <- interpTwo(x = 0, f = f, mean = wm, var = wv, at = mu.limit, gap = gap) * -0.5
      # curve(f, 0, 9)
      # abline(v = c(mu.limit, mu.limit - gap), lty = 3)
      # curve((x-wm)^2/wv, 0, 4, col = 2, add = TRUE)
      if (z1 == 0 || zn == 0) {
        converged <- TRUE
        exitcode <- 11
      } else {
        exitcode <- 12
      }
    }
    # Else: do nothing, classical ELR failure
  } else {  # The main 'good' loop: EL may proceed because the spanning condition holds
    # Any requests for self-concordance or other Taylor forms?
    do.taylor <- (!is.na(taylor.order)) && (any(is.finite(upper)) || any(is.finite(lower)))

    negz <- z < 0
    comp <- (ct / N - 1 - shift) / z
    min.lam <- suppressWarnings(max(comp[!negz]))
    max.lam <- suppressWarnings(min(comp[negz]))
    if (!is.finite(min.lam)) min.lam <- -max.lam
    if (!is.finite(max.lam)) max.lam <- -min.lam
    int <- c(min.lam, max.lam)
    int <- int + abs(int) * c(2, -2) * me # To avoid bad rounding
    con <- list(tol = me, maxiter = 100, trace = 0) # Wishing the strictest convergence
    con[names(uniroot.control)] <- uniroot.control

    # Taylor expansion beyond these point to the left and to the right
    # If Taylor expansion is requested by everything is infinite, request left-tail expansion
    if (is.finite(taylor.order) && all(!is.finite(lower))) lower <- ct / N
    if (length(lower) == 1) lower <-  rep(lower, n.final)
    if (length(upper) == 1) upper <- rep(upper, n.final)
    if (length(lower) != n.final || length(upper) != n.final)
      stop("'lower' and 'upper' must have length 1 or length(z).")
    # logELr <- function(lambda) sum(ct * logTaylor(1 + z * lambda + shift, eps = lower, M = upper, der = 0, order = taylor.order, drop = TRUE))
    dllik <- function(lambda) sum(ct * z * logTaylor(1 + z * lambda + shift, eps = lower, M = upper, der = 1, order = taylor.order)[, "deriv1"])
    # xs <- seq(int[1], int[2], length.out = 51)
    # ys <- sapply(xs, logELr)
    # ys1 <- sapply(xs, dllik)
    # plot(xs, ys)
    # plot(xs, ys1)

    lam.list <- tryCatch(brentZero(dllik, interval = int, tol = con$tol, extendInt = if (do.taylor) "yes" else "no",
                                 maxiter = con$maxiter, trace = con$trace),
                       error = function(e) return(NULL))
    # There can be only one kind of warning: maximum iterations reached
    if (!is.null(lam.list)) { # Some result with or without a warning as the second element of the list
      lam <- lam.list$root
      zlam1 <- 1 + z * lam + shift
      wvec <- ct * logTaylor(zlam1, der = 1, order = taylor.order, eps = lower, M = upper)[, "deriv1"]
      if (!SEL) wvec <- wvec / N
      if (return.weights) wts[nonz] <- wvec
      logelr <- -sum(ct * logTaylor(zlam1, der = 0, order = taylor.order, eps = lower, M = upper))
      # Empirical fix for nonsensical probabilities
      # This should not happen unless the spanning condition fails and the Taylor expansion is very inaccurate
      if (any(wvec < 0) && logelr > 0) logelr <- -logelr
      if (any(!is.finite(wvec))) exitcode <- 14

      # brentZero returns the number of iterations times -1 in case it exceeds the maximum number allowed
      converged <- lam.list$iter  >= 0
      estim.prec <- lam.list$estim.prec
      f.root <- lam.list$f.root
      exitcode <- 0

      if (abs(f.root) > sqrt(me)) exitcode <- 1
      # Relative tolerance check for boundary closeness
      int.len <- max.lam - min.lam
      if (min(abs(lam - max.lam), abs(lam - min.lam)) < boundary.tolerance * int.len)
        exitcode <- exitcode + 2
      if (do.taylor && abs(sum(wvec) - 1) > 0.0001) exitcode <- 13
    } else { # The original bad output stays intact, only the exit code updates
      exitcode <- 4
    }

  }


  if (return.weights && any(!is.finite(wts[nonz]))) exitcode <- 9

  msg <- c("FOC not met: |d/dlambda(EL)| > sqrt(Mach.eps)!", # 1
           "Lambda is very close to the boundary! (May happen with extremely small weights.)", # 2
           "Lambda is very close to the boundary and FOC not met: |d/dlambda(EL)| > sqrt(Mach.eps)!", # 3
           "Root finder returned an error!", # 4
           "mu is not strictly in the convex hull of z, and no Taylor expansion is used.", # 5
           "At least two points are needed for meaningful EL computations.", # 6
           "mu lies exactly on the boundary of the convex hull, log ELR(mu) := -Inf.", # 7
           "Observations with substantial counts are identical and equal to mu (degenerate sample).", # 8
           "mu lies exactly on the convex hull boundary, quadratic branch approximation used to obtain ELR(mu) > 0.", # 9
           "mu lies outside the convex hull, quadratic branch approximation used to obtain ELR(mu) > 0.", # 10
           "mu lies exactly on the convex hull boundary, Wald approximation used to obtain ELR(mu) > 0.", # 11
           "mu lies outside the convex hull, Wald approximation used to obtain ELR(mu) > 0.", # 12
           "mu is strictly in the convex hull of z, Taylor expansion yields sum(weights) != 1 and/or weights < 0.",  # 13
           "Lambda search succeeded but some probabilities are not finite (division by zero?)") # 14
  msg <- if (exitcode > 0) msg[exitcode] else "successful convergence within the first-order-condition tolerance"
  if (verbose && exitcode > 0) warning(msg)


  return(list(logelr = logelr, lam = lam, wts = wts,
              converged = converged, iter = iter,
              bracket = int, estim.prec = estim.prec, f.root = f.root,
              exitcode = exitcode, message = msg))
}


#' Euclidean likelihood for vectors of arbitrary dimensions
#'
#' @param z Numeric data vector.
#' @inheritParams weightedEL
#' @param chull.diag Logical: if \code{TRUE}, checks if there is a definite convex hull failure
#'   in at least one dimension (\code{mu} being smaller than the smallest or larger
#'   than the largest element). Note that it does not check if \code{mu} is strictly in the
#'   convex hull because this procedure is much slower and is probably unnecessary.
#'
#' @return A list with the same structure as that in [weightedEL()].
#' @seealso [weightedEL()]
#' @export
#'
#' @examples
#' set.seed(1)
#' z <- cbind(rnorm(10), runif(10))
#' colMeans(z)
#' a <- weightedEuL(z, return.weights = TRUE)
#' a$wts
#' sum(a$wts)  # Unity
#' colSums(a$wts * z)  # Zero
#  smoothemplik:::weightedEuL_cpp(z, return_weights = TRUE)
weightedEuL <- function(z, mu = NULL, ct = NULL, shift = NULL,
                        SEL = TRUE, n.orig = NULL,
                        weight.tolerance = NULL, trunc.to = 0,
                        return.weights = FALSE, verbose = FALSE, chull.diag = FALSE
) {
  if (is.null(dim(z)) || is.data.frame(z)) z <- as.matrix(z, rownames.force = TRUE)
  n <- nrow(z)
  k <- ncol(z)
  if (is.null(mu)) mu <- rep(0, k)
  if (is.null(ct)) ct <- rep(1, n)
  if (is.null(shift)) shift <- rep(0, n)
  if (is.null(n.orig)) n.orig <- n
  if (is.null(weight.tolerance))
    weight.tolerance <- if (!SEL) .Machine$double.eps^(1/3) else max(ct) * sqrt(.Machine$double.eps)
  ret <- weightedEuLCPP(z = z, mu = mu, ct = ct, shift = shift, n_orig = n.orig,
                        weight_tolerance = weight.tolerance, trunc_to = trunc.to, SEL = SEL,
                        return_weights = return.weights, verbose = verbose, chull_diag = chull.diag)
  if (return.weights && !is.null(nz <- rownames(z))) names(ret$wts) <- nz
  return(ret)
}

