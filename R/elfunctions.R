#' Uni-variate empirical likelihood via direct lambda search
#'
#' Empirical likelihood with counts to solve one-dimensional problems efficiently with Brent's root search algorithm.
#' Conducts an empirical likelihood ratio test of the hypothesis that the mean of \code{z} is \code{mu}.
#' The names of the elements in the returned list are consistent with the original R code in \insertCite{owen2017weighted}{smoothemplik}.
#'
#' @param z Numeric data vector.
#' @param mu Hypothesized mean of \code{z} in the moment condition.
#' @param ct Numeric count variable with non-negative values that indicates the multiplicity of observations.
#'   Can be fractional. Very small counts below the threshold \code{weight.tolerance} are zeroed.
#' @param shift The value to add in the denominator (useful in case there are extra Lagrange multipliers): 1 + lambda'Z + shift.
#' @param return.weights Logical: if TRUE, individual EL weights are computed and returned.
#'   Setting this to FALSE gives huge memory savings in large data sets, especially when smoothing is used.
#' @param SEL If \code{FALSE}, then the boundaries for the lambda search are based on the total sum of counts, like in vanilla empirical likelihood,
#' due to formula (2.9) in \insertCite{owen2001empirical}{smoothemplik}, otherwise according to Cosma et al. (2019, p. 170, the topmost formula).
#' @param weight.tolerance Weight tolerance for counts to improve numerical stability
#'   (similar to the ones in Art B. Owen's 2017 code, but adapting to the sample size).
#' @param boundary.tolerance Relative tolerance for determining when the lambda is not an interior
#'   solution because it is too close to the boundary. Corresponds to a fraction of the
#'   interval range length.
#' @param trunc.to Counts under \code{weight.tolerance} will be set to this value.
#'   In most cases, setting this to \code{0} or \code{weight.tolerance} is a viable solution of the zero-denominator problem.
#' @param chull.fail A character: what to do if the convex hull of \code{z} does not contain \code{mu}
#'   (spanning condition does not hold). \code{"taylor"} creates a Taylor approximation
#'   of the log-ELR function near the ends of the sample. \code{"wald"} smoothly transitions
#'   between the log-ELR function into -0.5 * the Wald statistic for the weighted mean of \code{z}.
#'   \code{"adjusted"} invokes the method of \insertCite{chen2008adjusted}{smoothemplik},
#'   where an extra observation is added to ensure that the convex hull contains the mean, and
#'   \code{"balanced"} calls the method of \insertCite{emerson2009calibration}{smoothemplik}
#'   and \insertCite{liu2010adjusted}{smoothemplik} with two extra points.
#' @param uniroot.control A list passed to the \code{brentZero}.
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
#' The actual tolerance of the lambda search in \code{brentZero} is
#' \eqn{2 |\lambda_{\max}| \epsilon_m + \mathrm{tol}/2}{2 * MachEps * l_max + tol/2},
#' where \code{tol} can be set in \code{uniroot.control} and
#' \eqn{\epsilon_m}{MachEps} is \code{.Machine$double.eps}.
#'
#' The sum of log-weights is maximised without Taylor expansion, forcing \code{mu} to be inside
#' the convex hull of \code{z}. If a violation is happening, consider using the \code{chull.fail} argument
#' or switching to Euclidean likelihood via [EuL()].
#'
#' @return A list with the following elements:
#'
#' \describe{
#'   \item{logelr}{Logarithm of the empirical likelihood ratio.}
#'   \item{lam}{The Lagrange multiplier.}
#'   \item{wts}{Observation weights/probabilities (of the same length as \code{z}).}
#'   \item{converged}{\code{TRUE} if the algorithm converged, \code{FALSE} otherwise (usually means that \code{mu} is not within the range of \code{z}, i.e. the one-dimensional convex hull of \code{z}).}
#'   \item{iter}{The number of iterations used (from \code{brentZero}).}
#'   \item{bracket}{The admissible interval for lambda (that is, yielding weights between 0 and 1).}
#'   \item{estim.prec}{The approximate estimated precision of lambda (from \code{brentZero}).}
#'   \item{f.root}{The value of the derivative of the objective function w.r.t. lambda at the root (from \code{brentZero}). Values \code{> sqrt(.Machine$double.eps)} indicate convergence problems.}
#'   \item{exitcode}{An integer indicating the reason of termination.}
#'   \item{message}{Character string describing the optimisation termination status.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [EL()]
#'
#' @examples
#' # Figure 2.4 from Owen (2001) -- with a slightly different data point
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' set.seed(1)
#' system.time(r1 <- replicate(40, EL(sample(earth, replace = TRUE), mu = 5.517)))
#' set.seed(1)
#' system.time(r2 <- replicate(40, EL0(sample(earth, replace = TRUE), mu = 5.517)))
#' plot(apply(r1, 2, "[[", "logelr"), apply(r1, 2, "[[", "logelr") - apply(r2, 2, "[[", "logelr"),
#'      bty = "n", xlab = "log(ELR) computed via dampened Newthon method",
#'      main = "Discrepancy between EL and EL0", ylab = "")
#' abline(h = 0, lty = 2)
#'
#' # Handling the convex hull violation differently
#' EL0(1:9, chull.fail = "none")
#' EL0(1:9, chull.fail = "taylor")
#' EL0(1:9, chull.fail = "wald")
#'
#' # Interpolation to well-defined branches outside the convex hull
#' mu.seq <- seq(-1, 7, 0.1)
#' wEL1 <- -2*sapply(mu.seq, function(m) EL0(1:9, mu = m, chull.fail = "none")$logelr)
#' wEL2 <- -2*sapply(mu.seq, function(m) EL0(1:9, mu = m, chull.fail = "taylor")$logelr)
#' wEL3 <- -2*sapply(mu.seq, function(m) EL0(1:9, mu = m, chull.fail = "wald")$logelr)
#' plot(mu.seq, wEL1)
#' lines(mu.seq, wEL2, col = 2)
#' lines(mu.seq, wEL3, col = 4)
#'
#' # Warning: depending on the compiler, the discrepancy between EL and EL0
#' # can be one million (1) times larger than the machine epsilon despite both of them
#' # being written in pure R
#' # The results from Apple clang-1400.0.29.202 and Fortran GCC 12.2.0 are different from
#' # those obtained under Ubuntu 22.04.4 + GCC 11.4.0-1ubuntu1~22.04,
#' # Arch Linux 6.6.21 + GCC 14.1.1, and Windows Server 2022 + GCC 13.2.0
#' out1 <- EL(earth, mu = 5.517)[1:4]
#' out2 <- EL0(earth, mu = 5.517, return.weights = TRUE)[1:4]
#' print(c(out1$lam, out2$lam), 16)
#'
#' # Value of lambda                                 EL                  EL0
#' # aarch64-apple-darwin20         -1.5631313955??????   -1.5631313957?????
#' # Windows, Ubuntu, Arch           -1.563131395492627   -1.563131395492627
#' @export
EL0 <- function(z, mu = NULL, ct = NULL, shift = NULL, return.weights = FALSE, SEL = FALSE,
                        weight.tolerance = NULL, boundary.tolerance = 1e-9, trunc.to = 0,
                        chull.fail = c("taylor", "wald", "adjusted", "balanced", "none"),
                        uniroot.control = list(),
                        verbose = FALSE
) {
  if (NCOL(z) > 1) stop("Only one-dimensional vectors or matrices are supported.")
  if (is.data.frame(z)) z <- as.matrix(z)
  if (!is.null(dim(z))) z <- drop(z)
  if (any(!is.finite(z))) stop("Non-finite observations (NA, NaN, Inf) are not welcome.")


  n <- length(z)

  if (is.null(mu)) mu <- 0
  if (is.null(ct)) ct <- rep(1, n)
  if (any(!is.finite(ct))) stop("Non-finite weights (NA, NaN, Inf) are not welcome.")
  if (min(ct) < 0) stop("Negative weights are present.")
  if (sum(ct) <= 0) stop("The total sum of EL weights must be positive.")
  n.orig <- n
  if (is.null(shift)) shift <- rep(0, n)
  if (is.null(weight.tolerance))
    weight.tolerance <- if (!SEL) .Machine$double.eps^(1/3) else max(ct) * sqrt(.Machine$double.eps)
  chull.fail <- match.arg(chull.fail)

  # If originally the weights were too small, too many points would be truncated
  # Warn if any non-zero weights are smaller than weight.tolerance
  closeto0 <- (abs(ct) < weight.tolerance)
  if (any(closeto0)) {
    if (verbose) warning(paste0("Counts closer to 0 than ", sprintf("%1.2e", weight.tolerance),
                                " have been replaced with ", if (trunc.to == 0) "0." else trunc.to))
    ct[closeto0 & ct > 0] <- trunc.to
  }

  if (return.weights) {
    wts <- numeric(n.orig)
    names(wts) <- names(z)
  } else {
    wts <- NULL
  }

  # Not all observations contribute meaningfully to the likelihood; tiny weights push lambda to the boundary
  nonz <- which(ct != 0)
  n.final <- length(nonz)
  if (n.final < n) {
    z <- z[nonz]
    ct <- ct[nonz]
    shift <- shift[nonz]
  }


  if (SEL) ct <- ct / sum(ct) # We might have truncated some weights, so re-normalisation is needed!
  # The denominator for EL with counts is the sum of total counts, and for SEL, it is the number of observations!
  n.denom <- if (SEL) n else sum(ct)
  if (n.denom <= 0) stop(paste0("Total weights after tolerance checks (", n.denom, ") must be positive (check the counts and maybe decrease 'weight.tolerance', which is now ", sprintf("%1.1e", weight.tolerance), ".\n"))

  # Enforcing the moment condition
  z <- z - mu

  # Finally, implement Adjusted and Balanced EL to maintain the sample average
  # Handling AEL and BAEL at the very end is correct because the formulae are simpler
  if (chull.fail %in% c("adjusted", "balanced")) {
    an <- computeBartlett(z) * 0.5  # Unweighted -- because there is no theory on weighted AEL
    if (an < 0) stop("EL0: Bartlett factor a_n < 0 -- please report this bug on GitHub.")
    zbar <- trimmed.weighted.mean(z, trim = 0.1, w = ct)
    if (chull.fail == "adjusted") {
      z <- c(z, -zbar * an)
      ct <- if (!SEL) c(ct, 1) else c(ct, 1)
      shift <- c(shift, 0)
    } else {
      z <- c(z, -zbar * an, 2 * mean(z) + zbar * an)
      ct <- if (!SEL) c(ct, 1, 1) else c(ct, 1, 1)
      if (SEL) ct <- ct / sum(ct)
      shift <- c(shift, 0, 0)
    }
  }

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
  exitcode <- 5L

  # Checking the spanning condition; 0 must be in the convex hull of z, that is, min(z) < 0 < max(z),
  # or some form extrapolation must be used to avoid this check an simply return a strongly negative
  # value (e.g. Euclidean L, Balanced EL etc.) or a Wald statistic instead of LR
  # NB: here, we declare SC failure even in the stronger sense: 'one must extrapolate'
  # due the zero being too close to the boundary
  zu <- sort(unique(z))
  # The values cannot be too close to each other because it may break numerical differentiation
  # Keeping only sufficiently distinct ones
  # TODO: add relative difference, too
  abs.diff <- c(1, abs(diff(zu)))
  zu <- zu[abs.diff > 64*.Machine$double.eps]
  l <- length(zu)
  if (length(zu) >= 2) {
    z12 <- zu[1:2]
    znn <- zu[(l-1):l]
  } else {
    z12 <- znn <- rep(zu[1], 2)
  }

  ExEL <- chull.fail %in% c("taylor", "wald")
  if (ExEL) {  # Extrapolation will be done -- the limits are required
    mu.llimit <- if (l > 2) mean(z12) else sum(z12*c(0.9, 0.1))
    mu.rlimit <- if (l > 2) mean(znn) else sum(znn*c(0.1, 0.9))
    if (chull.fail == "wald") {
      wm <- stats::weighted.mean(z, ct)  # Over-writing what had been calculated before
      wv <- stats::weighted.mean((z-wm)^2, ct) / sum(ct)
    }
    left.extrap  <- (mu.llimit > 0)
    right.extrap <- (mu.rlimit < 0)
    if (left.extrap && right.extrap) stop("Extrapolation must be one-sided.")

    # Optimisation 1: extrapolate only the right end
    if (left.extrap) {
      # All observations are to the right -- zero would be in the left branch
      z <- -z
      zu <- sort(unique(z))
      z1 <- min(z)
      zn <- max(z)
      z12 <- zu[1:2]
      znn <- zu[(l-1):l]
      mu.llimit <- if (l > 2) mean(z12) else sum(z12*c(0.9, 0.1))
      mu.rlimit <- if (l > 2) mean(znn) else sum(znn*c(0.1, 0.9))
      if (chull.fail == "wald") wm <- -wm
    }
  }

  spanning <- ((!ExEL) & (z1 <= 0 && zn >= 0)) || (ExEL && (mu.llimit <= 0 && mu.rlimit >= 0))

  if (n < 2) { # Codes > 5
    exitcode <- 6L
  } else if (z1 == zn) { # The sample is degenerate without variance, no extrapolation possible
    exitcode <- 8L
  } else if (chull.fail == "none" && (z1 == 0 || zn == 0)) {
    # mu is on the boundary -- special case (still in the book)
    logelr <- -Inf
    lam <- iter <- estim.prec <- f.root <- 0
    converged <- TRUE
    int <- c(0, 0)
    exitcode <- 7L
  } else if (!spanning) {  ######### TODO: test with 2 observations

    if (chull.fail == "taylor") {
      if (length(zu) < 2) stop("For Taylor extrapolation, at least two unique observations are required.")
      ma <- stats::median(abs(z))
      if (ma == 0) ma <- stats::IQR(z)
      if (ma == 0) ma <- stats::sd(z)
      if (ma == 0) stop("No variability in remaining 'z', calculations impossible.")
      if (mean(sign(zu)) == 0) {  # There are two unique points?
        stop("Error checking the signs of the data for correct extrapolation. Please report this bug.")
      }
      # Anchor the parabola at the end that is closest to 0 after the (possible) sign flip
      mu.limit <- if (abs(mu.llimit) < abs(mu.rlimit)) mu.llimit else mu.rlimit
      # Apply the parabola formula for 3 points using f, f', f'' from a 4-point stencil with a larger step size
      stepsize <- max(diff(znn)*0.01, ma*.Machine$double.eps^0.25)
      # Cap so that mu.limit - 3*stepsize >= z1 (stay inside data)
      max.step <- (mu.limit - z1)/3
      if (max.step <= 0) stop("EL0: cannot form left stencil (mu.limit <= z_min).")
      stepsize <- min(stepsize, 0.5*max.step)
      zgrid <- mu.limit + stepsize*(-3:0)
      w.fp <- c(-2, 9, -18, 11) / 6  # These numbers are obtained from the pnd package
      w.fpp <- c(-1, 4, -5, 2)    # pnd::fdCoef(deriv.order = 2, stencil = -3:0)

      llgrid <- vapply(zgrid, function(m) EL0(z = z, mu = m, ct = ct, shift = shift, SEL = SEL, chull.fail = "none")$logelr, numeric(1))
      fp <- sum(llgrid * w.fp) / stepsize
      fpp <- sum(llgrid * w.fpp) / stepsize^2
      # Check with
      # pnd::Grad(function(m) EL0(z = z, mu = m, ct = ct)$logelr, zgrid[1],
      #   elementwise = FALSE, vectorised = FALSE, multivalued = FALSE, h = 1e-5)
      abc <- getParabola(x = mu.limit, f = llgrid[4], fp = fp, fpp = fpp)
      parab  <- function(x) abc[1]*x^2  + abc[2]*x  + abc[3]
      logelr <- parab(0)  # Just c, the intercept
      # For z = 1:9
      # xgrid <- seq((z[1]+z[2])/2, (znn[1]+znn[2])/2, length.out = 51)
      # ygrid <- sapply(xgrid, function(m) EL0(z = z, mu = m, ct = ct, shift = shift, SEL = SEL)$logelr)
      # plot(xgrid, ygrid, xlim = range(z, 0) + c(-0.25, 0.25), ylim = c(min(ygrid)*1.5, 0), bty = "n")
      # points(zgrid, llgrid, col = 2, pch = 16)
      # xgrid2 <- seq(-0.1, 1.5, length.out = 31)
      # ygrid2 <- parab(xgrid2)
      # lines(xgrid2, ygrid2, col = 2)
      if (z1 == 0 || zn == 0) {
        converged <- TRUE
        exitcode <- 9L
      } else {
      exitcode <- 10L
      }
    } else if (chull.fail == "wald") {  # Wald approximation to ELR at 0 = squared t-stat
      if (length(zu) < 2) stop("For Wald extrapolation, at least two unique observations are required.")

      mu.limit <- if (abs(mu.llimit) < abs(mu.rlimit)) mu.llimit else mu.rlimit
      gap <- if (l > 2) abs(diff(znn)) * 0.5 else abs(diff(znn)) * 0.05

      # TODO: speed up the evaluation; extrapolate only where necessary; check the gap location
      # Extract info from the interpTwo function
      f <- function(mm) vapply(mm, function(m) -2*EL0(z = z, mu = m, ct = ct, shift = shift, SEL = SEL, chull.fail = "none")$logelr, numeric(1))
      logelr <- -0.5 * interpTwo(x = 0, f = f, mean = wm, var = wv, at = mu.limit, gap = gap)
      # curve(f, 0, 9)
      # abline(v = c(mu.limit, mu.limit - gap), lty = 3)
      # curve((x-wm)^2/wv, 0, 4, col = 2, add = TRUE)
      if (z1 == 0 || zn == 0) {
        converged <- TRUE
        exitcode <- 9L
      } else {
        exitcode <- 10L
      }
    }
    # Else: do nothing, classical ELR failure
  } else {  # The main 'good' loop: EL may proceed because the spanning condition holds
    negz <- z < 0
    comp <- (ct / n.denom - 1 - shift) / z
    min.lam <- suppressWarnings(max(comp[!negz]))
    max.lam <- suppressWarnings(min(comp[negz]))
    if (!is.finite(min.lam)) min.lam <- -max.lam
    if (!is.finite(max.lam)) max.lam <- -min.lam
    int <- c(min.lam, max.lam)
    int <- int + abs(int) * c(2, -2) * .Machine$double.eps # To avoid bad rounding
    con <- list(tol = .Machine$double.eps, maxiter = 100, trace = 0) # Wishing the strictest convergence
    con[names(uniroot.control)] <- uniroot.control

    # dllik <- function(lambda) sum(ct * z * logTaylor(1 + z * lambda + shift, lower = lower, upper = upper, der = 1, order = taylor.order))
    dllik <- function(lambda) sum(ct * z * dlog(1 + z * lambda + shift, d = 1))
    # xs <- seq(int[1], int[2], length.out = 51)
    # ys <- sapply(xs, logELr)
    # ys1 <- sapply(xs, dllik)
    # plot(xs, ys)
    # plot(xs, ys1)

    lam.list <- tryCatch(brentZero(dllik, interval = int, tol = con$tol, extendInt = "no",
                                   maxiter = con$maxiter, trace = con$trace),
                         error = function(e) return(NULL))
    # There can be only one kind of warning: maximum iterations reached
    if (!is.null(lam.list)) { # Some result with or without a warning as the second element of the list
      lam <- lam.list$root
      zlam1 <- 1 + z * lam + shift
      wvec <- ct * dlog(zlam1, d = 1)
      if (!SEL) wvec <- wvec / n.denom
      if (return.weights) wts[nonz] <- wvec
      logelr <- -sum(ct * dlog(zlam1, d = 0))
      # Empirical fix for nonsensical probabilities
      # This should not happen unless the spanning condition fails and the Taylor expansion is very inaccurate
      if (any(wvec < 0) && logelr > 0) logelr <- -logelr
      if (any(!is.finite(wvec))) exitcode <- 12L

      # brentZero returns the number of iterations times -1 in case it exceeds the maximum number allowed
      converged <- lam.list$iter  >= 0
      estim.prec <- lam.list$estim.prec
      f.root <- lam.list$f.root
      iter <- lam.list$iter
      exitcode <- 0L

      if (abs(f.root) > sqrt(.Machine$double.eps)) exitcode <- 1L
      # Relative tolerance check for boundary closeness
      int.len <- max.lam - min.lam
      if (min(abs(lam - max.lam), abs(lam - min.lam)) < boundary.tolerance * int.len)
        exitcode <- exitcode + 2L
      if (abs(sum(wvec) - 1) > 1e-6) exitcode <- 11L
    } else { # The original bad output stays intact, only the exit code updates
      exitcode <- 4L
    }

  }


  if (return.weights && any(!is.finite(wts[nonz]))) exitcode <- 9L

  msg <- c("FOC not met: |d/dlambda(EL)| > sqrt(Mach.eps)!", # 1
           "Lambda is very close to the boundary! (May happen with extremely small weights.)", # 2
           "Lambda is very close to the boundary and FOC not met: |d/dlambda(EL)| > sqrt(Mach.eps)!", # 3
           "Root finder returned an error!", # 4
           "mu is not strictly in the convex hull of z, and no Taylor expansion is used.", # 5
           "At least two points are needed for meaningful EL computations.", # 6
           "mu lies exactly on the boundary of the convex hull, log ELR(mu) := -Inf.", # 7
           "Observations with substantial counts are identical and equal to mu (degenerate sample).", # 8
           "mu lies exactly on the convex hull boundary, approximation used to obtain ELR(mu) > 0.", # 9
           "mu lies outside the convex hull, approximation used to obtain ELR(mu) > 0.", # 10
           "mu is strictly in the convex hull of z, sum(weights) != 1 and/or weights < 0.",  # 11
           "Lambda search succeeded but some probabilities are not finite (division by zero?)") # 12
  msg <- if (exitcode > 0L) msg[exitcode] else "successful convergence within the first-order-condition tolerance"
  if (verbose && exitcode > 0L) warning(msg)


  return(list(logelr = logelr, lam = lam, wts = wts,
              converged = converged, iter = iter,
              bracket = int, estim.prec = estim.prec, f.root = f.root,
              exitcode = exitcode, message = msg))
}


#' Self-concordant multi-variate empirical likelihood with counts
#'
#' Implements the empirical-likelihood-ratio test for the mean of the coordinates of \code{z}
#' (with the hypothesised value \code{mu}). The counts need not be integer;
#' in the context of local likelihoods, they can be kernel observation weights.
#'
#' @source This original code was written for \insertCite{owen2013self}{smoothemplik}
#' and [published online](https://artowen.su.domains/empirical/) by Art B. Owen
#' (March 2015, February 2017). The present version was rewritten in \code{Rcpp} and
#' slightly reworked to contain fewer inner functions and loops.
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
#' \code{alpha = 0.3} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).
#'
#' The argument names, except for \code{lambda.init}, are matching the original names in Art B. Owen's implementation.
#' The highly optimised one-dimensional counterpart, \code{EL0}, is designed to return a faster
#' and a more accurate solution in the one-dimensional case.
#'
#' @param z A numeric vector or a matrix with one data vector per column.
#' @param ct A numeric vector of non-negative counts.
#' @param mu Hypothesised mean, default (0 ... 0) in R^ncol(z)
#' @param lambda.init Starting lambda, default (0 ... 0)
#' @param SEL If \code{FALSE}, the default weight tolerance is \code{MachEps^(1/3)}, otherwise
#'   it is \code{MachEps^(1/2)} of the maximum count.
#' @param return.weights Logical: if \code{TRUE}, returns the empirical probabilities. Default is memory-saving (\code{FALSE}).
#' @param lower Lower cut-off for [logTaylor()], default \code{1/nrow(z)}
#' @param upper Upper cutoff for [logTaylor()], default \code{Inf}
#' @param order Positive integer such that the Taylor approximation of this order to log(x) is self-concordant; usually 4 or higher. Passed to [logTaylor()].
#' @param weight.tolerance Weight tolerance for counts to improve numerical stability
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
#' @seealso [logTaylor()], [EL0()]
#'
#' @examples
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' EL(earth, mu = 5.517, verbose = TRUE) # 5.517 is the modern accepted value
#'
#' # Linear regression through empirical likelihood
#' coef.lm <- coef(lm(mpg ~ hp + am, data = mtcars))
#' xmat <- cbind(1, as.matrix(mtcars[, c("hp", "am")]))
#' yvec <- mtcars$mpg
#' foc.lm <- function(par, x, y) {  # The sample average of this
#'   resid <- y - drop(x %*% par)   # must be 0
#'   resid * x
#' }
#' minusEL <- function(par) -EL(foc.lm(par, xmat, yvec), itermax = 10)$logelr
#' coef.el <- optim(c(mean(yvec), 0, 0), minusEL)$par
#' abs(coef.el - coef.lm) / coef.lm  # Relative difference
#'
#' # Likelihood ratio testing without any variance estimation
#' # Define the profile empirical likelihood for the coefficient on am
#' minusPEL <- function(par.free, par.am)
#'   -EL(foc.lm(c(par.free, par.am), xmat, yvec), itermax = 20)$logelr
#' # Constrained maximisation assuming that the coef on par.am is 3.14
#' coef.el.constr <- optim(coef.el[1:2], minusPEL, par.am = 3.14)$par
#' print(-2 * EL(foc.lm(c(coef.el.constr, 3.14), xmat, yvec))$logelr)
#' # Exceeds the critical value qchisq(0.95, df = 1)
#' @export
EL <- function(z, mu = NULL, ct = NULL, lambda.init = NULL, SEL = FALSE,
                       return.weights = FALSE, lower = NULL, upper = NULL,
                       order = 4L, weight.tolerance = NULL,
                       thresh = 1e-16, itermax = 100L, verbose = FALSE,
                       alpha = 0.3, beta = 0.8, backeps = 0) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  n <- nrow(z)
  d <- ncol(z)
  if (length(mu) == 0) mu <- rep(0, d)
  if (length(mu) != d) stop("The length of mu must be the same as the dimension of z.")
  if (length(lambda.init) == 0) lambda.init <- rep(0, d)
  if (length(lambda.init) != d) stop("The length of mu must be the same as the dimension of z.")

  if (is.null(lower)) lower <- rep(1/n, n)
  if (is.null(upper)) upper <- rep(Inf, n)

  if (is.null(weight.tolerance))
    weight.tolerance <- if (!SEL) .Machine$double.eps^(1/3) else max(ct) * sqrt(.Machine$double.eps)

  if (is.null(ct)) ct <- rep(1, n)
  if (min(ct) < 0) stop("Negative weights are not allowed.")
  if (SEL) ct <- ct / sum(ct)
  if (any(0 < ct & ct < weight.tolerance)) {
    if (verbose) warning(paste("Positive counts below", weight.tolerance, "have been replaced by zero."))
    ct[ct < weight.tolerance] <- 0
    if (SEL) ct <- ct / sum(ct)  # Re-normalising again
  }
  if (sum(ct) <= 0) stop("Total weight must be positive.")

  ELCPP(z = z, ct = ct, mu = mu, lambda_init = lambda.init,
        return_weights = return.weights, lower = lower, upper = upper,
        order = order, weight_tolerance = weight.tolerance, thresh = thresh,
        itermax = itermax, verbose = verbose,
        alpha = alpha, beta = beta, backeps = backeps)
}

#' Compute empirical likelihood on a trajectory
#'
#' @param z Passed to \code{EL}.
#' @param ct Passed to \code{EL}.
#' @param mu0 Starting point of trajectory
#' @param mu1 End point of trajectory
#' @param N Number of segments into which the path is split (i. e. \code{N+1} steps are used).
#' @param verbose Logical: report iteration data?
#' @param verbose.solver Logical: report internal iteration data from the optimiser? Very verbose.
#' @param ... Passed to \code{EL}.
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
#' EL(earth, mu = 5.1,  verbose = TRUE)
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

  lam0 <- NULL
  elr <- ec <- it <- nd <- gn <- numeric(N+1)
  m <- l <- matrix(NA, N+1, ncol = d)
  for (i in 0:N) {
    mui <- (i*mu1 + (N-i)*mu0) / N
    if (i > 0) lam0 <- if (all(is.finite(x$lam))) drop(x$lam) else NULL
    x <- EL(z = z, ct = ct, mu = mui, lambda.init = lam0, verbose = verbose.solver, ...)
    m[i+1, ] <- mui
    l[i+1, ] <- x$lam
    elr[i+1] <- x$logelr
    ec[i+1] <- x$exitcode
    it[i+1] <- x$iter
    nd[i+1] <- x$ndec
    gn[i+1] <- x$gradnorm
    if (verbose) cat("Point ", i, "/", N, ", ", if (x$exitcode != 0) "NOT " else "", "converged",
                     ", log(ELR) = ", x$logelr, "\n", sep = "")
  }
  ans <- data.frame(mu = m, logelr = elr, lambda = l, exitcode = ec, iter = it, ndec = nd, gradnorm = gn)
  return(ans)
}


#' Multi-variate Euclidean likelihood with analytical solution
#'
#' @param z Numeric data vector.
#' @param vt Numeric vector: non-negative variance weights for estimating the conditional
#'   variance of \code{z}. Probabilities are returned only for the observations where \code{vt > 0}.
#' @inheritParams EL0
#' @param chull.diag Logical: if \code{TRUE}, checks if there is a definite convex hull failure
#'   in at least one dimension (\code{mu} being smaller than the smallest or larger
#'   than the largest element). Note that it does not check if \code{mu} is strictly in the
#'   convex hull because this procedure is much slower and is probably unnecessary.
#'
#' The arguments \code{ct} and \code{vt} are responsible for smoothing of the moment function
#' and conditional variance, respectively. The objective function is
#' \deqn{\min_{p_{ij}} \frac1n \sum_{i=1}^n \sum_{j=1}^n \mathbb{I}_{ij} \frac{(p_{ij} -
#'   c_{ij})^2}{2v_{ij}}}{min_(p_ij) 1/n * sum_i sum_j I_ij (p_ij - c_ij)^2 / (2v_{ij})},
#' where \eqn{\mathbb{I}_{ij}}{I_ij} is 1 if \eqn{v_{ij} \ne 0}{v_ij != 0}.
#'
#' This estimator is numerically equivalent to the Sieve Minimum Distance estimator
#' of \insertCite{ai2003efficient}{smoothemplik} with kernel sieves, but this interface
#' provides more flexibility through the two sets of weights. If \code{ct} and
#' \code{vt} are not provided, their default value is set to 1, and the resulting
#' estimator is the CUE-GMM estimator: a quadratic form in which the unconditional
#' mean vector is weighted by the inverse of the unconditional variance.
#'
#' @return A list with the same structure as that in [EL()].
#' @seealso [EL()]
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' set.seed(1)
#' z <- cbind(rnorm(10), runif(10))
#' colMeans(z)
#' a <- EuL(z, return.weights = TRUE)
#' a$wts
#' sum(a$wts)  # Unity
#' colSums(a$wts * z)  # Zero
EuL <- function(z, mu = NULL, ct = NULL, vt = NULL, shift = NULL,
                        SEL = TRUE,
                        weight.tolerance = NULL, trunc.to = 0,
                        return.weights = FALSE, verbose = FALSE, chull.diag = FALSE
) {
  if (is.null(dim(z)) || is.data.frame(z)) z <- as.matrix(z, rownames.force = TRUE)
  n <- nrow(z)
  k <- ncol(z)
  if (is.null(mu)) mu <- rep(0, k)
  if (is.null(ct)) ct <- rep(1, n)
  if (is.null(vt)) vt <- rep(1, n)
  if (is.null(shift)) shift <- rep(0, n)
  n.orig <- n
  if (is.null(weight.tolerance))
    weight.tolerance <- if (!SEL) .Machine$double.eps^(1/3) else max(ct) * sqrt(.Machine$double.eps)
  ret <- EuLCPP(z = z, mu = mu, ct = ct, vt = vt, shift = shift, n_orig = n.orig,
                weight_tolerance = weight.tolerance, trunc_to = trunc.to, SEL = SEL,
                return_weights = return.weights, verbose = verbose, chull_diag = chull.diag)
  if (return.weights && !is.null(nz <- rownames(z))) names(ret$wts) <- nz
  return(ret)
}

