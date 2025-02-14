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
  c <- (x[2] * (x[1]^2 * y[3] - x[3]^2 * y[1]) + x[2]^2 * (x[3] * y[1] - x[1] * y[3]) + x[1] * x[3] * y[2] * (x[3] - x[1]))/ xprod
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
#'   to ensure self-concordance of the log-likelihoogetParabola(-1, -6, -1, 4)d.
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
#' @param truncto Counts under \code{weight.tolerance} will be set to this value.
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
#' weightedEL(1:9)
#' weightedEL(1:9, chull.fail = "wald")
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
                       truncto = 0, chull.fail = c("taylor", "wald", "none"),
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
                     if (truncto == 0) "0." else  paste0("+-", truncto, " of appropriate sign.")))
    ct[closeto0 & ct > 0] <- truncto
    ct[closeto0 & ct < 0] <- -truncto
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
  if (n.final < length(z)) {
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
  wv <- stats::weighted.mean((z-wm)^2, ct) / sum(ct)

  # Checking the spanning condition; 0 must be in the convex hull of z, that is, min(z) < 0 < max(z),
  # or some form of Taylor expansion must be used to avoid this check an simply return a strongly negative
  # value (e.g. Euclidean LL, quartic LL etc.) without w_i > 0 or a Wald statistic instead of LR
  spanning <- z1 < 0 && zn > 0
  if (length(z) < 2) { # Codes > 5
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
      l <- length(z)
      zo <- order(z)
      zsort <- z[zo]
      zu <- unique(zsort)
      if (length(zu) < 2) stop("Even with Taylor extrapolation, at least two unique observations are required.")
      z12 <- zu[1:2]
      znn <- zu[length(zu)-1:0]
      iqr <- stats::IQR(z)

      # With f(x + (0, 1, 2, 3)h), one can estimate the following derivatives
      # 2nd: weights 2 -5 4 -1, error O(h^2)
      # 1st: weights -11/6 3 -1.5 1/3, error O(h^3)
      # With f(x + (-1, 0, 1)h), it is the standard 1 -2 1 for f'' and -0.5 0.5 for f'
      # Step size for noisy functions: for 3 wrong digits, multiply the optimum by 10
      if (z1 >= 0) {  # All non-negative values, sample average > 0
        # Larger step size for numerically more stable parabola estimates
        # If the space between 2 points is too small, compute 1-sided derivatives
        mu.limit <- if (l > 2) mean(z12) else sum(z12*c(0.9, 0.1))
        if (diff(z12)*0.01 < iqr*0.001) {  # Points very close to one another
          zgrid <- mu.limit + iqr*0.001*(0:3)
          w.fp <- c(-11/6, 3, -1.5, 1/3)  # These numbers are obtained from the pnd package
          w.fpp <- c(2, -5, 4, -1)  # pnd::fdCoef(deriv.order = 2, stencil = 0:3)
        } else {
          zgrid <- mu.limit + diff(z12)*0.01*c(-1, 0, 1)
          w.fp <- c(-0.5, 0, 0.5)
          w.fpp <- c(1, -2, 1)
        }
      } else if (zn <= 0) {
        mu.limit <- if (l > 2) mean(znn) else sum(znn*c(0.1, 0.9))
        if (diff(znn)*0.01 < iqr*0.001) {
          zgrid <- mu.limit + iqr*0.001*(-3:0)
          w.fp <- c(-1/3, 1.5, -3, 11/6)
          w.fpp <- c(-1, 4, -5, 2)
        } else {
          zgrid <- mu.limit + diff(znn)*0.01*c(-1, 0, 1)
          w.fp <- c(-0.5, 0, 0.5)
          w.fpp <- c(1, -2, 1)
        }
      }

      llgrid <- sapply(zgrid, function(m) weightedEL(z = z, mu = m, ct = ct, shift = shift, SEL = SEL)$logelr)
      ss <- zgrid[2] - zgrid[1]
      fp <- sum(llgrid * w.fp) / ss
      fpp <- sum(llgrid * w.fpp) / ss^2
      # Check with
      # pnd::Grad(function(m) weightedEL(z = z, mu = m, ct = ct)$logelr, zgrid[1],
      #   elementwise = FALSE, vectorised = FALSE, multivalued = FALSE, h = 1e-5)
      if (z1 >= 0) {
        f <- if (length(zgrid) == 4) llgrid[1] else llgrid[2]
      } else if (zn <= 0) {
        f <- if (length(zgrid) == 4) llgrid[4] else llgrid[2]
      }
      abc <- getParabola(x = mu.limit, f, fp, fpp)
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
      logelr <- -0.5 * wm^2 / wv
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
                       warning = function(w) return(list(brentZero(dllik, interval = int, tol = con$tol, extendInt = if (do.taylor) "yes" else "no",
                                                                   maxiter = con$maxiter, trace = con$trace), w)),
                       error = function(e) return(NULL)) # If there is a warning, we still return the object
    if (!is.null(lam.list)) { # Some result with or without a warning as the second element of the list
      if ("warning" %in% class(lam.list[[2]])) {
        if (verbose) print(lam.list[[2]]) # The user should know that went wrong in uniroot()
        lam.list <- lam.list[[1]]
      }
      lam <- lam.list$root
      zlam1 <- 1 + z * lam + shift
      wvec <- ct * logTaylor(zlam1, der = 1, order = taylor.order, eps = lower, M = upper)[, "deriv1"]
      if (!SEL) wvec <- wvec / N
      if (return.weights) wts[nonz] <- wvec
      logelr <- -sum(ct * logTaylor(zlam1, der = 0, order = taylor.order, eps = lower, M = upper))
      # Empirical fix for nonsensical probabilities
      # This should not happen unless the spanning condition fails and the Taylor expansion is very inaccurate
      if (any(wvec < 0) && logelr > 0) logelr <- -logelr
      if (any(!is.infinite(wvec))) exitcode <- 14

      converged <- if ("warning" %in% class(lam.list[[2]])) FALSE else TRUE
      iter <- lam.list$iter
      if (iter == con$maxiter) converged <- FALSE

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


#' Smoothed Empirical Likelihood function value
#'
#' Evaluates SEL function for a given moment function at a certain parameter value.
#'
#' @param rho The moment function depending on parameters and data (and potentially other parameters). Must return a numeric vector.
#' @param theta A parameter at which the moment function is evaluated.
#' @param data A data object on which the moment function is computed.
#' @param sel.weights Either a matrix with valid kernel smoothing weights with rows adding up to 1,
#'   or a list of kernel weights for smoothing where the sum of each element is 1
#'   (must be returned by \code{sparseVectorToList}), or a function that computes
#'   the kernel weights based on the \code{data} argument passed to \code{...}.
#'   If \code{memory.saving} is \code{"partial"} or \code{"full"}, then it must
#'   be a function that computes the kernel weights for the data set.
#' @inheritParams weightedEL
#' @param trim A vector of trimming function values to multiply the output of \code{rho(...)} with. If NULL, no trimming is done.
#' @param minus If TRUE, returns SEL times -1 (for optimisation via minimisation).
#' @param parallel If TRUE, uses \code{parallel::mclapply} to speed up the computation.
#' @param cores The number of cores used by \code{parallel::mclapply}.
#' @param memory.saving A string. \code{"none"} implies no memory-saving tricks,
#' and the entire problem is processed in the computer memory at once (good for
#' sample sizes 2000 and below; if \code{sel.weights} is not provided or is a function,
#' the weight matrix / list is computed at once.). If \code{"full"}, then, the smoothed
#' likelihoods are computed in series, which saves memory but computes kernel weights at
#' every step of a loop, increasing CPU time; the SEL weights, normally found in the rows
#' of the \code{sel.weights} matrix, are computed on the fly. If \code{"partial"}, then,
#' the problem is split into \code{chunks} sub-problems with smaller weight matrices / lists.
#' If \code{parallel} is \code{TRUE}, parallelisation occurs within each chunk.
#' @param chunks The number of chunks into which the weight matrix is split.
#'   Only used if \code{memory.saving} is \code{"partial"}. If there are too many chunks
#'   (resulting in fewer than 2 observations per chunk), then it is treated as if \code{memory.saving} were \code{"full"}.
#' @param print.progress If \code{TRUE}, a progress bar is made to display the evaluation progress in case partial or full memory saving is in place.
#' @param bad.value Replace non-finite individual SEL values with this value.
#'   May be useful if the optimiser does not allow specific non-finite values (like L-BFGS-B).
#' @param attach.attributes If \code{"none"}, returns just the sum of expected likelihoods;
#' otherwise, attaches certain attributes for diagnostics:
#' \code{"ELRs"} for expected likelihoods,
#' \code{"residuals"} for the residuals (moment function values),
#' \code{"lam"} for the Lagrange multipliers lambda in the EL problems,
#' \code{"nabla"} for d/d(lambda)EL (should be close to zero because this must be true for any \code{theta}),
#' \code{"converged"} for the convergence of #' individual EL problems,
#' \code{"exitcode"} for the \code{weightedEL} exit codes (0 for success),
#' \code{"probabilities"} for the matrix of weights (very large, not recommended for sample sizes larger than 2000).
#' @param ... Passed to \code{rho}.
#'
#' @return A scalar with the SEL value and, if requested, attributes containing the diagnostic information attached to it.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(rlnorm(100))
#' # Heteroskedastic DGP
#' y <- abs(1 + 1*x + rnorm(100) * (1 + x + sin(x)))
#' mod.OLS <- lm(y ~ x)
#' rho <- function(theta, ...) y - theta[1] - theta[2]*x  # Moment fn
#' w <- kernelWeights(x, PIT = TRUE, bw = 0.25, kernel = "epanechnikov")
#' w <- w / rowSums(w)
#' image(x, x, w, log = "xy")
#' theta.vals <- list(c(1, 1), coef(mod.OLS))
#' SEL <- function(b, ...) smoothEmplik(rho = rho, theta = b, sel.weights = w)
#' sapply(theta.vals, SEL) # Smoothed empirical likelihood
#' # SEL maximisation
#' b.SEL <- optim(coef(mod.OLS), SEL, method = "BFGS",
#'                control = list(fnscale = -1, reltol = 1e-6))
#' print(b.SEL$par) # Closer to the true value (1, 1) than OLS
#' plot(x, y)
#' abline(1, 1, lty = 2)
#' abline(mod.OLS, col = 2)
#' abline(b.SEL$par, col = 4)
#'
#' # Now we start from (0, 0), for which the Taylor expansion is necessary
#' # because all residuals at this starting value are positive and the
#' # unmodified EL ratio for the test of equality to 0 is -Inf
#' smoothEmplik(rho=rho, theta=c(0, 0), sel.weights=w, chull.fail="none")
#' smoothEmplik(rho=rho, theta=c(0, 0), sel.weights=w)
#' # The next example is very slow; approx. 5 s
#' \dontrun{
#' # Experiment: a small bandwidth so that the spanning condition should fail
#' w <- kernelWeights(x, PIT = TRUE, bw = 0.10, kernel = "epanechnikov")
#' w <- w / rowSums(w)
#' b.SELt <- optim(c(0, 0), SEL, chull.fail = "taylor",  # This is more reliable
#'                 method = "BFGS", control = list(fnscale = -1))
#' b.SELw <- optim(c(0, 0), SEL, chull.fail = "wald",  # This is more reliable
#'                 method = "BFGS", control = list(fnscale = -1))
#' b0grid <- seq(-0.75, 4.5, length.out = 51)
#' b1grid <- seq(-0.75, 2.5, length.out = 51)
#' bgrid <- as.matrix(expand.grid(b0grid, b1grid))
#' selgrid <- unlist(parallel::mclapply(1:nrow(bgrid), function(i)
#'   smoothEmplik(rho, bgrid[i, ], sel.weights = w, chull.fail = "taylor"),
#'     mc.cores = parallel::detectCores()/2-1))
#' selgrid <- matrix(selgrid, nrow = length(b0grid))
#' probs <- c(0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1-10^seq(-4, -16, -2))
#' levs <- qchisq(probs, df = 2)
#' # levs <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
#' labs <- round(levs, 1)
#' cols <- rainbow(length(levs), end = 0.7, v = 0.7)
#' par(mar = c(4, 4, 2, 0) + .1)
#' contour(b0grid, b1grid, -2*(selgrid - max(selgrid, na.rm = TRUE)), levels = levs,
#'         labels = labs, col = cols, lwd = 1.5, bty = "n",
#'         main = "'Safe' likelihood contours", asp = 1)
#' image(b0grid, b1grid, log1p(-2*(selgrid - max(selgrid, na.rm = TRUE))))
#' # The narrow lines are caused by the fact that if two observations are close together
#' # at the edge, the curvature at that point is extreme
#' }
smoothEmplik <- function(rho, theta, data,
                         sel.weights = NULL,
                         chull.fail = c("taylor", "wald", "none"),
                         weight.tolerance = NULL,
                         trim = NULL, minus = FALSE,
                         parallel = FALSE, cores = 1,
                         memory.saving = c("none", "full", "partial"), chunks = 10, print.progress = FALSE,
                         bad.value = -Inf,
                         attach.attributes = c("none", "all", "ELRs", "residuals", "lam", "nabla", "converged", "exitcode", "probabilities"),
                         ...
) {
  chull.fail <- match.arg(chull.fail)
  # Constructing residuals
  rho.series <- rho(theta, data, ...)
  n <- length(rho.series)
  if (any("none" %in% attach.attributes)) attach.attributes <- "none"
  attach.probs <- ("probabilities" %in% attach.attributes) | isTRUE(attach.attributes == "all")

  # Since SEL is a non-parametric method and relies on smoothing with kernel weights, using large matrices for large problems
  # can be terribly inefficient. Instead, we split it into chunks.
  memory.saving <- memory.saving[1]
  taylor.order <- NA
  empliklist <- vector("list", n)
  if (memory.saving == "full") {
      if (!is.function(sel.weights)) stop("smoothEmplik: sel.weights must be a function that accepts an index and a data frame and returns a numeric vector of kernel weights.")
    if (print.progress) pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
      for (i in 1:n) {
        w <- suppressWarnings(as.numeric(sel.weights(i, data)))
        w <- w / sum(w)
        if (is.null(weight.tolerance)) weight.tolerance <- max(w) * sqrt(.Machine$double.eps)
        empliklist[[i]] <- weightedEL(rho.series, ct = w, mu = 0, SEL = TRUE,
                                      weight.tolerance = weight.tolerance, return.weights = attach.probs)
        if (print.progress && (i %% floor(n / 50) == 0)) utils::setTxtProgressBar(pb, i)
      }
    if (print.progress) close(pb)
    } else {
      if (memory.saving == "none") {
        chunks <- 1
        chunk.list <- list(1:n)
      } else if (memory.saving == "partial") {
        if (chunks > n/2) {
          chunks <- max(ceiling(n/1000), 2)
          warning(paste0("Too many chunks for memory saving (fewer than 2 observations per chunk); to reduce overhead, using chunks = ", ceiling(n/1000), "."))
        }
        chunk.labels <- cut(1:n, breaks = chunks, labels = 1:chunks)
        chunk.list <- split(1:n, chunk.labels)
      } else {
        stop("smoothEmplik: Wrong value of the memory.saving argument.")
      }
      if (print.progress) pb <- utils::txtProgressBar(min = 0, max = chunks, style = 3)
      for (k in 1:chunks) {
        inds <- chunk.list[[k]]
        if (k == 1 && memory.saving == "none" && (is.matrix(sel.weights) || is.list(sel.weights))) {
          w <- sel.weights
        } else {
          w <- suppressWarnings(sel.weights(chunk.list[[k]], data))
        }

        # If w is a sparse matrix, it can be converted into a list to save memory (helps in most cases)
        ELiFunc <- if (is.list(w)) {
          function(i) {
            ii <- w[[i]]$idx
            wi <- w[[i]]$ct
            wi <- wi / sum(wi)
            if (is.null(weight.tolerance)) weight.tolerance <- max(wi) * sqrt(.Machine$double.eps)
            weightedEL(rho.series[ii], ct = wi, mu = 0, SEL = TRUE, chull.fail = chull.fail,
                       weight.tolerance = weight.tolerance, return.weights = attach.probs)
          }
        } else if (is.matrix(w)) {
          function(i) {
            wi <- w[i, ]
            wi <- wi / sum(wi)
            if (is.null(weight.tolerance)) weight.tolerance <- max(wi) * sqrt(.Machine$double.eps)
            weightedEL(rho.series, ct = wi, mu = 0, SEL = TRUE, chull.fail = chull.fail,
                       weight.tolerance = weight.tolerance, return.weights = attach.probs)
          }
        } else {
          if (memory.saving == "none")
            stop("smoothEmplik: sel.weights must be a list or a matrix of weights for every observation of the data set, or a function returning a matrix or a list (favouring CPU over memory).") else
              stop("smoothEmplik: sel.weights must be a function returning a matrix or a list (favouring memory over CPU).")
        }
        empliklist[inds] <- if (parallel && cores > 1) parallel::mclapply(X = seq_along(inds), FUN = ELiFunc, mc.cores = cores) else lapply(seq_along(inds), ELiFunc)
        if (print.progress) utils::setTxtProgressBar(pb, k)
      }
      if (print.progress) close(pb)
    }

  if (is.null(trim)) trim <- rep(1, n)
  log.ELR.values <- unlist(lapply(empliklist, "[[", "logelr"))
  if (any(bad <- !is.finite(log.ELR.values))) log.ELR.values[bad] <- bad.value
  if (minus) log.ELR.values <- -log.ELR.values
  log.SELR <- sum(trim * log.ELR.values)
  ret <- log.SELR
  if (isTRUE(attach.attributes == "none")) return(ret)
  if ("ELRs" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "ELRs") <- log.ELR.values
  if ("residuals" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "residuals") <- rho.series
  if ("lam" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "lam") <- unlist(lapply(empliklist, "[[", "lam"))
  if ("nabla" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "nabla") <- unlist(lapply(empliklist, "[[", "f.root"))
  if ("converged" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "converged") <- unlist(lapply(empliklist, "[[", "converged"))
  if ("exitcode" %in% attach.attributes || isTRUE(attach.attributes == "all")) attr(ret, "exitcode") <- unlist(lapply(empliklist, "[[", "exitcode"))
  if (attach.probs) attr(ret, "probabilities") <- lapply(empliklist, "[[", "wts")
  return(ret)
}

#' Constrained smoothed empirical likelihood
#'
#' This function takes two vectors of parameters: free and fixed ones, so that SEL could be optimised with equality constraints.
#' The fixed parameter vector must have full length, and have NA's in places where the free parameters should go.
#' Then, the free parameter values will be inserted in places with NA, and the entire vector passed to smoothEmplik.
#'
#' @param par.free A numeric vector of *only the free parameters* passed to \code{rho}.
#'   In case of constrained optimisation, must be shorter than the number of parameters in the model.
#' @param par.fixed A numeric vector of the same length as \code{par.free}:
#'   numeric values should be in the places where the values are to be kept fixed,
#'   and NA where the free parameters should be. Use NULL or a vector of NA's for unrestricted optimisation.
#' @param ... Passed to \code{\link{smoothEmplik}}.
#'
#' @return A scalar with the constrained (or, if \code{par.fixed = NULL}, unconstrained SEL.)
#' @export
#'
#' @examples
constrSmoothEmplik <- function(par.free = NULL,
                               par.fixed = NULL,
                               ...
) {
  if (is.null(par.free)) { # No free parameters supplied
    theta <- par.fixed
  } else if (isTRUE(all(is.na(par.fixed))) || is.null(par.fixed)) { # All free parameters
    theta <- par.free
  } else { # Some free, some fixed parameters
    theta <- par.fixed
    theta[is.na(theta)] <- par.free
  }
  SEL <- smoothEmplik(theta = theta, ...)
  return(SEL)
}

#' Optimise SEL with or without constraints
#'
#' @param start.values Initial values for unrestricted parameters.
#' @param restricted.params Fixed parameter values used for constrained optimisation, hypothesis testing, power and size calculations etc.
#' @param verbose If TRUE, reports optimisation progress, otherwise remains silent.
#' @param optmethod A string indicating the optimisation method: "nlm" (the default,
#'   invoking \code{stats::nlm}) is slightly quicker, and if it fails, "BFGS" is used.
#' @param nlm.step.max Passed to \code{nlm} if \code{optmethod == "nlm"}; in case
#'   of convergence issues, can be reduced, but if it is too small and 5 \code{nlm} iterations
#'   are done with the maximum step size, optimisation is re-done via BFGS.
#' @param maxit Maximum number of numerical optimiser steps. If it has not converged in this number of steps, fail gracefully with a meaningfull return.
#' @param ... Passed to \code{\link{constrSmoothEmplik}} or, in case all parameters are fixed, \code{\link{smoothEmplik}}.
#'
#' @return A list with the following elements:
#' \describe{
#' \item{par}{The argmax of the SEL that was found by the optimiser.}
#' \item{value}{The SEL value corresponding to \code{par}.}
#' \item{restricted}{A logical vector of the same length as \code{par} indicating if the respective element was kept fixed during the optimisation.}
#' \item{xtimes}{A vector with two components indicating time (in seconds) it took to find the candidate initial value and the final optimum respectively.}
#' }
#' @export
#'
#' @examples
maximiseSEL <- function(start.values = NULL,
                        restricted.params = NULL,
                        verbose = FALSE,
                        optmethod = c("nlm", "BFGS", "Nelder-Mead"),
                        nlm.step.max = NULL,
                        maxit = 50,
                        ...) {
  optError <- function(e = NULL) return(list(code = 5))
  optmethod <- optmethod[1]
  tic0 <- Sys.time()

  # Case 1: If all parameters are restricted, just evaluate the SEL at the given point
  if (!is.null(restricted.params) && all(is.finite(restricted.params))) {
    SEL <- smoothEmplik(theta = restricted.params, ...)
    diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
    return(list(par = restricted.params, value = SEL, restricted = rep(TRUE, length(restricted.params)), code = 0, xtimes = c(initial = 0, opt = diff.opt)))
  }

  # Case 2: # If the model is estimated without restrictions; safeguarging against logical(0)
  if (isTRUE(all(is.na(restricted.params))) || is.null(restricted.params)) {
    restricted <- rep(FALSE, length(start.values))
  } else {
    restricted <- !is.na(restricted.params)
  }

  if (verbose) cat("Maximising SEL.\n")
  opt.controls <- list(trace = as.numeric(verbose), REPORT = 1, maxit = maxit, fnscale = -1)
  if (!is.null(list(...)$minus)) stop("Do not pass 'minus' to 'maximiseSEL'.")
  SELToOptim <- function(theta, ...) constrSmoothEmplik(par.free = theta, par.fixed = restricted.params, minus = FALSE, ...)
  if (optmethod == "nlm") {
    if (is.null(nlm.step.max)) nlm.step.max <- sqrt(sum(start.values^2)) / 3
    optim.SEL <- tryCatch(stats::nlm(p = start.values, f = function(x, ...) -SELToOptim(x, ...), print.level = verbose * 2, stepmax = nlm.step.max, ...), error = optError)
    if (optim.SEL$code %in% c(4, 5)) { # If there are optimisation issues
      warning(paste0("nlm exit code: ", optim.SEL$code, ", restarting with BFGS."))
      rm(optim.SEL) # Purging the faulty result
    }
  }
  if (!exists("optim.SEL")) {
    if (optmethod == "nlm") optmethod <- "BFGS" # Overriding the more fragile method
    optim.SEL <- tryCatch(stats::optim(par = start.values, fn = SELToOptim, control = opt.controls, method = optmethod, ...), error = optError)
  }

  diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
  if (any(restricted)) {
    thetahat <- restricted.params
    thetahat[is.na(thetahat)] <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  } else {
    thetahat <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  }
  SEL <- if (optmethod != "nlm") optim.SEL$value else -optim.SEL$minimum
  if (!isTRUE(abs(SEL) < 1e15)) {
    SEL <- NA
    thetahat <- rep(NA, length(thetahat))
  }
  names(thetahat) <- names(start.values)

  iters <- if (optmethod != "nlm") min(optim.SEL$counts, na.rm = TRUE) else optim.SEL$iterations
  code  <- if (optmethod != "nlm") optim.SEL$convergence else optim.SEL$code
  if (iters == maxit) {
    SEL <- NA
    thetahat <- rep(NA, length(thetahat))
    code <- 4
  }
  if (optmethod %in% c("BFGS", "L-BFGS-B")) {
    if (code == 1) {
      SEL <- NA
      thetahat <- rep(NA, length(thetahat))
      code <- 4
    }
  }
  if (verbose) print(paste0("Unconstrained optimisation finished in ", round(diff.opt, 1), " secs; SEL(", paste(round(thetahat, 2), collapse = ","), ")=", round(SEL, 2), "."))

  results <- list(par = thetahat, value = SEL, restricted = restricted, code = code, iterations = iters, xtimes = diff.opt)
  return(results)
}

#' Convert a weight vector to list
#'
#' This function saves memory (which is crucial in large samples) and allows one to speed up the code by minimising the number of
#' time-consuming subsetting operations and memory-consuming matrix multiplications. We do not want to rely on extra packages for
#' sparse matrix manipulation since the EL smoothing weights are usually fixed at the beginning, and need not be recomputed dynamically,
#' so we recommend applying this function to the rows of a matrix. In order to avoid numerical instability, the weights are trimmed
#' at \code{0.01 / length(x)}. Using too much trimming may cause the spanning condition to fail (the moment function values can have the same sign in some neighbourhoods).
#'
#' @param x A numeric vector or matrix (with many close-to-zero elements).
#' @param trim A trimming function that returns a threshold value below which the weights are ignored. In common applications, this function should tend to 0 as the length of \code{x} increases.
#' @param renormalise Logical: renormalise the sum of weights to one after trimming?
#'
#' @return A list with indices and values of non-zero elements.
#' @export
#'
#' @examples
#' set.seed(1)
#' m <- round(matrix(rnorm(100), 10, 10), 2)
#' m[as.logical(rbinom(100, 1, 0.7))] <- 0
#' sparseVectorToList(m[, 3])
#' sparseMatrixToList(m)
sparseVectorToList <- function(x, trim = NULL, renormalise = FALSE) {
  if (is.null(trim)) trim <- function(x) rep(2*.Machine$double.eps, length(x))
  idx <- which(x >= trim(x))
  x <- x[idx]
  if (renormalise) x <- x / sum(x)
  return(list(idx = idx, ct = x))
}

#' @rdname sparseVectorToList
#' @export
sparseMatrixToList <- function(x, trim = NULL, renormalise = FALSE)
  apply(x, 1, sparseVectorToList, trim = trim, renormalise = renormalise)

#' Construct memory-efficient weights for estimation
#'
#' This function constructs SEL weights with appropriate trimming for numerical stability and optional renormalisation so that the sum of the weights be unity
#'
#' @param x A numeric vector (with many close-to-zero elements).
#' @param bw A numeric scalar or a vector passed to `kernelWeights`.
#' @param trim A trimming function that returns a threshold value below which the weights are ignored. In common applications, this function should tend to 0 as the length of \code{x} increases.
#' @param renormalise Logical; passed to `sparseVectorToList`.
#' @param ... Other arguments pased to \code{kernelWeights}.
#'
#' @return A list with indices of large enough elements.
#' @export
#'
#' @examples
#' getSELWeights(1:5, bw = 2, kernel = "triangular")
getSELWeights <- function(x, bw = NULL, ..., trim = NULL, renormalise = TRUE) {
  sel.weights <- kernelWeights(x, bw = bw, ...)
  sel.weights <- sel.weights / rowSums(sel.weights)
  sel.weights <- apply(sel.weights, 1, sparseVectorToList, trim = trim, renormalise = renormalise)
  return(sel.weights)
}
