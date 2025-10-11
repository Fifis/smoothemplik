#' Uni-variate empirical likelihood via direct lambda search
#'
#' Empirical likelihood with counts to solve one-dimensional problems efficiently with Brent's root search algorithm.
#' Conducts an empirical likelihood ratio test of the hypothesis that the mean of \code{z} is \code{mu}.
#' The names of the elements in the returned list are consistent with the original R code
#' in \insertCite{owen2017weighted}{smoothemplik}.
#'
#' @param z A numeric vector containing the observations.
#' @param mu Hypothesised mean of \code{z} in the moment condition.
#' @param ct Numeric count variable with non-negative values that indicates the multiplicity of observations.
#'   Can be fractional. Very small counts below the threshold \code{weight.tolerance} are zeroed.
#' @param shift The value to add in the denominator (useful in case there are extra Lagrange multipliers):
#'   \eqn{1 + \lambda'Z + shift}{1 + lambda'Z + shift.}.
#' @param renormalise If \code{FALSE}, then uses the total sum of counts as the number of observations,
#'   like in vanilla empirical likelihood, due to formula (2.9) in \insertCite{owen2001empirical}{smoothemplik},
#'   otherwise re-normalises the counts to 1 according to \insertCite{cosma2019inference}{smoothemplik}
#'   (see p. 170, the topmost formula).
#' @param return.weights Logical: if \code{TRUE}, returns the empirical probabilities.
#'   Default is memory-saving (\code{FALSE}).
#' @param weight.tolerance Weight tolerance for counts to improve numerical stability
#'   (defaults to \code{sqrt(.Machine$double.eps)} times the maximum weight).
#' @param boundary.tolerance Relative tolerance for determining when lambda is not an interior
#'   solution because it is too close to the boundary. Unit: fraction of the feasble bracket length.
#' @param trunc.to Counts under \code{weight.tolerance} will be set to this value.
#'   In most cases, setting this to \code{0} (default) or \code{weight.tolerance}
#'   is a viable solution for the zero-denominator problem.
#' @param deriv Logical: if \code{TRUE}, computes and returns the first two derivatives of log-ELR w.r.t. \code{mu}.
#' @param log.control List of arguments passed to [logTaylor()].
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
#' where \code{tol} is \code{.Machine$double.eps} and
#' \eqn{\epsilon_m}{MachEps} is \code{.Machine$double.eps}.
#'
#' The sum of log-weights is maximised without Taylor expansion, forcing \code{mu} to be inside
#' the convex hull of \code{z}. If a violation is happening, consider using
#' \code{log.control(order = 4)} or switching to Euclidean likelihood via [EuL()].
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
#'   \item{deriv}{If requested, the first two derivatives of log-ELR w.r.t. \code{mu}}
#'   \item{exitcode}{An integer indicating the reason of termination.}
#'   \item{message}{Character string describing the optimisation termination status.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso [EL1()] for multi-variate EL based on minimisation w.r.t. lambda.
#'
#' @examples
#' # Figure 2.4 from Owen (2001) -- with a slightly different data point
#' earth <- c(
#'   5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
#'   5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
#' )
#' # Root searching (EL0) is faster than minimisation w.r.t. lambda (EL1)
#' set.seed(1)
#' system.time(r0 <- replicate(40, EL0(sample(earth, replace = TRUE), mu = 5.517)))
#' set.seed(1)
#' system.time(r1 <- replicate(40, EL1(sample(earth, replace = TRUE), mu = 5.517)))
#' plot(apply(r0, 2, "[[", "logelr"), apply(r1, 2, "[[", "logelr") - apply(r0, 2, "[[", "logelr"),
#'      bty = "n", xlab = "log(ELR) computed via dampened Newthon method",
#'      main = "Discrepancy between EL1 and EL0", ylab = "")
#' abline(h = 0, lty = 2)
#'
#' # Handling the convex hull violation differently
#' EL0(1:9)
#' EL0(1:9, log.control = list(order = 2))  # Warning + huge lambda
#' EL0(1:9, log.control = list(order = 4))  # Warning + huge lambda
#'
#' # Warning: depending on the compiler, the discrepancy between EL and EL0
#' # can be one million (1) times larger than the machine epsilon despite both of them
#' # being written in pure R
#' # The results from Apple clang-1400.0.29.202 and Fortran GCC 12.2.0 are different from
#' # those obtained under Ubuntu 22.04.4 + GCC 11.4.0-1ubuntu1~22.04,
#' # Arch Linux 6.6.21 + GCC 14.1.1, and Windows Server 2022 + GCC 13.2.0
#' out0 <- EL0(earth, mu = 5.517, return.weights = TRUE)[1:4]
#' out1 <- EL1(earth, mu = 5.517, return.weights = TRUE)[1:4]
#' print(c(out0$lam, out1$lam), 16)
#'
#' # Value of lambda                                EL0                 EL1
#' # aarch64-apple-darwin20          -1.5631313957?????  -1.5631313955?????
#' # Windows, Ubuntu, Arch           -1.563131395492627  -1.563131395492627
#' @export
EL0 <- function(z, mu = NULL, ct = NULL, shift = NULL, renormalise = FALSE, return.weights = FALSE,
                weight.tolerance = NULL, boundary.tolerance = 1e-9, trunc.to = 0, deriv = FALSE,
                log.control = list(order = NULL, lower = NULL, upper = NULL), verbose = FALSE
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
  if (is.null(shift)) shift <- rep(0, n)
  if (is.null(weight.tolerance)) weight.tolerance <- max(ct, 1) * sqrt(.Machine$double.eps)

  # If originally the weights were too small, too many points would be truncated
  # Warn if any non-zero weights are smaller than weight.tolerance
  closeto0 <- (abs(ct) < weight.tolerance)
  if (any(closeto0)) {
    if (verbose) warning(paste0("Counts closer to 0 than ", sprintf("%1.2e", weight.tolerance),
                                " have been replaced with ", if (trunc.to == 0) "0." else trunc.to))
    ct[closeto0 & ct > 0] <- trunc.to
  }

  if (return.weights) {
    wts <- numeric(n)
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

  # Re-normalising (optionally) only after weight truncation
  if (renormalise) ct <- ct / sum(ct)

  # The denominator for EL with counts is the sum of total counts, and with re-normalisation, it is the number of observations
  n.denom <- if (renormalise) n else sum(ct)
  if (n.denom <= 0) stop("Total weights after tolerance checks (", n.denom, ") must be positive (check the counts and maybe decrease 'weight.tolerance', which is now ", sprintf("%1.1e", weight.tolerance), ".")

  # Enforcing the moment condition
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
  exitcode <- 5L

  spanning <- z1 <= 0 && zn >= 0

  # If logarithm Taylor expansion is requested, prepare good defaults
  if (is.null(log.control$order) || !is.finite(log.control$order) || log.control$order == 0) log.control$order <- NA
  if (is.finite(log.control$order)) {
    if (floor(log.control$order/2) != log.control$order/2) stop("EL0: the extrapolation order must be an positive even number, typically 2 or 4.")
    if (is.null(log.control$lower)) log.control$lower <- if (!renormalise) 1/n else stats::median(ct)/n
    if (is.null(log.control$upper)) log.control$upper <- Inf
  }

  if (n < 2) { # Codes > 5
    exitcode <- 6L
  } else if (z1 == zn) { # The sample is degenerate without variance, no extrapolation possible
    exitcode <- 8L
  } else if (z1 == 0 || zn == 0) {
    # mu is on the boundary -- special case (still in the book)
    logelr <- -Inf
    lam <- iter <- estim.prec <- f.root <- 0
    converged <- TRUE
    int <- c(0, 0)
    exitcode <- 7L
  } else if ((!spanning) && is.na(log.control$order)) {
    # No Taylor extrapolation, no feasible solution
    # Do nothing, this is the default case
  } else {  # The main 'good' loop: EL may proceed because the spanning condition holds
    # Search bracket for lambda
    negz <- z < 0
    comp <- (ct / n.denom - 1 - shift) / z
    min.lam <- suppressWarnings(max(comp[!negz]))
    max.lam <- suppressWarnings(min(comp[negz]))
    if (!is.finite(min.lam)) min.lam <- -max.lam
    if (!is.finite(max.lam)) max.lam <- -min.lam
    int <- c(min.lam, max.lam)
    int <- int + abs(int) * c(2, -2) * .Machine$double.eps # To avoid bad rounding
    ur.ctrl <- list(tol = 1e-16, maxiter = 100, trace = 0)

    # Empirical log-likelihood and its analytical derivative
    llik  <- function(lambda) -sum(ct *     logTaylor(1 + z*lambda + shift, der = 0, order = log.control$order, lower = log.control$lower, upper = log.control$upper))
    dllik <- function(lambda)  sum(ct * z * logTaylor(1 + z*lambda + shift, der = 1, order = log.control$order, lower = log.control$lower, upper = log.control$upper))
    # Simpler: # dllik <- function(lambda) sum(ct * z * dlog(1 + z * lambda + shift, d = 1))

    # xseq <- seq(-0.1, 1, length.out = 301)
    # plot(xseq, logTaylor(xseq, order = 6, lower = 0.2)); abline(v = 0, lty = 3)

    # xs <- seq(int[1], int[2], length.out = 51)
    # ys <- sapply(xs,  llik)
    # dys <- sapply(xs, dllik)
    # plot(xs, ys)
    # plot(xs, dys)

    # Wishing the strictest convergence in the uniroot-like control
    if (spanning) {
      # Search with strict boundaries
      lam.list <- tryCatch(brentZero(dllik, interval = int, tol = .Machine$double.eps, extendInt = "no", maxiter = 100, trace = verbose*2),
                           error = function(e) return(NULL))
    } else {
      # Search with flexible boundaries
      lam.list <- tryCatch(brentZero(dllik, interval = int, tol = .Machine$double.eps, extendInt = "yes", maxiter = 100, trace = verbose*2),
                           error = function(e) return(e))
      exitcode <- 10L
      # The sum of weights is likely to be less than 1
      # If a lambda was not found in enough iterations, keep the wrong solution; the error code should be informative
    }

    # There can be only one kind of warning: maximum iterations reached
    if (!is.null(lam.list)) { # Some result with or without a warning as the second element of the list
      # Preparation for Owen's Taylor strategy
      lam <- lam.list$root
      wvec <- ct * logTaylor(1 + z*lam + shift, der = 1, order = log.control$order, lower = log.control$lower, upper = log.control$upper)
      if (!renormalise) wvec <- wvec / n.denom
      if (return.weights) wts[nonz] <- wvec
      logelr <- llik(lam)
      # Empirical fix for nonsensical probabilities
      # This should not happen unless the spanning condition fails and the Taylor expansion is very inaccurate
      if (any(wvec < 0) && logelr > 0) logelr <- -logelr
      if (any(!is.finite(wvec))) exitcode <- 12L

      # brentZero returns the number of iterations times -1 in case it exceeds the maximum number allowed
      converged <- lam.list$iter  >= 0 && lam.list$exitcode == 0
      estim.prec <- lam.list$estim.prec
      f.root <- lam.list$f.root
      iter <- lam.list$iter
      exitcode <- 0L

      if (abs(f.root) > sqrt(.Machine$double.eps)) exitcode <- 1L
      # Relative tolerance check for boundary closeness
      int.len <- max.lam - min.lam
      if (min(abs(lam - max.lam), abs(lam - min.lam)) < boundary.tolerance * int.len)
        exitcode <- exitcode + 2L
      if (abs(sum(wvec) - 1) > 1e-6 && exitcode != 10) exitcode <- 11L
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

  # Computing analytical derivatives
  if (deriv) {
    u <- 1 + z*lam + shift
    S0 <- sum(ct/u)
    S1 <- sum(ct/u^2)
    T1 <- sum(ct*z/u^2)
    T2 <- sum(ct*z^2/u^2)
    fp <- S0*lam
    fpp <- lam^2*S1 - (S0 - lam*T1)^2/T2
    deriv <- c(fp, fpp)
  } else {
    deriv <- NULL
  }

  return(list(logelr = logelr, lam = lam, wts = wts, converged = converged, iter = iter,
              bracket = int, estim.prec = estim.prec, f.root = f.root, deriv = deriv,
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
#' Implementation note: the EL solver also guarantees a descent direction; if the Newton step is non-descent or non-finite,
#' it falls back to steepest descent (negative gradient), which keeps the line search well-behaved.
#'
#' Tweak \code{alpha} and \code{beta} with extreme caution. See \insertCite{boyd2004convex}{smoothemplik},
#' pp. 464--466 for details. It is necessary that \code{0 < alpha < 1/2} and \code{0 < beta < 1}.
#' \code{alpha = 0.3} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).
#'
#' The argument names, except for \code{lambda.init}, are matching the original names in Art B. Owen's implementation.
#' The highly optimised one-dimensional counterpart, [EL0()], is designed to return a faster
#' and a more accurate solution in the one-dimensional case.
#'
#' @param z A numeric vector or a matrix with one data vector per column.
#' @param mu Hypothesised mean, default \code{(0 ... 0)} in \eqn{R^{\mathrm{ncol}(z)}}{R^ncol(z)}.
#' @param ct Numeric count variable with non-negative values that indicates the multiplicity of observations.
#' @param shift The value to add in the denominator (useful in case there are extra
#'   Lagrange multipliers): \eqn{1 + \lambda'Z + shift}{1 + lambda'Z + shift}.
#' @param lambda.init Starting lambda, default \code{(0 ... 0)}. Improves speed and accuracy
#'   in sequential problems if supplied from the previous iteration.
#' @param renormalise If \code{FALSE}, then uses the total sum of counts as the number of observations,
#'   like in vanilla empirical likelihood, due to formula (2.9) in \insertCite{owen2001empirical}{smoothemplik},
#'   otherwise re-normalises the counts to 1 according to \insertCite{cosma2019inference}{smoothemplik}
#'   (p. 170, the topmost formula).
#' @param return.weights Logical: if \code{TRUE}, returns the empirical probabilities.
#'   Default is memory-saving (\code{FALSE}).
#' @param lower Lower cut-off for [logTaylor()], default \code{1/NROW(z)}.
#' @param upper Upper cut-off for [logTaylor()], default \code{Inf}.
#' @param order Positive even integer such that the Taylor approximation of this order to
#'   \eqn{\log x}{log(x)} is self-concordant; usually 4 or 2. Passed to [logTaylor()].
#' @param weight.tolerance Weight tolerance for counts to improve numerical stability
#'   (defaults to \code{sqrt(.Machine$double.eps)} times the maximum weight).
#' @param deriv Logical: if \code{TRUE}, computes and returns the first two directional derivatives
#'   of log-ELR w.r.t. \code{mu} in the direction of the hypothesised value.
#' @param thresh Target tolerance on the squared Newton decrement: loop stops when \code{decr^2 <= thresh}.
#'   (If \code{verbose} is \code{TRUE}, decrement itself is printed.)
#' @param itermax Maximum number of outer iterations of the damped Newton method (seems ample).
#' @param alpha Backtracking line search Armijo parameter: acceptance of a decrease in function value
#'   by \eqn{\alpha f}{ALPHA*f} of the prediction based on the linear extrapolation. Smaller makes acceptance easier.
#' @param beta Backtracking step shrinkage factor in \code{[0, 1]}. 0.1 corresponds to a very crude search,
#'   0.8 corresponds to a less crude search.
#' @param backeps Backtrack threshold, a small slack added to Armijo RHS: the search can miss by this much.
#'   Accept if \eqn{f(x+tp) \le f(x)+\alpha t g'p + \mathrm{backeps}}{f(x+tp) <= f(x) + alpha\*t\*g'p + backeps}.
#'   Consider setting it to \code{1e-10} if backtracking seems to be failing due to round-off.
#' @param gradtol Gradient tolerance: stop if \code{||g|| <= gradtol}.
#' @param steptol Step tolerance: stop if the relative size is tiny: \code{||x2-x1||/max(1, ||x2||) < ftol}.
#' @param ftol Function change tolerance: stop if the relative function-value change is less than \code{ftol}.
#' @param stallmax Stop if both \code{rel_step <= steptol} and \code{rel_f <= ftol} hold for this many consecutive iterations.
#' @param verbose Logical: print output diagnostics?
#'
#' @return A list with the following values:
#' \describe{
#'     \item{logelr}{Log of empirical likelihood ratio (equal to 0 if the hypothesised mean is equal to the sample mean)}
#'     \item{lam}{Vector of Lagrange multipliers}
#'     \item{wts}{Observation weights/probabilities (vector of length n)}
#'     \item{deriv}{Length-2 vector: directional first and second derivatives along the ray toward mu (if \code{deriv = TRUE})}
#'     \item{converged}{\code{TRUE} if algorithm converged. \code{FALSE} usually means that mu is not in the convex hull of the data. Then, a very small likelihood is returned (instead of zero).}
#'     \item{iter}{Number of iterations taken.}
#'     \item{ndec}{Newton decrement (see Boyd & Vandenberghe).}
#'     \item{gradnorm}{Norm of the gradient of log empirical likelihood.}
#' }
#'
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
#' EL1(earth, mu = 5.517, verbose = TRUE) # 5.517 is the modern accepted value
#'
#' # Linear regression through empirical likelihood
#' coef.lm <- coef(lm(mpg ~ hp + am, data = mtcars))
#' xmat <- cbind(1, as.matrix(mtcars[, c("hp", "am")]))
#' yvec <- mtcars$mpg
#' foc.lm <- function(par, x, y) {  # The sample average of this
#'   resid <- y - drop(x %*% par)   # must be 0
#'   resid * x
#' }
#' minusEL <- function(par) -EL1(foc.lm(par, xmat, yvec), itermax = 10)$logelr
#' coef.el <- optim(c(26, -0.06, 5.3), minusEL, control = list(maxit = 100))$par
#' abs(coef.el - coef.lm) / coef.lm  # Relative difference
#'
#' # Likelihood ratio testing without any variance estimation
#' # Define the profile empirical likelihood for the coefficient on am
#' minusPEL <- function(par.free, par.am)
#'   -EL1(foc.lm(c(par.free, par.am), xmat, yvec), itermax = 20)$logelr
#' # Constrained maximisation assuming that the coef on par.am is 3.14
#' coef.el.constr <- optim(coef.el[1:2], minusPEL, par.am = 3.14)$par
#' print(-2 * EL1(foc.lm(c(coef.el.constr, 3.14), xmat, yvec))$logelr)
#' # Exceeds the critical value qchisq(0.95, df = 1)
#' @export
EL1 <- function(z, mu = NULL, ct = NULL, shift = NULL, lambda.init = NULL, renormalise = FALSE,
                return.weights = FALSE, lower = NULL, upper = NULL,
                order = NA, weight.tolerance = NULL, deriv = FALSE,
                thresh = 1e-30, itermax = 100L, verbose = FALSE,
                alpha = 0.3, beta = 0.8, backeps = 0,
                gradtol = 1e-12, steptol = 1e-12, ftol = 1e-14, stallmax = 5) {
  if (is.null(dim(z))) z <- matrix(z, ncol = 1)
  n <- nrow(z)
  d <- ncol(z)
  if (length(mu) == 0) mu <- rep(0, d)
  if (is.null(shift)) shift <- 0
  if (length(mu) != d) stop("The length of mu must be the same as the dimension of z.")
  if (length(lambda.init) == 0) lambda.init <- rep(0, d)
  if (length(lambda.init) != d) stop("The length of mu must be the same as the dimension of z.")

  if (is.null(lower)) lower <- rep(1/n, n)
  if (is.null(upper)) upper <- rep(Inf, n)

  if (is.null(weight.tolerance)) weight.tolerance <- max(ct, 1) * sqrt(.Machine$double.eps)

  if (is.null(ct)) ct <- rep(1, n)
  if (min(ct) < 0) stop("Negative weights are not allowed.")
  if (renormalise) ct <- ct / sum(ct)
  if (any(0 < ct & ct < weight.tolerance)) {
    if (verbose) warning(paste("Positive counts below", weight.tolerance, "have been replaced by zero."))
    ct[ct < weight.tolerance] <- 0
    if (renormalise) ct <- ct / sum(ct)  # Re-normalising again
  }
  if (sum(ct) <= 0) stop("Total weight must be positive.")

  ret <- ELCPP(z = z, ct = ct, mu = mu, shift = shift, lambda_init = lambda.init,
               return_weights = return.weights, lower = lower, upper = upper,
               order = order, weight_tolerance = weight.tolerance, deriv = deriv,
               thresh = thresh, itermax = itermax, verbose = verbose,
               alpha = alpha, beta = beta, backeps = backeps, grad_tol = gradtol,
               step_tol = steptol, f_tol = ftol, stallmax = stallmax)
  if (ret$exitcode > 0) ret$logelr <- -Inf
  return(ret)
}

#' Compute empirical likelihood on a trajectory
#'
#' @param z Passed to [EL1()].
#' @param ct Passed to [EL1()].
#' @param mu0 Starting point of trajectory
#' @param mu1 End point of trajectory
#' @param N Number of segments into which the path is split (i. e. \code{N+1} steps are used).
#' @param order Passed to [EL1()]. It is highly advised to avoid using \code{NA}
#'   (no extrapolation) because the lambda search may fail with unmodified logarithm.
#' @param verbose Logical: report iteration progress?
#' @param ... Passed to [EL1()].
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
#' EL1(earth, mu = 5.1,  verbose = TRUE)
#' logELR <- ctracelr(earth, mu0 = 5.1, mu1 = 5.65, N = 55, verbose = TRUE)
#' hist(earth, breaks = seq(4.75, 6, 1/8))
#' plot(logELR[, 1], exp(logELR[, 2]), bty = "n", type = "l",
#'      xlab = "Earth density", ylab = "ELR")
#'
#' # Two-dimensional trajectory
#' set.seed(1)
#' xy <- matrix(rexp(200), ncol = 2)
#' logELR2 <- ctracelr(xy, mu0 = c(0.5, 0.5), mu1 = c(1.5, 1.5), N = 100)
ctracelr <- function(z, ct = NULL, mu0, mu1, N = 5, order = 4, verbose = FALSE, ...) {
  if (is.vector(z)) z <- matrix(z, ncol = 1)
  d <- ncol(z)

  lam0 <- NULL
  elr <- ec <- it <- nd <- gn <- numeric(N+1)
  m <- l <- matrix(NA, N+1, ncol = d)
  for (i in 0:N) {
    mui <- (i*mu1 + (N-i)*mu0) / N
    if (i > 0) lam0 <- if (all(is.finite(x$lam))) drop(x$lam) else NULL
    x <- EL1(z = z, ct = ct, mu = mui, lambda.init = lam0, order = order, ...)
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
#'
#' @details
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
#' @return A list with the same structure as that in [EL1()].
#' @seealso [EL1()]
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
                weight.tolerance = NULL, trunc.to = 0, renormalise = TRUE,
                return.weights = FALSE, verbose = FALSE
) {
  if (is.null(dim(z)) || is.data.frame(z)) z <- as.matrix(z, rownames.force = TRUE)
  n <- nrow(z)
  k <- ncol(z)
  if (is.null(mu)) mu <- rep(0, k)
  if (is.null(ct)) ct <- rep(1, n)
  if (is.null(vt)) vt <- rep(1, n)
  if (is.null(shift)) shift <- rep(0, n)
  n.orig <- n
  if (is.null(weight.tolerance)) weight.tolerance <- max(ct, 1) * sqrt(.Machine$double.eps)
  ret <- EuLCPP(z = z, mu = mu, ct = ct, vt = vt, shift = shift, n_orig = n.orig,
                weight_tolerance = weight.tolerance, trunc_to = trunc.to, renormalise = renormalise,
                return_weights = return.weights, verbose = verbose)
  if (return.weights && !is.null(nz <- rownames(z))) names(ret$wts) <- nz
  return(ret)
}


#' Unified empirical likelihood wrapper
#'
#' @description
#' Call \code{EL0()}, \code{EL1()}, or \code{EuL()} through a single interface.
#' If extrapolation is requested, switch to dedicated functions.
#' Anything method-specific goes into \code{EL.args}.
#'
#' @param type Character: one of \code{c("auto", "EL1", "EL0", "EuL")}. If \code{"auto"},
#'   uses \code{"EL1"} for multi-variate data and \code{"EL0"} for uni-variate.
#' @inheritParams EL1
#' @param chull.fail Character: \code{"none"} calls the original EL (which may return
#'   \code{-Inf} in case of a convex-hull violation), \code{"taylor"} calls [ExEL1()],
#'   \code{"wald"} calls [ExEL2()], \code{"adjusted"} adds one pseudo-observation as in
#'   \insertCite{chen2008adjusted}{smoothemplik}, \code{"adjusted2"} adds one (in 1D) or
#'   two (2D+) pseudo-observations with improved coverage rate according to
#'   \insertCite{liu2010adjusted}{smoothemplik}, and \code{"balanced"} adds two
#'   pseudo-observations according to \insertCite{emerson2009calibration}{smoothemplik}.
#' @param ... Named extra arguments passed to the selected back-end (e.g. \code{order},
#'   \code{itermax}, \code{lambda.init}, \code{vt}, \code{trunc.to}, \code{boundary.tolerance}, ...).
#'
#' @return A list with either the return value of the selected back-end  or (for extrapolation
#'   methods) at least the \code{logelr} list value and extrapolation attributes.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # EL0 with extras:
#' EL(type = "EL0", z = 1:9, mu = 4, boundary.tolerance = 1e-8)
#' # EL1 with a custom order and iteration cap:
#' set.seed(1)
#' x <- cbind(rnorm(30), runif(30)-0.5)
#' EL(type = "EL1", z = x, mu = c(0, 0), order = 4, itermax = 50, return.weights = TRUE)
#' # EuL with vt and truncation:
#' set.seed(1)
#' EL(type = "EuL", z = x, vt = runif(NROW(x)), weight.tolerance = 0.1, trunc.to = 0.1)
#'
#' # Extrapolated variants
#' set.seed(1)
#' EL(type = "EL0", z = 1:9, mu = 12, chull.fail = "taylor", exel.control = list(xlim = c(2, 8)))
#' EL(type = "EL1", z = 1:9, mu = 12, chull.fail = "wald", exel.control = list(fmax = 10))
#' x <- matrix(runif(20), ncol = 2)
#' EL(x, mu = c(0, 0), chull.fail = "adjusted")
#' EL(x, mu = c(0, 0), chull.fail = "adjusted2")
#' EL(x, mu = c(0, 0), chull.fail = "balanced")
#'
#' @export
EL <- function(z, ct = NULL, mu = NULL, shift = NULL,
               type = c("auto", "EL1", "EL0", "EuL"),
               chull.fail = c("none", "taylor", "wald", "adjusted", "adjusted2", "balanced"),
               renormalise = FALSE, return.weights = FALSE, weight.tolerance = NULL,
               verbose = FALSE, ...) {

  type <- match.arg(type)
  if (type == "auto") type <- if (NCOL(z) == 1) "EL0" else "EL1"
  chull.fail <-  match.arg(chull.fail)
  chull.fail <- tolower(chull.fail)  # Preventing Wald and Taylor

  if (type == "EuL" & !renormalise) renormalise <- TRUE

  if (is.null(ct)) ct <- rep(1, NROW(z))
  zct <- ct == 0
  if (any(zct)) {
    z <- if (!is.null(dim(z))) z[!zct, , drop = FALSE] else z[!zct]
    ct <- ct[!zct]
    if (!is.null(shift)) shift <- shift[!zct]
  }

  # Extrapolated EL requires a default mu -- providing 0
  if (chull.fail %in% c("taylor", "wald") && is.null(mu)) mu <- rep(0, NCOL(z))

  # Common args shared by all back-ends, without the NULLs
  base_args <- list(z = z, ct = ct, mu = mu, shift = shift, renormalise = renormalise,
    return.weights = return.weights, weight.tolerance = weight.tolerance, verbose = verbose)
  base_args <- base_args[!vapply(base_args, is.null, logical(1))]

  # Collect method-specific extras from the ellipsis
  dots <- list(...)
  if (length(dots) > 0) {
    if (is.null(names(dots)) || any(names(dots) == "")) stop("All extra arguments in `...` must be named.")
  }
  call.args <- c(base_args, dots)


  if (chull.fail == "none") {
    res <- switch(
      type,
      EL1 = do.call(EL1, call.args),
      EL0 = do.call(EL0, call.args),
      EuL = do.call(EuL, call.args)
    )
  } else if (chull.fail == "taylor") {
    call.args$type <- type
    res <- switch(
      type,
      EL1 = do.call(ExEL1, call.args),
      EL0 = do.call(ExEL1, call.args),
      EuL = stop("Euclidean likelihood needs no extrapolation. Call it without the 'chull.fail' argument.")
    )
    res <- list(logelr = res)
  } else if (chull.fail == "wald") {
    call.args$type <- type
    res <- switch(
      type,
      EL1 = do.call(ExEL2, call.args),
      EL0 = do.call(ExEL2, call.args),
      EuL = stop("Euclidean likelihood needs no extrapolation. Call it without the 'chull.fail' argument.")
    )
    res <- list(logelr = res)
  } else if (chull.fail == "adjusted") {
    if (is.null(dim(z))) z <- as.matrix(z)
    an <- max(1, log(nrow(z))/2)
    zm <- apply(z, 2, stats::weighted.mean, w = ct)
    point1 <- -zm * an
    call.args$z <- rbind(z, point1)
    call.args$ct <- c(ct, mean(ct))
    res <- switch(
      type,
      EL1 = do.call(EL1, call.args),
      EL0 = do.call(EL0, call.args),
      EuL = stop("Euclidean likelihood needs no extrapolation. Call it without the 'chull.fail' argument.")
    )
    attr(res, "point1") <- point1
  } else if (chull.fail == "adjusted2") {
    if (is.null(dim(z))) z <- as.matrix(z)
    b  <- bartlettFactor(z)
    b1 <- attr(b, "components")[1]
    b2 <- attr(b, "components")[2]
    zm <- apply(z, 2, stats::weighted.mean, w = ct)
    if (NCOL(z) == 1) {
      point1 <- -zm*b/2
      point2 <- NULL
    } else {
      point1 <- -zm*b1/2
      point2 <-  zm*b2/2
    }
    call.args$z <- rbind(z, point1, point2)
    call.args$ct <- if (NCOL(z) == 1) c(ct, mean(ct)) else c(ct, rep(mean(ct), 2))
    res <- switch(
      type,
      EL1 = do.call(EL1, call.args),
      EL0 = do.call(EL0, call.args),
      EuL = stop("Euclidean likelihood needs no extrapolation. Call it without the 'chull.fail' argument.")
    )
    attr(res, "point1") <- point1
    attr(res, "point2") <- point2
  } else if (chull.fail == "balanced") {  # Emerson & Owen 2009
    # TODO: NOT TESTED
    if (is.null(dim(z))) z <- as.matrix(z)
    zbar <- apply(z, 2, mean, trim = 0.05)
    zm <- colMeans(z)
    zl <- sqrt(sum(zm^2))  # Length of the mean of z
    V <- stats::var(z)
    u <- zm / zl
    cu <- 1 / drop(sqrt(t(u) %*% solve(V) %*% u))
    s <- 1.6  # From Emerson & Owen (2009)
    delta <- s*cu / zl
    point1 <- -zbar * delta
    point2 <- (2+delta)*zbar
    call.args$z <- rbind(z, point1, point2)
    call.args$ct <- c(ct, rep(mean(ct), 2))
    res <- switch(
      type,
      EL1 = do.call(EL1, call.args),
      EL0 = do.call(EL0, call.args),
      EuL = stop("Euclidean likelihood needs no extrapolation. Call it without the 'chull.fail' argument.")
    )
    attr(res, "point1") <- point1
    attr(res, "point2") <- point2
  }


  attr(res, "method") <- type
  res
}

