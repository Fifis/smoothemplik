#' Brent's local minimisation
#'
#' @param f A function to be minimised on an interval.
#' @param interval A length-2 vector containing the end-points of the search interval.
#' @param lower Scalar: the lower end point of the search interval. Not necessary if \code{interval} is provided.
#' @param upper Scalar: the upper end point of the search interval. Not necessary if \code{interval} is provided.
#' @param tol Small positive scalar: stopping criterion. The search stops when the
#'   distance between the current candidate and the midpoint of the bracket is smaller than
#'   the dynamic threshold \code{2 * (sqrt(DBL_EPSILON) * abs(x) + tol)}
#' @param maxiter Positive integer: the maximum number of iterations.
#' @param trace Integer: 0, 1, or 2. Amount of tracing information on the optimisation progress
#'   printed. \code{trace = 0} produces no output, \code{trace = 1} reports the starting
#'   and final results, and \code{trace = 2} provides detailed iteration-level output.
#'
#' @details
#' This is an adaptation of the implementation by John Burkardt (currently available at
#' [https://people.math.sc.edu/Burkardt/m_src/brent/brent.html](https://people.math.sc.edu/Burkardt/m_src/brent/brent.html)).
#'
#' This function is similar to \code{local_min} or \code{R_zeroin2}-style logic, but with the
#' following additions: the number of iterations is tracked, and the algorithm stops when the
#' standard Brent criterion is met or if the maximum iteration count is reached.
#' The code stores the approximate final bracket width in \code{estim.prec}, like in [uniroot()].
#' If the minimiser is pinned to an end point, \code{estim.prec = NA}.
#'
#' There are no preliminary iterations, unlike [brentZero()].
#'
#' TODO: add preliminary iterations.
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{root}{Location of the minimum.}
#'   \item{f.root}{Function value at the minimuim location.}
#'   \item{iter}{Total iteration count used.}
#'   \item{estim.prec}{Estimate of the final bracket size.}
#' }
#'
#' @export
#'
#' @examples
#' f <- function (x) (x - 1/3)^2
#' brentMin(f, c(0, 1), tol = 0.0001)
#' brentMin(function(x) x^2*(x-1), lower = 0, upper = 10, trace = 1)
brentMin <- function(f, interval, lower = NA_real_, upper = NA_real_, tol = 1e-8,
                     maxiter = 200L, trace = 0L) {
  if (missing(interval)) interval <- c(lower, upper)
  brentMinCPP(f = f, interval = interval, lower = lower, upper = upper,
              tol = tol, maxiter = maxiter, trace = trace)
}

#' Brent's local root search with extended capabilities
#'
#' @param f The function for which the root is sought.
#' @param interval A length-2 vector containing the end-points of the search interval
#' @param lower Scalar: the lower end point of the search interval. Not necessary if \code{interval} is provided.
#' @param upper Scalar: the upper end point of the search interval. Not necessary if \code{interval} is provided.
#' @param f_lower Scalar: same as f(upper). Passing this value saves time if f(lower) is slow to compute and is known.
#' @param f_upper Scalar: same as f(lower).
#' @param extendInt Character:
#'   \describe{
#'     \item{\code{"no"}}{Do not extend the interval (default).}
#'     \item{\code{"yes"}}{Attempt to extend both ends until a sign change is found.}
#'     \item{\code{"upX"}}{Assumes the function is increasing around the root and extends upward if needed.}
#'     \item{\code{"downX"}}{Assumes the function is decreasing around the root and extends downward if needed.}
#'     \item{\code{"right"}}{Attempt to extend the upper (right) end until a sign change is found.}
#'     \item{\code{"left"}}{Attempt to extend the lower (left) end until a sign change is found.}
#'   }
#'   This behavior mirrors that of [uniroot()].
#' @param tol Small positive scalar: convergence tolerance. The search stops when the bracket size is smaller than
#'   \code{2 * .Machine$double.eps * abs(x) + tol}, or if the function evaluates to zero at the candidate root.
#' @param maxiter Positive integer: the maximum number of iterations before stopping.
#' @param trace Integer: 0, 1, or 2. Controls the verbosity of the output.
#'   \code{trace = 0} produces no output, \code{trace = 1} reports the starting and final results,
#'   and \code{trace = 2} provides detailed iteration-level output.
#'
#'
#' @returns A list with the following elements:
#' \describe{
#'   \item{root}{Location of the root.}
#'   \item{f.root}{Function value at the root.}
#'   \item{iter}{Total iteration count used.}
#'   \item{init.it}{Number of initial \code{extendInt} iterations if there were any; NA otherwise.}
#'   \item{estim.prec}{Estimate of the final bracket size.}
#'   \item{exitcode}{0 for success, 1 for maximum initial iteration limit, 2 for maximum main iteration limit.}
#' }
#' @export
#'
#' @examples
#' f <- function (x, a) x - a
#' str(uniroot(f, c(0, 1), tol = 0.0001, a = 1/3))
#' uniroot(function(x) cos(x) - x, lower = -pi, upper = pi, tol = 1e-9)$root
#'
#' # New capabilities: extending only one end of the interval
#' f <- function(x) x^2 - 1  # The roots are -1 and 1
#' brentZero(f, c(2, 3), extendInt = "left")
#' brentZero(f, c(2, 3), extendInt = "yes")
#' brentZero(f, c(2, 3), extendInt = "upX")
#' brentZero(f, c(0, 0.5), extendInt = "downX")  # This one finds the left crossing
#'
#' # This function is faster than the base R uniroot, and this is the primary
#' # reason why it was written in C++
#' system.time(replicate(1000, { shift <- runif(1, 0, 2*pi)
#'   uniroot(function(x) cos(x+shift) - x, lower = -pi, upper = pi)
#' }))
#' system.time(replicate(1000, { shift <- runif(1, 0, 2*pi)
#'   brentZero(function(x) cos(x+shift) - x, lower = -pi, upper = pi)
#' }))
#' # Roughly twice as fast
brentZero <-  function(f, interval, lower = NA_real_, upper = NA_real_,
                       f_lower = NULL, f_upper = NULL, extendInt = "no",
                       tol = 1e-8, maxiter = 500L, trace = 0L) {
  if (missing(interval)) interval <- c(lower, upper)
  brentZeroCPP(f = f, interval = interval, lower = lower, upper = upper,
               f_lower = f_lower, f_upper = f_upper, extendInt = extendInt,
               tol = tol, maxiter = maxiter, trace = trace)
}


