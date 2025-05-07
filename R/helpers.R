# Interpolate monotonically from a lower to a higher parabola function
# defined by the equation (x - mean)^2 / var
# after a certain point over an interval of the chosen length
interpToHigher <- function(x, f, mean, var, at, gap) {
  fa <- function(x) (x-mean)^2/var
  incr <- at > mean
  if (incr) {
    xleft  <- at + c(0, 0.05, 0.1)*gap
    xright <- at + c(0.9, 0.95, 1)*gap
    if (all(x >= max(xright))) return(fa(x))  # Speed-up: far away = extrapolate for sure
    fleft <- f(xleft)
    fright <- fa(xright)
    s <- stats::splinefun(c(xleft, xright), c(fleft, fright), method = "monoH.FC")
    ifelse(x < at, f(x), ifelse(x < at + gap, s(x), fa(x)))
  } else {
    xleft  <- at - c(0.9, 0.95, 1)*gap
    xright <- at - c(0.1, 0.05, 0)*gap
    if (all(x <= min(xleft))) return(fa(x))
    fleft <- fa(xleft)
    fright <- f(xright)
    s <- stats::splinefun(c(xleft, xright), c(fleft, fright), method = "monoH.FC")
    ifelse(x > at, f(x), ifelse(x > at - gap, s(x), fa(x)))
  }
}

interpToLower <- function(x, f, mean, var, at, gap) {
  fa <- function(x) (x-mean)^2/var
  incr <- at > mean
  if (incr) {
    xleft  <- at + c(0, 0.01, 0.02)*gap
    at2 <- sqrt(f(at) * var) + mean
    xright <- at2 + c(0.98, 0.99, 1)*gap
    if (all(x >= max(xright))) return(fa(x))
    fleft <- f(xleft)
    fright <- fa(xright)
    s <- stats::splinefun(c(xleft, xright), c(fleft, fright), method = "monoH.FC")
    ifelse(x < at, f(x), ifelse(x < at2 + gap, s(x), fa(x)))
  } else {
    xright <- at - c(0.02, 0.01, 0)*gap
    at2 <- -sqrt(f(at) * var) + mean
    xleft  <- at2 - c(0.98, 0.99, 1)*gap
    if (all(x <= min(xleft))) return(fa(x))
    fleft <- fa(xleft)
    fright <- f(xright)
    s <- stats::splinefun(c(xleft, xright), c(fleft, fright), method = "monoH.FC")
    ifelse(x > at, f(x), ifelse(x > at2 - gap, s(x), fa(x)))
  }
}

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
#' `splinefun(..., method = "monoH.FC")`.  Anchor points on each side of the
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
#' f <- Vectorize(function(m) -2*weightedEL(xx, mu = m, ct = w, chull.fail = "none")$logelr)
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
    interpToHigher(x, f, mean, var, at, gap)
  else
    interpToLower(x, f, mean, var, at, gap)
}
