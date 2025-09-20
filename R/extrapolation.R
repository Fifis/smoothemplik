#' Extrapolated EL of the first kind (Taylor expansion)
#'
#' @param z Passed to \code{EL0}/\code{EL1}.
#' @param mu Passed to \code{EL0}/\code{EL1}.
#' @param type If \code{"EL0"}, uses uni-variate [EL0()] for calculations; same for \code{"EL1"}.
#' @param exel.control A list with the following elements: \code{xlim} -- if \code{"auto"}, uses a quick boundary detection,
#'   otherwise should be a length-two numeric vector; \code{fmax} -- maximum allowed chi-squared statistic value for a thorough
#'   root search with probability \code{p} and degrees of freedom \code{df}.
#' @param ... Also passed to \code{EL0}/\code{EL1}.
#'
#' @returns A numeric vector of log-ELR statistic of the same length as \code{mu}.
#' @export
#'
#' @examples
#' z <- c(1, 4, 5, 5, 6, 6)
#' ExEL1(z, 0.5, ct = 1:6)
#'
#' xseq <- seq(0, 7, 0.2)
#' plot(xseq, -2*ExEL1(z, mu = xseq, ct = 1:6))
#' abline(v = c(1.2, 5.8), h = qchisq(0.99, 1), lty = 3)
#'
#' # User-defined 'good' interval
#' ctrl0 <- list(xlim = c(-1, 8)); ctrl1 <- list(xlim = c(2.5, 5.5))
#' plot(xseq, -2*ExEL1(z, xseq, ct = 1:6, exel.control = ctrl0), bty = "n")
#' lines(xseq, -2*ExEL1(z, xseq, ct = 1:6, exel.control = ctrl1), col = 3)
#' abline(v = ctrl1$xlim, lty = 3)
#'
#' # Root searching
#' ctrl2 <- list(fmax = qchisq(0.99, 1))
#' plot(xseq, -2*ExEL1(z, xseq, ct = 1:6, exel.control = ctrl0), bty = "n")
#' lines(xseq, -2*ExEL1(z, xseq, ct = 1:6, exel.control = ctrl2), col = 3)
#' abline(h = qchisq(0.99, 1), lty = 3)
#'
#' # With EL1 vs. EL0 -- very little discrepancy
#' xseq <- seq(0.8, 1.4, length.out = 101)
#' plot(xseq, -2*ExEL1(z, xseq, ct = 1:6, exel.control = ctrl0), bty = "n")
#' lines(xseq, -2*ExEL1(z, xseq, ct = 1:6, type = "EL0"), col = 3)
#' lines(xseq, -2*ExEL1(z, xseq, ct = 1:6, type = "EL1"), col = 2, lty = 2, lwd = 2)
ExEL1 <- function(z, mu, type = c("EL0", "EL1"), exel.control = list(xlim = "auto", fmax = NA, p = 0.999, df = 1), ...) {
  type <- match.arg(type)
  ell <- list(z = z, ...)
  ct <- ell$ct
  if (is.null(ct)) ct <- ell$ct <- rep(1, length(z))
  ell$mu <- NULL
  ell$return.weights <- FALSE  # To save memory
  ell$renormalise <- FALSE
  ell$deriv <- NULL
  fl <- switch(type,
               EL0 = function(m) do.call(EL0, c(mu = m, deriv = TRUE, ell))[c("logelr", "lam", "deriv")],
               EL1 = function(m) do.call(EL1, c(mu = m, deriv = TRUE, ell))[c("logelr", "lam", "deriv")])
  fl0 <- switch(type,  # Economical version that does not compute derivatives
               EL0 = function(m) do.call(EL0, c(mu = m, deriv = FALSE, ell))[c("logelr", "lam", "deriv")],
               EL1 = function(m) do.call(EL1, c(mu = m, deriv = FALSE, ell))[c("logelr", "lam", "deriv")])
  n <- length(z)
  wm <- stats::weighted.mean(z, ct)  # Used in determining which end to process

  # Evaluation range to determine which mu are acceptable
  # A value must be computed for all cases
  zu <- sort(unique(z))
  nu <- length(zu)
  dzu <- diff(zu)
  if (nu >= 10) {
    rgap <- stats::median(dzu)
  } else if (nu >= 5) {
    rgap <- stats::median(dzu) * ((n-4.5)/5)  # Starting at 0.1 median, ending at 0.9 median
  } else if (nu == 2) {
    rgap <- dzu*0.05  # Same as 0.05 median
  } else if (nu == 3) {
    rgap <- mean(dzu)*0.1  # Mean = median for 2 gaps
  } else if (nu == 4) {
    rgap <- mean(sort(dzu)[-1])*0.1  # Trimmed mean
  } else {
    stop("There must be at least two unique observations for extrapolation.")
  }
  xmax0 <- zu[length(zu)] - rgap
  xmin0 <- zu[1] + rgap
  if (xmin0 >= xmax0) stop("ExEL1: wrong limit order (left ", xmin0, ", right ", xmax0,
                           "). This should never be the case; please report this bug.")

  if (identical(exel.control$xlim, "auto") && is.na(exel.control$fmax)) {
    # No user xlim supplied, no user fmax supplied
    xmin <- xmin0
    xmax <- xmax0
  } else if (is.numeric(exel.control$xlim)) {
    xmin <- min(exel.control$xlim)
    xmax <- max(exel.control$xlim)
  } else {  # xmax and xmin must be searched for
    fmax <- exel.control$fmax
    if (!is.finite(fmax)) fmax <- stats::qchisq(exel.control$p, df = exel.control$df)
    # Fast and slightly inaccurate defaults for a quick search
    xmax <- if (any(mu >= wm)) brentZero(function(x) -2*fl0(x)$logelr - fmax, sort(c(wm, xmax0)),  extendInt = "upX", maxiter = 50)$root else Inf
    xmin <- if (any(mu <= wm)) brentZero(function(x) -2*fl0(x)$logelr - fmax, sort(c(xmin0, wm)), extendInt = "downX", maxiter = 50)$root else -Inf
  }

  i.left  <- mu < xmin
  i.right <- mu > xmax
  i.mid   <- !(i.left | i.right)

  # Analytical derivatives require lambda as well
  fleft <- fmid <- fright <- aleft <- aright <- NULL
  if (any(i.mid)) {
    fmidlist <- lapply(mu[i.mid], fl0)
    fmid <- unlist(lapply(fmidlist, "[[", "logelr"))
  }
  if (any(i.left)) {
    fminlist <- fl(xmin)
    aleft <- getParabola(xmin, fminlist$logelr, fminlist$deriv[1], fminlist$deriv[2])
    fleft <- aleft[1]*mu[i.left]^2 + aleft[2]*mu[i.left] + aleft[3]
  }
  if (any(i.right)) {
    fmaxlist <- fl(xmax)
    aright <- getParabola(xmax, fmaxlist$logelr, fmaxlist$deriv[1], fmaxlist$deriv[2])
    fright <- aright[1]*mu[i.right]^2 + aright[2]*mu[i.right] + aright[3]
  }
  logelr <- c(fleft, fmid, fright)

  attr(logelr, "xlim") <- c(xmin, xmax)
  attr(logelr, "parabola.left") <- aleft
  attr(logelr, "parabola.right") <- aright
  logelr
}

