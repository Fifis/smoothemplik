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
#'
#' # Compating ExEL2 vs ExEL1
#' xseq <- seq(-7, 10.5, 0.1)
#' xl <- range(xseq)
#' # No extrapolation
#' a0 <- ExEL1(-4:4, mu = xseq, ct = 1:9, exel.control = list(xlim = xl))
#' a1 <- ExEL1(-4:4, mu = xseq, ct = 1:9)
#' a2 <- ExEL2(-4:4, mu = xseq, ct = 1:9)
#' v1 <- attr(a1, "xlim")
#' v2 <- c(attr(a2, "bridge.left")[c("x1", "x2")], attr(a2, "bridge.right")[c("x1", "x2")])
#'
#' plot(xseq, a0, ylim = c(-300, 0), xlim = xl, main = "ExEL splices",
#'   bty = "n", xlab = "mu", ylab = "logELR(mu)")
#' lines(xseq, a1, col = 2, lwd = 2)
#' lines(xseq, a2, col = 4, lwd = 2)
#' abline(v = v2, lty = 3)
#' lines(xseq, attr(a2, "parabola.coef") * (xseq - attr(a2, "parabola.centre"))^2, lty = 2)
#' legend("topright", c("Taylor", "Wald", "ax^2"), col = c(2, 4, 1), lwd = c(2, 2, 1), lty = c(1, 1, 2))
#'
#' plot(xseq[-1], diff(a1)/diff(xseq[1:2]), col = 2, type = "l", lwd = 2,
#'   main = "Derivatives of ExEL splice", bty = "n", ylim = c(-100, 100),
#'   xlab = "mu", ylab = "d/dmu logELR(mu)")
#' lines(xseq[-1], diff(a2)/diff(xseq[1:2]), col = 4, lwd = 2)
#' abline(v = c(v1, v2), lty = 3, col = "#00000055")
#' legend("topright", c("Taylor", "Wald"), col = c(2, 4), lwd = 2)
ExEL1 <- function(z, mu, type = c("EL0", "EL1"),
                  exel.control = list(xlim = "auto", fmax = NA, p = 0.999, df = 1), ...) {
  type <- match.arg(type)
  ell <- list(z = z, ...)
  ct <- ell$ct
  if (is.null(ct)) ct <- ell$ct <- rep(1, length(z))
  ell$mu <- NULL
  ell$return.weights <- FALSE  # To save memory
  ell$renormalise <- FALSE
  ell$deriv <- NULL
  fl <- switch(type,
               EL0 = function(m) do.call(EL0, c(mu = m, deriv = TRUE, ell))[c("logelr", "deriv")],
               EL1 = function(m) do.call(EL1, c(mu = m, deriv = TRUE, ell))[c("logelr", "deriv")])
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

#' @rdname ExEL1
#' @export
ExEL2 <- function(z, mu, type = c("EL0", "EL1"),
                  exel.control = list(xlim = "auto", fmax = NA, p = 0.999, df = 1), ...) {
  type <- match.arg(type)
  ell <- list(z = z, ...)
  ct <- ell$ct
  if (is.null(ct)) ct <- ell$ct <- rep(1, length(z))
  ell$mu <- NULL
  ell$return.weights <- FALSE  # To save memory
  ell$renormalise <- FALSE
  ell$deriv <- NULL
  f <- switch(type,  # Lightweight, no derivatives
              EL0 = function(m) do.call(EL0, c(mu = m, deriv = TRUE, ell))$logelr,
              EL1 = function(m) do.call(EL1, c(mu = m, deriv = TRUE, ell))$logelr
  )
  ffp <- switch(type,
                EL0 = function(m) unlist(do.call(EL0, c(mu = m, deriv = TRUE, ell))[c("logelr", "deriv")]),
                EL1 = function(m) unlist(do.call(EL1, c(mu = m, deriv = TRUE, ell))[c("logelr", "deriv")])
                )
  n <- length(z)
  nc <- sum(ct)
  wm <- stats::weighted.mean(z, ct)  # Used in determining which end to process
  wv <- stats::weighted.mean((z-wm)^2, ct)
  a <- nc / wv  # Coefficient on the Wald parabola

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

  W  <- function(x) -0.5 * a*(x-wm)^2  # Wald parabola that should match log-ELR
  Wp <- function(x) -a*(x-wm)

  # xseq <- seq(min(z) , max(z), length.out = 101)
  # plot(xseq, sapply(xseq, f))
  # lines(xseq, sapply(xseq, W), col = 2)
  # plot(xseq, sapply(xseq, f) - sapply(xseq, W))

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
    xmax <- if (any(mu >= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(wm, xmax0)),  extendInt = "upX", maxiter = 20)$root else Inf
    xmin <- if (any(mu <= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(xmin0, wm)), extendInt = "downX", maxiter = 20)$root else -Inf
  }

  # Condition: -2*logELR(x1) > Wald(xmaxx1) --> logELR(x1) < W(x1) (ELR must be stronger)
  if (f(xmax) >= W(xmax)) {
    # xmax must be increased
    fdif <- function(x) f(x) - W(x)
    int  <- c(0.9999*xmax + 0.0001*max(z), 0.9*xmax + 0.1*max(z))
    xmax <- brentZero(fdif, interval = int, extendInt = "downX", tol = 1e-8, maxiter = 20)$root
  }
  if (f(xmin) <= W(xmin)) {
    # xmin must be decreased, i.e. shifted to the left to allow ELR to drop
    fdif <- function(x) f(x) - W(x)
    int  <- c(0.9*xmin + 0.1*min(z), 0.9999*xmin + 0.0001*min(z))
    xmin <- brentZero(fdif, interval = int, extendInt = "upX", tol = 1e-8, maxiter = 20)$root
    # Moving slightly closer to the outer wilds
    xmin <- 0.99*xmin + 0.01*min(z)
  }

  i.left  <- mu < xmin
  i.right <- mu > xmax
  i.mid   <- !(i.left | i.right)

  G <- function(t, x1) {  # Declaration before possible invocation
    x2 <- x1 + t
    ffd <- ffp(x1)  # Returns f and f'
    f1 <- ffd[1]  # f(x1)
    fp1 <- ffd[2] # f'(x1)
    W2  <- W(x2)  # W(x2)
    Wp2 <- Wp(x2) # W'(x2)
    et <- exp(t)
    a2 <-  (Wp2 - fp1) / (et - 1)
    fp1*t + a2 * (et-t-1) - (W2 - f1)
  }

  # Robust span for initial brackets
  tspan <- max(stats::IQR(z)/(qnorm(0.75)*2), stats::sd(z))

  fexleft <- fmid <- fexright <- aleft <- aright <- NULL

  if (any(i.mid)) {
    fmid <- sapply(mu[i.mid], f)
  }

  if (any(i.right)) {
    ffpright <- ffp(xmax)
    fright   <- ffpright[1]      # f(x1)
    fpright  <- ffpright[2]      # f'(x1)
    tr       <- wm - xmax - fpright/a
    tmax     <- min(tr - 1e-6, tspan)
    tpos     <- brentZero(function(t) G(t, x1 = xmax), c(1e-6, tmax), extendInt = "upX")$root

    x2  <- xmax + tpos
    et  <- exp(tpos)
    a2  <- (Wp(x2) - fpright) / (et - 1)    # <-- correct
    a1  <- fpright - a2
    a0  <- fright - a1*xmax - a2
    hright  <- function(x) a0 + a1*x + a2*exp(x - xmax)

    i.buf  <- i.right & (mu <= x2)
    i.wald <- mu > x2
    fexright <- numeric(length(mu))
    if (any(i.buf))  fexright[i.buf]  <- hright(mu[i.buf])
    if (any(i.wald)) fexright[i.wald] <- W(mu[i.wald])
    fexright <- fexright[i.right]
    aright   <- c(a0 = a0, a1 = a1, a2 = a2, a3 = xmax, x1 = xmax, x2 = x2)
  }

  if (any(i.left)) {
    ffpleft <- ffp(xmin)
    fleft   <- ffpleft[1]
    fpleft  <- ffpleft[2]
    tl      <- wm - xmin - fpleft/a
    tmin    <- max(-tspan, tl + 1e-6)
    tneg    <- brentZero(function(t) G(t, x1 = xmin), c(tmin, -1e-6), extendInt = "upX")$root

    x2  <- xmin + tneg
    et  <- exp(tneg)
    a2  <- (Wp(x2) - fpleft) / (et - 1)     # <-- correct
    a1  <- fpleft - a2
    a0  <- fleft - a1*xmin - a2
    hleft  <- function(x) a0 + a1*x + a2*exp(x - xmin)

    i.buf  <- i.left & (mu >= x2)
    i.wald <- mu < x2
    fexleft <- numeric(length(mu))
    if (any(i.buf))  fexleft[i.buf]  <- hleft(mu[i.buf])
    if (any(i.wald)) fexleft[i.wald] <- W(mu[i.wald])
    fexleft <- fexleft[i.left]
    aleft   <- c(a0 = a0, a1 = a1, a2 = a2, a3 = xmin, x1 = xmin, x2 = x2)
  }

  # curve(G(x, x1 = xmax), 0, tmax+5)
  # curve(G(x, x1 = xmin), tmax - tspan, tmax)

  logelr <- c(fexleft, fmid, fexright)

  attr(logelr, "xlim") <- c(xmin, xmax)
  attr(logelr, "parabola.coef") <- -0.5*a
  attr(logelr, "parabola.centre") <- wm
  attr(logelr, "bridge.left") <- aleft
  attr(logelr, "bridge.right") <- aright
  logelr
}

