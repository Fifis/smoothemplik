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
#' # Comparing ExEL2 vs ExEL1 with bridges containing exp(x)
#' z <- -4:4
#' ct <- 9:1
#' xseq <- seq(-7, 10.5, 0.1)
#' xl <- range(xseq)
#' a0 <- ExEL1(z, mu = xseq, ct = ct, exel.control = list(xlim = c(-11, 11)))
#' a1 <- ExEL1(z, mu = xseq, ct = ct)
#' a2 <- ExEL2(z, mu = xseq, ct = ct)
#' v1 <- attr(a1, "xlim")
#' v2 <- c(attr(a2, "bridge.left")[c("x1", "x2")], attr(a2, "bridge.right")[c("x1", "x2")])
#'
#' plot(xseq, a0, ylim = c(-300, 0), xlim = xl, main = "ExEL splices",
#'   bty = "n", xlab = "mu", ylab = "logELR(mu)")
#' lines(xseq, a1, col = 2, lwd = 2)
#' lines(xseq, a2, col = 4, lwd = 2)
#' abline(v = v2, lty = 3)
#' lines(xseq, attr(a2, "parabola.coef") * (xseq - attr(a2, "parabola.centre"))^2, lty = 2)
#' legend("topright", c("Taylor", "Wald", "ax^2"),
#'        col = c(2, 4, 1), lwd = c(2, 2, 1), lty = c(1, 1, 2))
#'
#' dx <- diff(xseq[1:2])
#' plot(xseq[-1], diff(a1)/dx, col = 2, type = "l", lwd = 2,
#'   main = "Derivatives of ExEL splice", bty = "n", ylim = c(-100, 100),
#'   xlab = "mu", ylab = "d/dmu logELR(mu)")
#' lines(xseq[-1], diff(a2)/dx, col = 4, lwd = 2)
#' abline(v = c(v1, v2), lty = 3, col = "#00000055")
#' legend("topright", c("Taylor", "Wald"), col = c(2, 4), lwd = 2)
#'
#' # Multivariate extension
#' set.seed(1)
#' X <- cbind(rchisq(30, 3), rchisq(30, 3))
#' ct <- runif(30)
#' -2*ExEL1(X, mu = c(-1, -1),  ct = ct)  # Outside the hull
#' -2*ExEL2(X, mu = c(-1, -1),  ct = ct)
ExEL1 <- function(z, mu, type = c("auto", "EL0", "EL1"),
                  exel.control = list(xlim = "auto", fmax = NA, p = 0.999, df = NA), ...) {
  type <- match.arg(type)
  if (type == "auto") type <- if (NCOL(z) > 1) "EL1" else "EL0"

  d <- NCOL(z)
  if (d > 1) type <- "EL1"  # Against hard-coding EL0

  ell <- list(z = z, ...)
  ct <- ell$ct
  n <- NROW(z)
  if (is.null(ct)) ct <- ell$ct <- rep(1, n)
  z <- if (d == 1) z[ct != 0] else z[ct != 0, , drop = FALSE]
  ct <- ct[ct != 0]
  # Cleaning the dots from undesirable user arguments
  ell$mu <- NULL
  ell$return.weights <- FALSE  # To save memory
  ell$renormalise <- FALSE
  ell$deriv <- NULL

  def.ctl <- list(xlim = "auto", fmax = NA_real_, p = 0.999, df = NA_integer_)
  exel.control <- if (is.null(exel.control)) def.ctl else utils::modifyList(def.ctl, exel.control)
  fmax   <- exel.control$fmax
  p  <- if (is.na(exel.control$p))  0.999 else exel.control$p
  df <- if (is.na(exel.control$df)) d else exel.control$df
  if (!is.finite(fmax)) fmax <- stats::qchisq(p = p, df = df)  # exel.control may still hold the NA

  fl <- switch(type, # f value + derivatives in a list; ct is already inside ell
               EL0 = function(m) do.call(EL0, c(list(mu = m, deriv = TRUE), ell))[c("logelr", "deriv")],
               EL1 = function(m) do.call(EL1, c(list(mu = m, deriv = TRUE), ell))[c("logelr", "deriv")])
  f <- switch(type,
              EL0 = function(m) do.call(EL0, c(list(mu = m, deriv = FALSE), ell))$logelr,
              EL1 = function(m) do.call(EL1, c(list(mu = m, deriv = FALSE), ell))$logelr)

  wm <- if (d == 1) stats::weighted.mean(z, ct) else drop(crossprod(ct, z)) / sum(ct)

  if (d > 1) {
    if (!identical(exel.control$xlim, "auto")) stop("ExEL1 (multivariate): only xlim='auto' is supported.")
    # Allow mu as a vector (length d) or a matrix (rows = targets)
    mumat <- if (is.null(dim(mu))) matrix(mu, ncol = d) else {
      if (ncol(mu) != d) stop("mu must have d columns.")
      mu
    }

    # Scalar objective along a ray (for Brent): g(t) = -2 f_v(t) - fmax

    # Evaluate each target
    logelr <- apply(mumat, 1, function(murow) {
      vv <- murow - wm
      t <- sqrt(sum(vv^2))
      if (t == 0) return(f(murow))  # At centre: return true EL
      v <- vv / t  # Direction between wm and mu

      # Find tmax with Brent
      tzero <- function(t) -2*f(wm + t*v) - fmax
      tmax <- tryCatch(brentZero(tzero, c(0, max(t, 1e-6)), extendInt = "right", maxiter = 100)$root, error = function(e) NA_real_)
      # tseq <- seq(0, t, length.out = 51)
      # plot(tseq, sapply(tseq, tzero)); abline(h = 0, lty = 2)

      if (!is.finite(tmax)) {
        warning("ExEL1: could not find the desired cut-off point. Trying shrinkage to find at least a rough finite value.")
        tmax <- t
        for (i in 1:20) {
          tmax <- tmax * 0.5
          tzero.safe <- tryCatch(tzero(tmax), error = function(e) return(NA_real_))
          if (is.finite(tzero.safe)) {
            tmax <- tryCatch(brentZero(tzero, c(0, tmax), extendInt = "right", maxiter = 100)$root, error = function(e) NA_real_)
            if (is.finite(tmax)) break  # Else shrink again and retry
          }
        }
        if (i == 20) stop("ExEL1: could not find the desired cut-off point. Report this bug to GitHub.")
      }

      if (t <= tmax) return(f(murow))  # No extrapolation, true EL

      # Analytic directional derivatives
      mumax <- wm + tmax * v
      flmax  <- fl(mumax)

      # Taylor branch along the ray
      a <- getParabola(tmax, flmax$logelr, flmax$deriv[1], flmax$deriv[2])
      return(a[1]*t^2 + a[2]*t + a[3])
    })

    return(as.numeric(logelr))
  } else {
    # Evaluation range to determine which mu are acceptable
    # A value must be computed for all cases
    zu <- sort(unique(z))
    nu <- length(zu)
    dzu <- diff(zu)
    if (nu >= 10) {
      rgap <- stats::median(dzu)
    } else if (nu >= 5) {
      rgap <- stats::median(dzu) * ((nu-4.5)/5)  # Starting at 0.1 median, ending at 0.9 median
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
      # Fast and slightly inaccurate defaults for a quick search
      xmax <- if (any(mu >= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(wm, xmax0)),  extendInt = "upX", maxiter = 50)$root else Inf
      xmin <- if (any(mu <= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(xmin0, wm)), extendInt = "downX", maxiter = 50)$root else -Inf
    }

    i.left  <- mu < xmin
    i.right <- mu > xmax
    i.mid   <- !(i.left | i.right)

    # Analytical derivatives require lambda as well
    fleft <- fmid <- fright <- aleft <- aright <- NULL
    if (any(i.mid)) {
      fmid <- sapply(mu[i.mid], f)
    }
    if (any(i.left)) {
      fminlist <- fl(xmin)
      fp  <- fminlist$deriv[1]
      fpp <- fminlist$deriv[2]
      if (type == "EL1") fp <- sign(xmin - wm) * fp  # Convert directional derivative (along ray) to d/dx if it came from EL1
      aleft <- getParabola(xmin, fminlist$logelr, fp, fpp)
      fleft <- aleft[1]*mu[i.left]^2 + aleft[2]*mu[i.left] + aleft[3]
    }
    if (any(i.right)) {
      fmaxlist <- fl(xmax)
      fp  <- fmaxlist$deriv[1]
      fpp <- fmaxlist$deriv[2]
      if (type == "EL1") fp <- sign(xmax - wm) * fp
      aright <- getParabola(xmax, fmaxlist$logelr, fmaxlist$deriv[1], fmaxlist$deriv[2])
      fright <- aright[1]*mu[i.right]^2 + aright[2]*mu[i.right] + aright[3]
    }
    logelr <- c(fleft, fmid, fright)

    attr(logelr, "xlim") <- c(xmin, xmax)
    attr(logelr, "parabola.left") <- aleft
    attr(logelr, "parabola.right") <- aright
    return(logelr)
  }
}


#' @rdname ExEL1
#' @export
ExEL2 <- function(z, mu, type = c("auto", "EL0", "EL1"),
                  exel.control = list(xlim = "auto", fmax = NA, p = 0.999, df = NA), ...) {
  type <- match.arg(type)
  if (type == "auto") type <- if (NCOL(z) > 1) "EL1" else "EL0"

  d <- NCOL(z)
  if (d > 1) type <- "EL1"  # Against hard-coding EL0

  ell <- list(z = z, ...)
  ct <- ell$ct
  if (is.null(ct)) ct <- ell$ct <- rep(1, NROW(z))
  z <- if (d == 1) z[ct != 0] else z[ct != 0, , drop = FALSE]
  ct <- ct[ct != 0]
  ell$mu <- NULL
  ell$return.weights <- FALSE  # To save memory
  ell$renormalise <- FALSE
  ell$deriv <- NULL

  f <- switch(type,  # Lightweight, no derivatives
              EL0 = function(m) do.call(EL0, c(list(mu = m, deriv = FALSE), ell))$logelr,
              EL1 = function(m) do.call(EL1, c(list(mu = m, deriv = FALSE), ell))$logelr
  )
  fl <- switch(type,
                EL0 = function(m) do.call(EL0, c(list(mu = m, deriv = TRUE), ell))[c("logelr", "deriv")],
                EL1 = function(m) do.call(EL1, c(list(mu = m, deriv = TRUE), ell))[c("logelr", "deriv")]
                )

  def.ctl <- list(xlim = "auto", fmax = NA_real_, p = 0.999, df = NA_integer_)
  exel.control <- if (is.null(exel.control)) def.ctl else utils::modifyList(def.ctl, exel.control)
  fmax   <- exel.control$fmax
  p  <- if (is.na(exel.control$p))  0.999 else exel.control$p
  df <- if (is.na(exel.control$df)) d else exel.control$df
  if (!is.finite(fmax)) fmax <- stats::qchisq(p = p, df = df)  # exel.control may still hold the NA

  if (d > 1) {
    if (!identical(exel.control$xlim, "auto")) stop("ExEL2 (multivariate): only xlim='auto' is supported.")

    # Weighted mean and weighted covariance for Wald tails
    nc <- sum(ct)
    wm <- as.numeric(colSums(ct * z) / nc)
    Zc <- sweep(z, 2, wm, FUN = "-")
    Sig <- (t(Zc) %*% (ct * Zc)) / nc  # Weighted Var(Z)
    # Small ridge for numerical safety (does not matter on log scale if tiny)
    if (!isTRUE(all(is.finite(Sig)))) stop("Non-finite covariance.")
    eig_min <- min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values)
    if (eig_min <= 0) {
      Sig <- Sig + diag(max(1e-12, 1e-8 - eig_min), nrow(Sig))
    }

    Wv  <- function(t, a) -0.5 * a * t^2
    Wpv <- function(t, a) -a * t

    # Accept mu as vector or matrix (rows = hypotheses)
    mumat <- if (is.null(dim(mu))) matrix(mu, ncol = d) else {
      if (ncol(mu) != d) stop("mu must have d columns.")
      mu
    }
    nQ <- nrow(mumat)

    # Numerical stability for a single ray
    # Helper: build a single splice for a given ray v
    build_splice <- function(v) {
      v <- as.numeric(v / sqrt(sum(v^2)))
      a <- nc / drop(t(v) %*% Sig %*% v)

      # Robust scale & reach alona1g the ray (from data)
      proj   <- drop(Zc %*% v)                 # weighted projections
      sd_v   <- sqrt(sum(ct * proj^2) / nc)    # SD of projection (weights ct)
      iqr_v  <- stats::mad(proj)
      t_hi0  <- max(0.01*sd_v, min(sd_v, iqr_v))

      # tseq <- seq(0, 1, length.out = 101)
      # fseq <- sapply(tseq, function(t) f(wm + t * v))
      # plot(tseq, fseq)
      # lines(tseq, -0.5*a*tseq^2, col = 2)
      # abline(h = -0.5*qchisq(0.999, 2), lty = 2)
      #
      # linfun <- function(x) f1 + fp1*(x-t1)
      # lines(tseq, linfun(tseq), col = 3)

      # (1) cut radius t1 : solve -2 f(wm + t v) = fmax
      g  <- function(t) -2 * f(wm + t * v) - fmax
      # curve(g(x), 0, t_hi0)
      # Add trace=2 here because this is one of the most error-prone parts
      t1 <- tryCatch(brentZero(g, c(0, t_hi0), extendInt = "right", maxiter = 100)$root, error = function(e) NA_real_)
      if (!is.finite(t1)) stop("ExEL2 splice: could not locate cut t1 (root of -2*f = fmax). Try a different fmax or check EL convergence.")

      f1 <- f(wm + t1 * v)
      fp1 <- fl(wm + t1 * v)$deriv[1]  # Directional derivative f'_v at the cut (ELCPP returns it already)
      W1 <- Wv(t1, a)

      left_branch <- (fp1 > 0)  # The dome of log-ELR is growing
      # Choice of bridge (see univariate logic):
      #   use +exp when (left & f1>W1) OR (right & f1<W1); else use -exp
      use_pos_exp <- (left_branch && (f1 > W1)) || (!left_branch && (f1 < W1))

      # Bridge length tau : solve G(tau) = 0 with tau > 0
      if (use_pos_exp) {
        # h(t) = a0 + a1*t + a2*exp(t - t1)
        G <- function(tau) {
          t2 <- t1 + tau
          et <- exp(tau)
          a2 <- (Wpv(t2, a) - fp1) / (et - 1)
          fp1 * tau + a2 * (et - tau - 1) - (Wv(t2, a) - f1)
        }
      } else {
        # h(t) = a0 + a1*t + a2*exp(-(t - t1))
        G <- function(tau) {
          t2 <- t1 + tau
          e  <- exp(-tau)
          d  <- 1 - e
          W2 <- Wv(t2, a)
          Wp2 <- Wpv(t2, a)
          a1 <- (Wp2 - fp1 * e) / d
          a1 * (tau - d) + fp1 * d - (W2 - f1)
        }
      }

      # bracket tau robustly to the right
      # tseq <- seq(0, t_hi0, length.out = 51)
      # gseq <- sapply(tseq, G)
      # plot(tseq, gseq)
      tau <- tryCatch(suppressWarnings(brentZero(G, c(1e-6, max(1e-2, 0.1*t_hi0)), extendInt = "right", maxiter = 100)$root), error = function(e) NA_real_)
      if (!is.finite(tau) || tau <= 1.01e-6) {
        for (i in 1:30) {
          # warning("ExEL2 splice: root for tau failed after bracketing. Starting extrapolating earlier by 80%.")
          t1 <- t1 * 0.8
          f1 <- f(wm + t1 * v)
          fp1 <- fl(wm + t1 * v)$deriv[1]  # Directional derivative f'_v at the cut (ELCPP returns it already)
          W1 <- Wv(t1, a)
          tau <- tryCatch(suppressWarnings(brentZero(G, c(1e-6, max(1e-2, 0.1*t_hi0)), extendInt = "right", maxiter = 100)$root), error = function(e) NA_real_)
          if (is.finite(tau) && tau >= 1.01e-6) break
        }
        if (i == 30) {
          warning("Could not find an appropriate configuration -- returning an immediate Wald transition.")
          tau <- 1e-8
        }
      }

      t2 <- t1 + tau

      if (use_pos_exp) {
        et <- exp(tau)
        a2 <- (Wpv(t2, a) - fp1) / (et - 1)
        a1 <- fp1 - a2
        a0 <- f1 - a1 * t1 - a2
      } else {
        e  <- exp(-tau)
        d  <- 1 - e
        a1 <- (Wpv(t2, a) - fp1 * e) / d
        a2 <- a1 - fp1
        a0 <- f1 - a1 * t1 - a2
      }
      # cat(sprintf("ray: a=%.6g, t1=%.6g, t2=%.6g, f1=%.6g, fp1=%.6g\n", a, t1, t2, f1, fp1))
      list(a = a, v = v, t1 = t1, t2 = t2, f1 = f1, fp1 = fp1, a0 = a0, a1 = a1, a2 = a2, bridge = if (use_pos_exp) "pos" else "neg")
    }

    # Are all queries on (approximately) the same ray?
    vv   <- sweep(mumat, 2, wm, FUN = "-")
    tlen <- sqrt(rowSums(vv^2))
    same_ray <- FALSE
    if (nQ >= 2 && all(is.finite(tlen))) {
      idx   <- which.max(tlen)
      ref_v <- vv[idx, ]
      if (tlen[idx] > 0) {
        ref_v <- ref_v / tlen[idx]              # unit reference
        dir   <- vv
        if (any(tlen > 0)) dir[tlen > 0, ] <- sweep(dir[tlen > 0, , drop = FALSE], 1, tlen[tlen > 0], "/")
        cosang <- as.numeric(dir[tlen > 0, , drop = FALSE] %*% ref_v)
        cosang <- pmin(1, pmax(-1, cosang))
        same_ray <- all(cosang > 0) && max(acos(cosang)) < 1e-6
      }
    }

    if (same_ray) {
      v  <- as.numeric(vv[which.max(tlen), ] / tlen[which.max(tlen)])
      sp <- build_splice(v)

      eval_on_ray <- function(t) {
        if (t <= sp$t1) return(f(wm + t * sp$v))
        if (t <= sp$t2) {
          if (sp$bridge == "pos") {
            return(sp$a0 + sp$a1 * t + sp$a2 * exp(t - sp$t1))
          } else if (sp$bridge == "neg") {
            return(sp$a0 + sp$a1 * t + sp$a2 * exp(-(t - sp$t1)))
          }
        }
        return(Wv(t, sp$a))
      }

      out <- vapply(tlen, eval_on_ray, numeric(1))
      # cat(sprintf("same-ray: t in [%.6g, %.6g]; t1=%.6g, t2=%.6g\n", min(tlen), max(tlen), sp$t1, sp$t2))
      return(as.numeric(out))
    }

    # Fallback: different directions — solve independently (previous behaviour)
    logelr <- apply(mumat, 1, function(murow) {
      vv <- murow - wm
      t  <- sqrt(sum(vv^2))
      if (t == 0) return(f(murow))
      v  <- vv / t
      sp <- build_splice(v)
      if (t <= sp$t1) return(f(murow))
      if (t <= sp$t2) {
        if (sp$bridge == "pos") return(sp$a0 + sp$a1 * t + sp$a2 * exp(t - sp$t1))
        return(sp$a0 + sp$a1 * t + sp$a2 * exp(-(t - sp$t1)))
      }
      Wv(t, sp$a)
    })
    return(as.numeric(logelr))
  } else {   ################################################## Univariate
    n <- NROW(z)
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
      rgap <- stats::median(dzu) * ((nu-4.5)/5)  # Starting at 0.1 median, ending at 0.9 median
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
    # abline(v = c(xmin0, xmax0), lty = 2)
    # plot(xseq, sapply(xseq, f) - sapply(xseq, W))

    if (identical(exel.control$xlim, "auto") && is.na(exel.control$fmax)) {
      # No user xlim supplied, no user fmax supplied
      xmin <- xmin0
      xmax <- xmax0
    } else if (is.numeric(exel.control$xlim)) {
      xmin <- min(exel.control$xlim)
      xmax <- max(exel.control$xlim)
    } else {  # xmax and xmin must be searched for
      # Fast and slightly inaccurate defaults for a quick search
      xmax <- if (any(mu >= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(wm, xmax0)),  extendInt = "upX", maxiter = 20)$root else Inf
      xmin <- if (any(mu <= wm)) brentZero(function(x) -2*f(x) - fmax, sort(c(xmin0, wm)), extendInt = "downX", maxiter = 20)$root else -Inf
    }

    # Transition to a lower function
    Gposexp <- function(x1, t) {
      x2 <- x1 + t
      ffd <- fl(x1)   # List with $logelr and $deriv
      f1  <- ffd$logelr
      fp1 <- ffd$deriv[1]
      W2  <- W(x2)  # W(x2)
      Wp2 <- Wp(x2) # W'(x2)
      et <- exp(t)
      a2 <-  (Wp2 - fp1) / (et - 1)
      fp1*t + a2 * (et-t-1) - (W2 - f1)
    }
    # Transition to a higher function
    Gnegexp <- function(x1, t) {
      ffd <- fl(x1)
      f1 <- ffd$logelr
      fp1 <- ffd$deriv[1]
      x2  <- x1 + t
      W2  <- W(x2);  Wp2 <- Wp(x2)
      e   <- exp(-t)
      d   <- 1 - e
      a1  <- (Wp2 - fp1*e) / d
      a1*(t - d) + fp1*d - (W2 - f1)   # = 0 at the solution
    }
    # Bridges
    hposexp <- function(x, x1, a) a[1] + a[2]*x + a[3]*exp(x - x1)
    hnegexp <- function(x, x1, a) a[1] + a[2]*x + a[3]*exp(-(x - x1))

    i.left  <- mu < xmin
    i.right <- mu > xmax
    i.mid   <- !(i.left | i.right)

    # Robust span for initial brackets
    # TODO: revise
    tspan <- max(stats::IQR(z)/(stats::qnorm(0.75)*2), stats::sd(z))

    fexleft <- fmid <- fexright <- aleft <- aright <- NULL

    if (any(i.mid)) {
      fmid <- sapply(mu[i.mid], f)
    }

    if (any(i.right)) {
      x1 <- xmax
      for (i in 1:20) {  # Fail-safe
        ffpright <- fl(x1)
        f1   <- ffpright$logelr      # f(x1)
        fp1  <- ffpright$deriv[1]      # f'(x1)
        if (type == "EL1") fp1 <- sign(x1 - wm) * fp1
        if (W(x1) > f1) {  # Type-1 bridge with exp(x)
          tr       <- wm - x1 - fp1/a
          tmax     <- min(tr - 1e-6, tspan)
          H <- function(t) Gposexp(x1 = x1, t = t)
          # xseq <- seq(0, 20, length.out = 101)
          # gseq <- sapply(xseq, H)
          # plot(xseq, gseq)
          troot <- tryCatch(suppressWarnings(brentZero(H, c(1e-6, tmax), extendInt = "right", maxiter = 100)), error = function(e) NULL)
          if (is.null(troot) || troot$iter >= 100) {
            a2 <- Inf  # For failure
            if (is.null(troot)) tspan <- tspan * 0.8
          } else {
            tpos <- troot$root
            x2  <- x1 + tpos
            et  <- exp(tpos)
            a2  <- (Wp(x2) - fp1) / (et - 1)
            a1  <- fp1 - a2
            a0  <- f1 - a1*x1 - a2
            h  <- function(x) hposexp(x, x1 = x1, a = c(a0, a1, a2))
          }
        } else {  # Type-2 bridge with exp(-x)
          # Wald below EL (transition to "higher"); search x2 > x1 for the root of Gnegexp
          H <- function(t) Gnegexp(x1 = x1, t = t)
          # xseq <- seq(x1, x1+2, length.out = 101)
          # plot(xseq, H(xseq))
          tmax <- max(1e-2, 0.1*tspan)
          troot <- tryCatch(suppressWarnings(brentZero(H, c(1e-6, tmax), extendInt = "right", maxiter = 100)), error = function(e) NULL)
          if (is.null(troot) || troot$iter >= 100) {
            a2 <- Inf
            if (is.null(troot)) tspan <- tspan * 0.8
          } else {
            x2   <- x1 + troot$root
            emt  <- exp(-troot$root)
            a2 <- (Wp(x2) - fp1) / (1 - emt)
            a1 <- fp1 + a2
            a0 <- f1 - a1*x1 - a2
            h  <- function(x) hnegexp(x, x1 = x1, a = c(a0, a1, a2))
          }
        }
        # xseq <- seq(min(0, min(z)), max(0, 0.2, max(z)), length.out = 301)
        # plot(xseq, sapply(xseq, f), ylim = c(-120, 0))
        # lines(xseq, sapply(xseq, W), col = 2)
        # lines(xseq, h(xseq), col = 4)
        # abline(v = c(x1, x2), lty = 2)
        if (a2 < 0) {
          break
        } else {
          # message("Wrong exponent coefficient sign, bridge lost concavity; attempt ", i, " to shift x1 20% closer to the mean.")
          x1 <- 0.8*x1 + 0.2*wm
          x2 <- x1
          a0 <- a1 <- a2 <- NA
        }
      }
      failure <- i >= 20

      i.buf  <- i.right & (mu <= x2)
      i.wald <- i.right & (mu > x2)
      fexright <- numeric(length(mu))
      # if (failure) warning("Could not find a suitable bridge (possible with severe outliers). Returning the usual Wald.")
      if (any(i.buf))  fexright[i.buf]  <- h(mu[i.buf])
      if (any(i.wald)) fexright[i.wald] <- W(mu[i.wald])
      fexright <- fexright[i.right]
      aright   <- c(a0 = a0, a1 = a1, a2 = a2, x1 = x1, x2 = x2)
    }

    if (any(i.left)) {
      x1 <- xmin
      for (i in 1:20) {  # Fail-safe
        ffpleft <- fl(x1)
        f1   <- ffpleft$logelr
        fp1  <- ffpleft$deriv[1]
        if (type == "EL1") fp1 <- sign(x1 - wm) * fp1

        if (W(x1) > f1) {
          # Type-2 bridge with exp(-x)
          # Wald below EL (transition to higher) -- search x2 < x1 for the root of Gnegexp
          H <- function(t) Gnegexp(x1 = x1, t = t)
          # xseq <- seq(-10, 0, length.out = 101)
          # plot(xseq, H(xseq))
          tmin <- min(-1e-6, -0.1*tspan)
          troot <- tryCatch(suppressWarnings(brentZero(H, c(tmin, -1e-6), extendInt = "left", maxiter = 100)), error = function(e) NULL)
          if (is.null(troot) || troot$iter >= 100) {
            a2 <- Inf
            x2 <- x1
            if (is.null(troot)) tspan <- tspan * 0.8  # The only possible error is 'large span, infinite H'
          } else {
            tneg <- troot$root
            x2   <- x1 + tneg
            emt  <- exp(tneg)
            a1   <- (Wp(x2)*emt - fp1) / (emt - 1)
            a2 <- a1 - fp1
            a0 <- f1 - a1 * x1 - a2
            h  <- function(x) hnegexp(x, x1 = x1, a = c(a0, a1, a2))
          }
        } else {  # Type-1 bridge with a0 + a1*x + a2*exp(x-x1)
          tl      <- wm - x1 - fp1/a
          tmin    <- max(-tspan, tl + 1e-6)
          H <- function(t) Gposexp(x1 = x1, t = t)
          # xseq <- seq(-10, 0, length.out = 101)
          # hseq <- sapply(xseq, H)
          # plot(xseq, hseq)
          troot <- tryCatch(suppressWarnings(brentZero(H, c(tmin, -1e-6), extendInt = "left", maxiter = 100)), error = function(e) NULL)
          if (is.null(troot) || troot$iter >= 100) {
            a2 <- Inf
            if (is.null(troot)) tspan <- tspan * 0.8
          } else {
            tneg    <- troot$root
            x2  <- x1 + tneg
            et  <- exp(tneg)
            a2  <- (Wp(x2) - fp1) / (et - 1)
            a1  <- fp1 - a2
            a0  <- f1 - a1*x1 - a2
            h  <- function(x) hposexp(x, x1 = x1, a = c(a0, a1, a2))
          }
        }
        # xseq <- seq(min(0, min(z)), max(0, 0.2, max(z)), length.out = 301)
        # plot(xseq, sapply(xseq, f), ylim = c(-20, 0))
        # lines(xseq, sapply(xseq, W), col = 2)
        # lines(xseq, h(xseq), col = 4)
        # abline(v = c(x1, x2), lty = 2)
        if (a2 < 0) {
          break
        } else {
          # message("Wrong exponent coefficient sign, bridge lost concavity; attempt ", i, " to shift x1 20% closer to the mean.")
          x1 <- 0.8*x1 + 0.2*wm
          x2 <- x1
          a0 <- a1 <- a2 <- NA
        }
      }
      failure <- i >= 20

      i.buf  <- i.left & (mu >= x2)
      i.wald <- i.left & (mu < x2)
      fexleft <- numeric(length(mu))
      # if (failure) warning("Could not find a suitable left bridge (possible with severe outliers). Returning the usual Wald.")
      if (any(i.buf))  fexleft[i.buf]  <- h(mu[i.buf])
      if (any(i.wald)) fexleft[i.wald] <- W(mu[i.wald])
      fexleft <- fexleft[i.left]
      aleft   <- c(a0 = a0, a1 = a1, a2 = a2, x1 = x1, x2 = x2)
    }

    logelr <- c(fexleft, fmid, fexright)

    attr(logelr, "xlim") <- c(xmin, xmax)
    attr(logelr, "parabola.coef") <- -0.5*a
    attr(logelr, "parabola.centre") <- wm
    attr(logelr, "bridge.left") <- aleft
    attr(logelr, "bridge.right") <- aright
    logelr
  }

}

