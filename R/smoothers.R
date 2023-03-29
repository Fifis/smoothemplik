#' Probability integral transform
#'
#' @param x A numeric vector of data points.
#' @param xgrid A numeric vector. If supplied, then the transformed function at the grid points different from \code{x} takes values
#' equidistant between themselves and the ends of the interval to which they belong.
#'
#' @return A numeric vector of values strictly between 0 and 1 of the same length as \code{xgrid} (or \code{x}, if \code{xgrid} is \code{NULL}).
#' @export
#'
#' @examples
#' set.seed(2)
#' x1 <- c(4, 3, 7, 10, 2, 2, 7, 2, 5, 6)
#' x2 <- sample(c(0, 0.5, 1, 2, 2.5, 3, 3.5, 10, 100), 25, replace = TRUE)
#' l <- length(x1)
#' pit(x1)
#'
#' plot(pit(x1), ecdf(x1)(x1), xlim = c(0, 1), ylim = c(0, 1), asp = 1)
#' abline(v = seq(0.5 / l, 1 - 0.5 / l, length.out = l), col = "#00000044", lty = 2)
#' abline(v = c(0, 1))
#' points(pit(x1, x2), ecdf(x1)(x2), pch = 16, col = "#CC000088", cex = 0.9)
#' abline(v = pit(x1, x2), col = "#CC000044", lty = 2)
#'
#' x1 <- c(1, 1, 3, 4, 6)
#' x2 <- c(0, 2, 2, 5.9, 7, 8)
#' pit(x1)
#' pit(x1, x2)
#'
#' set.seed(1)
#' l <- 10
#' x1 <- rlnorm(l)
#' x2 <- sample(c(x1, rlnorm(10)))
#' plot(pit(x1), ecdf(x1)(x1), xlim = c(0, 1), ylim = c(0, 1), asp = 1)
#' abline(v = seq(0.5 / l, 1 - 0.5 / l, length.out = l), col = "#00000044", lty = 2)
#' abline(v = c(0, 1))
#' points(pit(x1, x2), ecdf(x1)(x2), pch = 16, col = "#CC000088", cex = 0.9)
pit <- function(x, xgrid = NULL) {
  x.transformed <- stats::ecdf(x)(x) - 0.5 / length(x)
  if (is.null(xgrid)) return(x.transformed) else {
    if (isTRUE(all.equal(x, xgrid, tolerance = .Machine$double.eps))) return(x.transformed)
    x.uniq.sorted <- sort(unique(x))
    ecdf.uniq.sorted <- sort(x.transformed[!duplicated(x)])
    xgrid.uniq.sorted <- sort(unique(xgrid))
    xgrid.cut <- cut(xgrid, breaks = c(-Inf, x.uniq.sorted, Inf))
    exact <- xgrid %in% x.uniq.sorted
    xgrid.cut[xgrid %in% x.uniq.sorted] <- NA # Exact matches will be dealt with at the last step
    xgrid.list <- split(xgrid, xgrid.cut)
    ecdf.uniq01 <- c(0, ecdf.uniq.sorted, 1)
    xgrid.list.unique <- lapply(xgrid.list, function(x) sort(unique(x)))
    xgrid.spaced <- lapply(1:(length(xgrid.list.unique)), function(i) {
      xg <- xgrid.list.unique[[i]]
      l <- length(xg)
      if (l > 0) return(seq(ecdf.uniq01[i], ecdf.uniq01[i+1], length.out = l + 2)[c(-1, -l-2)]) else return(NULL)
    })
    xgrid.spaced <- unlist(xgrid.spaced)
    xgrid.uniq.sorted.noorig <- xgrid.uniq.sorted[!(xgrid.uniq.sorted %in% x)]
    ret <- xgrid.spaced[match(xgrid, xgrid.uniq.sorted.noorig)]
    ret[exact] <- x.transformed[match(xgrid[exact], x)]
    return(ret)
  }
}

kernelWeightsOneR <- function(x, # A numeric vector
                              xgrid, # A numeric vector
                              bw, # A scalar
                              gaussian = TRUE) {
  kernelFun <- function(x) if (gaussian) stats::dnorm(x) else 0.75 / sqrt(5) * (1 - x^2 / 5) * as.numeric(abs(x) <= sqrt(5))
  K <- kernelFun(outer(xgrid, x, "-") / bw)
  return(K)
}

kernelDensityOneR <- function(x, # A numeric vector
                              xgrid, # A numeric vector
                              bw, # A scalar
                              gaussian = TRUE) {
  nx <- length(x)
  K <- kernelWeightsOneR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
  dens <- rowSums(K) / (nx * bw)
  return(dens)
}

kernelSmoothOneR <- function(x, # A numeric vector
                             y, # A numeric vector of the same length
                             xgrid, # A numeric vector
                             bw, # A scalar
                             gaussian = TRUE,
                             LOO = FALSE) {
  nx <- length(x)
  K <- kernelWeightsOneR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
  if (LOO) diag(K) <- 0
  num <- rowSums(sweep(K, 2, y, "*"))
  den <- rowSums(K)
  muhat <- num / den
  return(muhat)
}

kernelWeightsMultiR <- function(x, # A matrix with d columns
                                xgrid, # A matrix with d columns
                                bw, # A vector of length d
                                gaussian = TRUE) {
  kernelFun <- function(x) if (gaussian) stats::dnorm(x) else 0.75 / sqrt(5) * (1 - x^2 / 5) * as.numeric(abs(x) <= sqrt(5))
  nx <- nrow(x)
  ng <- nrow(xgrid)
  d <- ncol(xgrid)
  pk <- matrix(1, nrow = ng, ncol = nx) # Accumulating the product kernel
  for (i in 1:d) {
    K <- kernelFun(outer(xgrid[, i], x[, i], "-") / bw[i])
    pk <- pk * K
  }
  return(pk)
}

kernelDensityMultiR <- function(x, # A matrix with d columns
                                xgrid, # A matrix with d columns
                                bw, # A vector of length d
                                gaussian = TRUE) {
  nx <- nrow(x)
  pk <- kernelWeightsMultiR(x=x, xgrid=xgrid, bw=bw, gaussian=gaussian)
  dens <- rowSums(pk) / (nx * prod(bw))
  return(dens)
}

kernelSmoothMultiR <- function(x,
                               y,
                               xgrid,
                               bw,
                               gaussian = TRUE,
                               LOO = FALSE) {
  nx <- nrow(x)
  pk <- kernelWeightsMultiR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
  if (LOO) diag(pk) <- 0
  num <- rowSums(sweep(pk, 2, y, "*"))
  den <- rowSums(pk)
  muhat <- num / den
  return(muhat)
}

#' Get kernel-based weights
#'
#' @param x A vector or a matrix of data points involved in weight computation.
#' @param xgrid A vector or a matrix of data points for which the weights of every observation in \code{x} are computed. If \code{NULL}, then
#' \code{x} itself is used as the grid.
#' @param bw Kernel bandwidth. Since it is the crucial parameter in many applications, throws a warning if not supplied, and then, Silverman's
#' rule of thumb (via \code{bw.row()}) is applied to every dimension of \code{x}.
#' @param gaussian If TRUE, the normal (Gaussian) product kernel with full support is used, otherwise Epanechnikov.
#' @param PIT If TRUE, the Probability Integral Transform (PIT) is applied to all columns of \code{x} via \code{ecdf} in order to map all values
#'   into the [0, 1] range. May be an integer vector of indices of columns to which the PIT should be applied.
#' @param method "R" or "CPP". If "CPP", then faster Rcpp functions are invoked, otherwise use a pure R implementation.
#'
#' Note that if \code{pit = TRUE}, then the kernel-based weights become nearest-neighbour weights (i.e. not much different from the ones used
#' internally in the built-in \code{loess} function) since the distances now depend on the ordering of data, not the values per se.
#'
#' @return A matrix
#' @export
#'
#' @examples
kernelWeights <- function(x,
                          xgrid = NULL,
                          bw = NULL,
                          gaussian = TRUE,
                          PIT = FALSE,
                          method = "CPP"
) {
  one.dimensional <- is.vector(x)
  if (is.null(xgrid)) xgrid <- x
  if (one.dimensional) {
    if (PIT) {
      x <- pit(x)
      xgrid <- pit(x = x, xgrid = xgrid)
    }
    if (is.null(bw)) {
      bw <- bw.rot(x)
      warning(paste0("No bandwidth supplied, using Silverman's one-dimensional rule of thumb: bw = ", round(bw, 5), "."))
    }
    if (length(bw) != 1) stop("The bandwidth for one-dimensional smoothing must be a scalar.")
    if (!is.vector(xgrid)) stop("The grid for one-dimensional smoothing must be a vector.")
    result <- switch(method,
                     R   = kernelWeightsOneR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
                     CPP = kernelWeightsOneCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
    )
  } else {
    d <- ncol(x)
    if (!is.matrix(xgrid)) stop("x and xgrid must be both matrices.")
    if (d != ncol(xgrid)) stop("x and xgrid must be have the same number of columns (i.e. the same dimension).")
    if (PIT) {
      for (i in 1:d) {
        x[, i] <- pit(x[, i])
        xgrid[, i] <- pit(x = x[, i], xgrid = xgrid[, i])
      }
    }
    if (is.null(bw)) {
      bw <- bw.rot(x)
      warning("No bandwidth supplied, using Silverman's multi-dimensional rule of thumb: bw = (", paste(round(bw, 5), collapse = ", "), ").")
    }
    if (length(bw) == 1) bw <- rep(bw, d)
    if (length(bw) != d) stop("The number of data dimensions is not the same as the number of bandwidths!")
    result <- switch(method,
                     R   = kernelWeightsMultiR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
                     CPP = kernelWeightsMultiCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
    )
  }
  return(result)
}


#' Kernel density estimation
#'
#' @param x A numeric vector or numeric matrix of predictors.
#' @param xgrid A numeric vector or numeric matrix with \code{ncol(xgrid) = ncol(x)} of points at which the density is estimated.
#' @param bw Bandwidth: a scalar or a vector of the same length as \code{ncol(x)}.
#' @param gaussian If TRUE, the normal (Gaussian) product kernel with full support is used, otherwise Epanechnikov.
#' @param method "R" or "CPP". If "CPP", then faster Rcpp functions are invoked, otherwise use a pure R implementation.
#' @param return.grid Logical: if true, then returns \code{xgrid} and appends the estimated density as the last column.
#'
#' @return A vector of density estimates evaluated at the grid points or, if \code{return.gric}, a matrix with the density in the last column.
#' @export
#'
#' @examples
kernelDensity <- function(x,
                          xgrid = NULL,
                          bw = NULL,
                          gaussian = TRUE,
                          method = "CPP",
                          return.grid = FALSE
) {
  one.dimensional <- is.vector(x)
  if (one.dimensional) {
    if (is.null(bw)) {
      bw <- bw.rot(x)
      warning(paste0("No bandwidth supplied, using Silverman's one-dimensional rule of thumb: bw = ", round(bw, 5), "."))
    }
    if (length(bw) != 1) stop("The bandwidth for one-dimensional smoothing must be a scalar.")
    if (is.null(xgrid)) xgrid <- x
    if (!is.vector(xgrid)) stop("The grid for one-dimensional smoothing must be a vector.")
    result <- switch(method,
      R = kernelDensityOneR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
      CPP = kernelDensityOneCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
    )
  } else {
    d <- ncol(x)
    if (is.data.frame(x)) x <- as.matrix(x)
    if (is.null(xgrid)) xgrid <- x
    if (is.data.frame(xgrid)) xgrid <- as.matrix(xgrid)
    if (d != ncol(xgrid)) stop("x and xgrid must be have the same number of columns (i.e. the same dimension).")
    if (is.null(bw)) {
      bw <- apply(x, 2, bw.rot)
      warning("No bandwidth supplied, using Silverman's one-dimensional rule of thumb in every dimension: bw = (", paste(round(bw, 5), collapse = ", "), ").")
    }
    if (length(bw) == 1) bw <- rep(bw, d)
    if (length(bw) != d) stop("The number of data dimensions is not the same as the number of bandwidths!")
    result <- switch(method,
      R = kernelDensityMultiR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
      CPP = kernelDensityMultiCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian)
    )
  }
  if (return.grid) result <- cbind(xgrid, result)
  return(result)
}


#' Local kernel smoother
#'
#' @param x A vector or a matrix of exogenous explanatory variables.
#' @param y A vector of dependent variable values.
#' @param xgrid A vector or a matrix with the same number of columns as \code{x} contatining the values at which to evaluate the smoother.
#' @param bw A scalar or a vector of smoothing bandwidths.
#' @param gaussian If \code{TRUE}, Gaussian kernel is used, otherwise Epanechnikov.
#' @param LOO If \code{TRUE}, the leave-one-out estimator is returned.
#' @param method \code{"R"} or \code{"CPP"}. The latter requires \code{Rcpp} compilation.
#' @param degree Integer: 0 for locally constant estimator (Nadaraya---Watson), 1 for locally linear (Cleveland's LOESS), 2 for locally quadratic (use with care)
#' @param standardise If \code{TRUE}, then the variables are centred to have mean 0 and re-scaled to have standard deviation 1.
#' @param trim Trimming function for small weights to speed up local regression (if \code{degree} is 1 or 2).
#' @param robust.iterations The number of robustifying iterations (due to Cleveland, 1979).
#' @param return.grid If \code{TRUE}, prepends \code{xgrid} to the return results.
#'
#' Standardisation is recommended for the purposes of numerical stability (sometimes \code{lm()} might choke when the dependent variable takes very large absolute values and its square is used).
#'
#' @return A vector of predicted values or, if \code{return.grid} is \code{TRUE}, a matrix with the predicted values in the last column.
#' @export
#'
#' @examples
kernelSmooth <- function(x,
                         y,
                         xgrid = NULL,
                         bw = NULL,
                         gaussian = TRUE,
                         LOO = FALSE,
                         method = "CPP",
                         degree = 0,
                         standardise = TRUE,
                         trim = function(x) 0.01 / length(x),
                         robust.iterations = 0,
                         return.grid = FALSE
) {
  one.dimensional <- is.vector(x)
  if (robust.iterations > 0 & !is.null(xgrid)) {
    warning("kernelSmooth: robust LOESS requested but a custom xgrid provided! Ignoring it.")
    xgrid <- NULL
  }
  if (LOO & !is.null(xgrid)) {
    warning("kernelSmooth: Leave-One-Out estimation requested, but a custom xgrid passed! Ignoring it.")
    xgrid <- NULL
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.null(xgrid)) xgrid <- x
  if (one.dimensional) {
    if (length(y) != length(x)) stop("x and y lengths differ.")
    if (is.null(bw)) {
      bw <- bw.rot(x)
      warning(paste0("No bandwidth supplied, using Silverman's one-dimensional rule of thumb: bw = ", round(bw, 5), "."))
    }
    if (length(bw) != 1) stop("The bandwidth for one-dimensional smoothing must be a scalar.")
    if (!is.vector(xgrid)) stop("The grid for one-dimensional smoothing must be a vector.")
    if (degree == 0) {
      result <- switch(method,
                       R   = kernelSmoothOneR(x = x, y = y, xgrid = xgrid, bw = bw, gaussian = gaussian, LOO = LOO),
                       CPP = kernelSmoothOneCPP(x = x, y = y, xgrid = xgrid, bw = bw, gaussian = gaussian, LOO = LOO))
    } else {
      K <- switch(method,
                  R = kernelWeightsOneR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
                  CPP = kernelWeightsOneCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian))
    }
  } else {
    if (length(y) != nrow(x)) stop("x and y lengths differ.")
    d <- ncol(x)
    if (is.data.frame(xgrid)) xgrid <- as.matrix(xgrid)
    if (d != ncol(xgrid)) stop("x and xgrid must be have the same number of columns (i.e. the same dimension).")
    if (is.null(bw)) {
      bw <- apply(x, 2, bw.rot)
      warning("No bandwidth supplied, using Silverman's one-dimensional rule of thumb in every dimension: bw = (", paste(round(bw, 5), collapse = ", "), ").")
    }
    if (length(bw) == 1) bw <- rep(bw, d)
    if (length(bw) != d) stop("The number of data dimensions is not the same as the number of bandwidths!")
    if (degree == 0) {
      result <- switch(method,
                       R   = kernelSmoothMultiR(x = x, y = y, xgrid = xgrid, bw = bw, gaussian = gaussian, LOO = LOO),
                       CPP = kernelSmoothMultiCPP(x = x, y = y, xgrid = xgrid, bw = bw, gaussian = gaussian, LOO = LOO))
    } else {
      K <- switch(method,
                  R = kernelWeightsMultiR(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian),
                  CPP = kernelWeightsMultiCPP(x = x, xgrid = xgrid, bw = bw, gaussian = gaussian))
    }
  }
  if (degree %in% 1:2) {
    if (LOO) diag(K) <- 0
    K <- K / rowSums(K)

    if (standardise) {
      if (one.dimensional) {
        m <- mean(x)
        s <- stats::sd(x)
        xs <- (x - m) / s
        xgrids <- (xgrid - m) / s
      } else {
        m <- colMeans(x)
        s <- apply(x, 2, stats::sd)
        xs <- sweep(sweep(x, 2, m, "-"), 2, s, "/")
        xgrids <- sweep(sweep(xgrid, 2, m, "-"), 2, s, "/")
      }
    } else {
      xs <- x
      xgrids <- xgrid
    }
    if (degree == 2) {
      xs <- cbind(xs, xs^2)
      xgrids <- cbind(xgrids, xgrids^2)
    }
    if (is.vector(xs)) { # For simplicity, now we need to use matrices in weighs OLS
      xs <- matrix(xs, ncol = 1)
      xgrids <- matrix(xgrids, ncol = 1)
    }
    coefhat <- t(apply(K, 1, function(w) {
      nonzw <- sparseVectorToList(w, trim = trim)
      # If there are no non-zero neighbours or one point has too much influence, declare failure, declare failure
      if (length(nonzw$idx) == 1 & any(nonzw$ct > 0.999)) return(rep(NA, d + 1))
      wreg <- sqrt(nonzw$ct)
      xw <- cbind(1, xs[nonzw$idx, ]) * wreg
      yw <- y[nonzw$idx] * wreg
      return(stats::.lm.fit(xw, yw)$coefficients)
    }))
    if (any(is.na(coefhat))) {
      bad.rows <- which(apply(coefhat, 1, function(x) any(is.na(x))))
      stop(paste0("Local weighted fit numerical problems: one point has a huge kernel weight > 0.999.
Potentially no neighbours with weights > trimming value, terminating to avoid a singular fit.
Problematic xgrid row indices: ", paste0(bad.rows, collapse = ", "), ".
Try increasing the bandwidth to get more data points in the vicinity!"))
    }
    result <- rowSums(coefhat * cbind(1, xgrids))
    if (robust.iterations > 0) {
      for (i in 1:robust.iterations) {
        resid <- y - result
        ss <- stats::median(abs(resid))
        resid <- resid / ss / 6
        deltak <- (1 - resid^2)^2 # Bisquare weights
        deltak[abs(resid) > 1] <- 0
        coefhat <- t(apply(K, 1, function(w) {
          nonzw <- sparseVectorToList(w, trim = trim)
          wreg <- sqrt(nonzw$ct) * sqrt(deltak[nonzw$idx])
          xw <- cbind(1, xs[nonzw$idx, ]) * wreg
          yw <- y[nonzw$idx] * wreg
          return(stats::.lm.fit(xw, yw)$coefficients)
        }))
        result <- rowSums(coefhat * cbind(1, xgrids))
      }
    }
  }

  if (return.grid) result <- cbind(xgrid, result)
  return(result)
}

#' Density and/or kernel regression estimator with conditioning on discrete variables
#'
#' @param x A vector or a matrix/data frame of discrete explanatory variables (exogenous). Non-integer values are fine because the data are split into bins defined by interactions of these variables.
#' @param y Optional: a vector of dependent variable values.
#' @param compact Logical: return unique values instead of full data with repeated observations?
#' @param fun A function that computes a statistic of \code{y} inside every category defined by \code{x}.
#'
#' @return A list with \code{x}, density estimator (\code{fhat}) and, if \code{y} was provided, regression estimate.
#' @export
#'
#' @examples
kernelDiscreteDensitySmooth <- function(x,
                                        y = NULL,
                                        compact = FALSE,
                                        fun = mean
) {
  if (is.matrix(x)) x <- as.data.frame(x)
  if (is.data.frame(x)) x <- as.integer(interaction(x, drop = TRUE, lex.order = TRUE))
  if (!is.vector(x)) stop("x must be a numeric vector, matrix, or a data frame!")
  n <- length(x)
  if (compact) {
    xtab <- table(x)
    fhat <- unname(xtab / n)
    xvals <- as.numeric(names(xtab))
  } else {
    fhat <- stats::ave(x, x, FUN = function(a) length(a) / n)
  }
  if (!is.null(y)) { # Smoothing y on x given the function
    if (!is.vector(y)) stop("x and y must be numeric vectors.")
    if (length(x) != length(y)) stop("x and y must be of equal length.")
    if (compact) {
      ys <- split(y, x)
      muhat <- lapply(ys, fun)
      return(list(x = xvals, y = unname(unlist(muhat)), fhat = fhat))
    } else {
      muhat <- stats::ave(y, x, FUN = fun)
      return(list(x = x, y = muhat, fhat = fhat))
    }
  } else {
    if (compact) return(x = xvals, fhat = fhat) else return(x = x, fhat = fhat)
  }
}

#' Density with conditioning on discrete and continuous variables
#'
#' @param x A numeric vector or matrix.
#' @param by A variable containing unique identifiers of discrete categories.
#' @param xgrid A numeric vector or a numeric matrix.
#' @param xgrid.by A variable containing unique identifiers of discrete categories of the grid.
#' @param bw Smoothing bandwidth(s) for continuous variables.
#' @param gaussian Passed to \code{kernelDensity},
#'
#' @return A numeric vector of the density estimate of the same length as \code{nrow(xgrid)}.
#' @export
#'
#' @examples
kernelMixedDensity <- function(x,
                               by,
                               xgrid = NULL,
                               xgrid.by = NULL,
                               bw = NULL,
                               gaussian = TRUE
) {
  if (is.null(xgrid)) xgrid <- x
  by <- as.character(by)
  xtab <- table(by)
  xgridtab <- table(xgrid.by)
  if (! names(xgridtab) %in% names(xtab))
  if (any(xtab<2)) stop("Some of the categories have only one observation, the distribution is degenerate. At least 2 observarions per category are needed.")
  n <- sum(xtab)
  res <- numeric(n)
  for (v in unique(names(xtab))) {
    s <- by==v
    prob <- sum(s)/n
    xsub <- x[s, ]
    xgridsub <- xgrid[s, ]
    res[s] <- kernelDensity(x=xsub, xgrid = xgridsub, bw = bw, gaussian = gaussian) * prob
  }
  return(res)
}

#' Density cross-validation
#'
#' @param x A numeric vector or matrix.
#' @param bw Candidate bandwidth values: scalar, vector, or a matrix (with columns corresponding to columns of \code{x}).
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#' @param gaussian Passed to \code{kernelWeights}.
#' @param method Passed to \code{kernelWeights}.
#'
#' @return A numeric vector of the same length as \code{bw} or \code{nrow(bw)}.
#'
#' @export
#' @examples
DCV <- function(x, bw, same = FALSE, gaussian = TRUE, method = "CPP") {
  if (!gaussian) stop("So far, I have programmed only Gaussian DCV.")
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (is.data.frame(bw)) bw <- as.matrix(bw)
  if (one.dim) {
    n <- length(x)
    many.bw <- (length(bw) > 1)
  } else {
    n <- nrow(x)
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else bw <- list(bw) # If there is only one bw, make it a list
  }
  CV <- function(b) {
    # A sub-function to compute the CV for one BW, parallelisable
    if (any(b <= 0)) return(Inf)
    if (!one.dim & length(b) == 1) b <- rep(b, ncol(x))
    K0 <- if (gaussian) stats::dnorm(0)^length(b) else 0.75^length(b) # Product kernel at 0
    KK <- kernelWeights(x = x, xgrid = x, bw = b*sqrt(2), gaussian = gaussian, method = method) # Easy Gaussian convolution!
    term1 <- sum(KK) / (n^2 * prod(b))
    # Computing the LOO estimator efficiently: fhat_i(x) = n/(n-1) * fhat(x) - 1/((n-1)*b^s) * K((X[i] - x)/b)
    fhat <- kernelDensity(x = x, xgrid = x, bw = b, gaussian = gaussian, method = method)
    fhat.LOO <- (n * fhat - K0 / (prod(b))) / (n - 1)
    term2 <- -2 * mean(fhat.LOO)
    return(term1 + term2)
  }
  CV.values <- lapply(bw, CV)
  return(unlist(CV.values))
}

#' Least-squares cross-validation function for the Nadaraya-Watson estimator
#'
#' @param x A numeric vector or matrix.
#' @param y A vector of dependent variable values.
#' @param bw Candidate bandwidth values: scalar, vector, or a matrix (with columns corresponding to columns of \code{x}).
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#' @param degree Passed to \code{kernelSmooth}.
#' @param gaussian Passed to \code{kernelSmooth}.
#' @param robust.iterations Passed to \code{kernelSmooth}.
#' @param method Passed to \code{kernelSmooth}.
#'
#' @return A numeric vector of the same length as \code{bw} or \code{nrow(bw)}.
#' @export
#'
#' @examples
LSCV <- function(x, y, bw, same = FALSE, degree = 0, gaussian = TRUE, robust.iterations = 0, method = "CPP") {
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (is.data.frame(bw)) bw <- as.matrix(bw)
  if (one.dim) {
    n <- length(x)
    many.bw <- (length(bw) > 1)
  } else {
    n <- nrow(x)
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else bw <- list(bw) # If there is only one bw, make it a list
  }
  ASE <- function(b) {
    if (any(b <= 0)) return(Inf)
    muhat_i <- kernelSmooth(x = x, y = y, bw = b, gaussian = gaussian, degree = degree, method = method, LOO = TRUE, robust.iterations = robust.iterations)
    m <- mean((y - muhat_i)^2)
    if (!is.finite(m)) m <- Inf
    return(m)
  }
  ASE.values <- lapply(bw, ASE)
  return(unlist(ASE.values))
}

#' Bandwidth Selectors for Kernel Density Estimation
#'
#' Finds the optimal bandwidth by minimising the density cross-valication or least-squares criteria.
#' Remember that since usually, the CV function is highly non-linear, the return value should be taken with a grain of salt.
#' With non-smooth kernels (such as uniform), it will oftern return the local minimum after starting from a reasonable value.
#' The user might want to standardise the input matrix \code{x} by column (divide by some estimator of scale, like \code{sd}
#' or \code{IQR}) and examine the behaviour of the CV criterion as a function of unique bandwidth (\code{same} argument).
#' If it seems that the optimum is unique, then they may proceed by multiplying the bandwidth by the scale measure,
#' and start the search for the optimal bandwidth in multiple dimensions.
#'
#' @param x A numeric vector or numeric matrix.
#' @param y A numeric vector of responses (dependent variable) if \code{CV = "LSCV"}.
#' @param gaussian Passed to \code{kernelDensity} or to \code{kernelSmooth}.
#' @param robust.iterations Passed to \code{kernelSmooth} if \code{CV = "LSCV"}.
#' @param method Passed to \code{kernelDensity} or to \code{kernelSmooth}.
#' @param degree Passed to \code{kernelSmooth} if \code{CV = "LSCV"}.
#' @param start.bw Initial value for bandwidth search.
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#' @param CV "DCV" for density cross-validation, "LSCV" for least-squares cross-validation.
#' @param opt.fun Optimiser name (can be a custom function).
#' @param ret.fun A function that extract the optimal parameter from the output of the optimiser.
#' @param par.name.in.opt The name of the initial parameter vector in the output of the optimiser.
#' @param fun.name.in.opt The name of the criterion in the output of the optimiser (to be minimised).
#' @param report.code Logical: print out the optimiser return code for diagnostics?
#' @param ... Other parameters passed to the optimiser.
#'
#' @return Numeric vector or scalar of the optimal bandwidth.
#' @export
#'
#' @examples
bw.CV <- function(x, y = NULL, gaussian = TRUE, robust.iterations = 0, method = "CPP", degree = 0, start.bw = NULL, same = FALSE, CV = c("DCV", "LSCV"),
                  opt.fun = c("nlm", "optim", "nlminb", "optimise"),
                  ret.fun = NULL,
                  par.name.in.opt = NULL,
                  fun.name.in.opt = NULL,
                  report.code = FALSE,
                  ...) {
  opt.fun <- opt.fun[1]
  CV <- CV[1]
  one.dim <- is.vector(x) # Are our data one-dimensional?
  opt <- get(opt.fun)
  f.to.min <- if (CV == "DCV") function(b) DCV(x = x, bw = b, gaussian = gaussian, same = same, method = method) else if (CV == "LSCV") function(b) LSCV(x = x, y = y, bw = b, same = same, degree = degree, gaussian = gaussian, robust.iterations = robust.iterations, method = method) else stop("bw.CV: 'CV' should be either 'DCV' or 'LSCV'!")
  if (is.null(ret.fun)) ret.fun <- switch(opt.fun, nlm = function(x) x[["estimate"]],
                                          optim = function(x) x[["par"]],
                                          optimise = function(x) x[["minimum"]],
                                          nlminb = function(x) x[["par"]])
  if (is.null(start.bw)) start.bw <- bw.rot(x)
  if (same) start.bw <- mean(start.bw)
  opt.result <- switch(opt.fun,
                       nlm = suppressWarnings(stats::nlm(f = f.to.min, p = start.bw, ...)),
                       optim = stats::optim(par = start.bw, fn = f.to.min, ...),
                       optimise = stats::optimise(f = f.to.min, ...),
                       nlminb = stats::nlminb(start = start.bw, objective = f.to.min, lower = 1e-8, ...))
  if (report.code) switch(opt.fun,
         nlm = message(paste0("nlm exit code ", opt.result$code, ", ||gradient||=", round(sqrt(sum(opt.result$gradient^2)), 6), ", done in ", opt.result$iterations, " iterations.")),
         optim = message(paste0("optim exit code ", opt.result$convergence, ", done in (", paste0(opt.result$counts, collapse = ", "), ") iterations.")),
         optimise = message(paste0("optimise does not return useful convergence information.")),
         nlminb = message(paste0("nlminb exit code ", opt.result$convergence, ", done in ", paste0(opt.result$iterations, collapse = ", "), " iterations.")))
  if (is.null(opt.result)) {
    arg.list <- list()
    if (!is.null(par.name.in.opt)) arg.list[[par.name.in.opt]] <- start.bw
    if (!is.null(fun.name.in.opt)) arg.list[[fun.name.in.opt]] <- f.to.min
    arg.list <- c(arg.list, list(...))
    opt.result <- do.call(opt, args = arg.list)
  }
  return(ret.fun(opt.result))
}


#' Basic univatiate kernel functions
#'
#' Computes 5 most
#'  popular kernel functions of orders 2, 4, and 6 with the potential of returning an analytical convolution kernel for density cross-validation.
#'
#' @param x A numeric vector of values at which to compute the kernel function.
#' @param kernel Kernel type: uniform, Epanechnikov, triangular, quartic, or Gaussian.
#' @param order Kernel order. 2nd-order kernels are always non-negative. kth-order kernels have all moments from 1 to (k-1) equal to zero, which is achieved by having some negative values.
#' @param rescale Logical: rescale to unit variance? If \code{TRUE}, ensures that, for the chosen kernel
#' name, the second-order kernel integrates to 1:
#' \eqn{\int_{-\infty}^{+\infty} x^2 k(x) = \sigma^2_k = 1}.
#' This is useful because in this case, the constant \code{k_2} in formulæ 3.12 and 3.21
#' from \insertCite{silverman1986density;textual}{smoothemplik} is equal to 1.
#' @param convolution Logical: return the convolution kernel? (Useful for density cross-validation.)
#'
#' @return A numeric vector of the same length as input.
#' @importFrom Rdpack reprompt
#' @export
#' @references
#' \insertCite{silverman1986density}{smoothemplik}
#'
#'
#'
#' @examples
#' ks <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian"); names(ks) <- ks
#' os <- c(2, 4, 6); names(os) <- paste0("o", os)
#' cols <- c("#000000CC", "#0000CCCC", "#CC0000CC", "#00AA00CC", "#BB8800CC")
#' put.legend <- function() legend("topright", legend = ks, lty = 1, col = cols, bty = "n")
#' xgrid <- seq(-4, 4, length.out = 301)
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(0, 1.1),
#'   xlab = "", ylab = "", main = "Unscaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], rescale = FALSE), col = cols[i])
#' par(mfrow = c(1, 2))
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.1, 0.8), xlab = "", ylab = "",
#'   main = "4th-order scaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = 4), col = cols[i])
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.25, 1.2), xlab = "", ylab = "",
#'   main = "6th-order scaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = 6), col = cols[i])
#' par(mfrow = c(1, 1))
#' plot(NULL, NULL, xlim = range(xgrid), ylim = c(-0.25, 1.4), xlab = "", ylab = "",
#'   main = "Convolution kernels", bty = "n"); put.legend()
#' for (i in 1:5) {
#'   for (j in 1:3) lines(xgrid, kernelFun(xgrid, kernel = ks[i], order = os[j],
#'   convolution = TRUE), col = cols[i], lty = j)
#' }; legend("topleft", c("2nd order", "4th order", "6th order"), lty = 1:3, bty = "n")
#'
#' # All kernels integrate to correct values; we compute the moments
#' mom <- Vectorize(function(k, o, m, c) integrate(function(x) x^m * kernelFun(x, k, o,
#'   rescale = FALSE, convolution = c), lower = -Inf, upper = Inf)$value)
#' for (m in 0:6) {
#'   cat("\nComputing integrals of x^", m, " * f(x). \nSimple unscaled kernel:\n", sep = "")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = FALSE)), 4))
#'   cat("Convolution kernel:\n")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = TRUE)), 4))
#' }
#'
kernelFun <- function(x,
                      kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                      order = c(2, 4, 6),
                      rescale = TRUE,
                      convolution = FALSE
) {
  order <- order[1]
  if (!(order %in% c(2, 4, 6))) stop("Only kernels of orders 2, 4, 6 have been implemented.")
  kernel <- kernel[1]
  adj.factor <- switch(kernel,
    uniform = sqrt(3),
    triangular = sqrt(6),
    epanechnikov = sqrt(5),
    quartic = sqrt(7),
    gaussian = 1
  )
  if (!rescale) adj.factor <- 1
  x <- x / adj.factor
  x <- abs(x)
  if (!convolution) {
    abc <- if (order == 2) c(1, 0, 0) else if (order == 4) {
      switch(kernel,
             uniform = c(NA, NA, NA),
             triangular = c(12, -30, 0) / 7,
             epanechnikov = c(15, -35, 0) / 8,
             quartic = c(7, -21, 0) / 4,
             gaussian = c(3, -1, 0) / 2
      )
    } else if (order == 6) {
      switch(kernel,
             uniform = c(NA, NA, NA),
             triangular = c(1635, -10500, 11970) / 683,
             epanechnikov = c(175, -1050, 1155) / 64,
             quartic = c(315, -2310, 3003) / 128,
             gaussian = c(15, -10, 1) / 8
      )
    }
    k <- switch(kernel,
      uniform = 1/2 * (x < 1) * (order == 2) + (-1/6 * (x < 1) + 4/3 * (x < 0.5)) * (order == 4) + (1/20 * (x < 1) - 9/20 * (x < 2/3) + 9/4 * (x < 1/3)) * (order == 6),
      triangular = (1 - x) * (x < 1),
      epanechnikov = 3/4 * (1 - x^2) * (x < 1),
      quartic = 15/16 * (1 - x^2)^2 * (x < 1),
      gaussian = stats::dnorm(x)
    )
    if (order > 2 & kernel != "uniform") {
      polynomial <- abc[1] + abc[2] * x^2 + abc[3] * x^4
      k <- k * polynomial
    }
  } else { # Convolution kernel
    if (order == 2) {
      k <- switch(kernel,
                  uniform = 1 / 4 * (2 - x) * (x < 2),
                  triangular = 1 / 6 * ((3 * x^3 - 6 * x^2 + 4) * (x <= 1) + (8 - 12 * x + 6 * x^2 - x^3) * (x > 1 & x < 2)),
                  epanechnikov = 3 / 160 * (2 - x)^3 * (x^2 + 6*x + 4) * (x < 2),
                  quartic = 5 / 3584 * (2 - x)^5 * (16 + (2*x + x^2)*(20 + 8*x + x^2)) * (x < 2),
                  gaussian = stats::dnorm(x, sd = sqrt(2))
      )
    } else if (order == 4) {
    k <- switch(kernel,
                uniform = 1/36*(64*(1 - x)*(x < 1) - 16*(1*(2*x < 1) + (3/2 - x)*(1/2 <= x & x < 3/2)) + (2 - x)*(x<2)),
                triangular = 3/343 * ((152 + x^2*(-616 + x*(238 + x*(560 + x*(-322 + 5*x*(-14 + 9*x))))))*(x<1) + (2-x)^3 * (30 + x*(-4 + x*(-74 + 5*x*(4 + 3*x))))*(x>=1 & x<2)),
                epanechnikov = 1/2048 * (5*(2-x)^3*(64 + x*(96 + x*(-144 + x*(-160 + x*(48 + 7*x*(6 + x))))))) * (x < 2),
                quartic = 1/2342912 * (35*(2 - x)^5*(2944 + x*(7360 + x*(-1440 + x*(-18320 + x*(-9896 + 9*x*(560 + x*(536 + 15*x*(10 + x))))))))) * (x < 2),
                gaussian = 1/64 * stats::dnorm(x, sd = sqrt(2)) * (108 - 28*x^2 + x^4)
    )
  } else if (order == 6) {
    k <- switch(kernel,
                uniform = 1/400*(956 - 2107*x)*(x<1/3) + 1/400*(680 - 1279*x)*(x>=1/3 & x<2/3) + (-61/40 + 41/25*x)*(x>=2/3 & x<1) + (1/2 - 77/200*x)*(x>=1 & x<4/3) + 1/400*(-28 + 17*x)*(x>=4/3 & x<5/3) + (2-x)/400*(x>=5/3 & x<2),
                triangular = 1/466489*((15/22)*(1353180 + x^2*(-10733690 + x*(3663605 + 2*x*(11535524 + x*(-6217288 + x*(-8435812 + 3*x*(1846130 + 133*x*(4400 + x*(-3168 + 19*x*(-22 + 15*x))))))))))*(x<1) + (-(15/22))*(-2 + x)^3*(240311 + 2*x*(-24496 + x*(-770780 + x*(257608 + x*(1018338 + 133*x*(-3072 + x*(-2180 + 57*x*(8 + 5*x))))))))*(x>=1 & x<2)),
                epanechnikov = -((105*(-61440 + x^2*(465920 + x*(-232960 + x*(-838656 + x*(640640 + x*(439296 - 486720*x + 84448*x^3 - 9828*x^5 + 495*x^7)))))))/3407872) * (x < 2),
                quartic = (-((315*(-2 + x)^5*(380928 + x*(952320 + x*(-1739776 + x*(-6254080 + x*(478464 + x*(9024512 + x*(2918912 + x*(-3982464 + x*(-2272576 + 143*x*(1000 + x*(3012 + 91*x*(10 + x)))))))))))))/1853882368)) * (x < 2),
                gaussian = stats::dnorm(x, sd = sqrt(2)) * (36240 - 19360*x^2 + 2312*x^4 - 88*x^6 + x^8)/16384
    )
  }
  }
  k <- k / adj.factor
  return(k)
}

