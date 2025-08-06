#' Silverman's rule-of-thumb bandwidth
#'
#' A fail-safe function that would return a nice Silverman-like bandwidth
#'   suggestion for data for which the standard deviation might be NA or 0.
#'
#' \eqn{\Sigma = \mathrm{\mathop{diag}}(\sigma^2_k)}{\Sigma = diag(\sigma^2_k)}
#' with
#' \eqn{\det\Sigma = \prod_k \sigma^2_k}{det \Sigma = prod(\sigma^2_k)}
#' and
#' \eqn{\Sigma^{-1} = \mathrm{\mathop{diag}}(1/\sigma^{2}_k)}{\Sigma = diag(1/\sigma^2_k)}).
#' Then, the formula 4.12 in Silverman (1986) depends only on \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta}.
#' \eqn{\alpha = \mathrm{\mathop{diag}}(\sigma^2_k)}{\Sigma = diag(\sigma^2_k)}
#' (which depend only on the kernel and are fixed for a multivariate normal), and on the L2-norm of the
#' second derivative of the density. The (i, i)th element of the Hessian of multi-variate normal
#' (\eqn{\phi(x_1, \ldots, x_d) = \phi(X)}{\phi(x_1, ..., x_d) = \phi(X)}) is
#' \eqn{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}.
#'
#' @param x A numeric vector without non-finite values.
#' @param kernel A string character: \code{"gaussian"}, \code{"uniform"},
#'   \code{"triangular"}, \code{"epanechnikov"}, or  \code{"quartic"}.
#' @param na.rm Logical: should missing values be removed? Setting it to TRUE
#'   may cause issues because variable-wise removal of NAs may return a
#'   bandwidth that is inappropriate for the final data set for which it is
#'   suggested.
#' @param robust Logical: safeguard against extreme observations? If TRUE, uses
#'   \code{min(sd(x), IQR(x)/1.34)} to estimate the spread.
#' @param discontinuous Logical: if the true density is discontinuous (i.e. has
#'   jumps), then, the formula for the optimal bandwidth for density estimation changes.
#'
#' @details
#' The rule-of-thumb bandwidth is obtained under the assumption that the true
#' density is multivariate normal with zero covariances
#' (i.e. a diagonal variance-covariance matrix). For details,
#' see \insertCite{silverman1986density}{smoothemplik}.
#'
#' @return A numeric vector of bandwidths that are a reasonable start optimal non-parametric density estimation of \code{x}.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' set.seed(1); bw.rot(stats::rnorm(100)) # Should be 0.3787568 in R version 4.0.4
#' set.seed(1); bw.rot(matrix(stats::rnorm(500), ncol = 10)) # 0.4737872 ... 0.7089850
#' @export
bw.rot <- function(x, kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                   na.rm = FALSE, robust = TRUE, discontinuous = FALSE) {
  kernel <- kernel[1]
  if (!(kernel %in% c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"))) stop("bw.rot: Wrong kernel type.")
  if (any(is.na(x))) {
    if (na.rm) warning("bw.rot: the input data contain missing values. Do something about them because proper analysis is impossible with NA in 'x'.") else
      stop("bw.rot: There are missing values in the data, but non-parametric methods rely on data with finite numeric values only.")
  }
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) x <- matrix(x, ncol = 1)
  d <- ncol(x)
  n <- nrow(x)
  s <- apply(x, 2, function(x) stats::sd(x, na.rm = na.rm))
  if (robust) {
    sr <- apply(x, 2, function(x) stats::IQR(x, na.rm = na.rm)) / 1.34898
    # Some IQRs can be zeros if the variable is dummy; we want to specifically avoid that
    gt0 <- sr > 0
    if (any(gt0)) s[gt0] <- pmin(s[gt0], sr[gt0])
  }
  vk <- switch(kernel, gaussian = 1, uniform = 1/3, triangular = 1/6, epanechnikov = 1/5, quartic = 1/7) # Variance of the kernel
  rk <- switch(kernel, gaussian = 1/sqrt(4*pi), uniform = 1/2, triangular = 2/3, epanechnikov = 3/5, quartic = 5/7)^d # Roughness of the kernel
  rdnorm2 <- (0.5*d + 0.25*d^2) / (2*sqrt(pi))^d
  p <- 1 / (d+4)
  AK <- (d*rk / vk^2 / rdnorm2)^p # (4.15 from Silverman, 1986)

  if (!discontinuous) {
    if (any(!is.finite(s))) {
      stop("bw.rot: Could not compute the bandwidth; check your data -- most likely it has fewer than 2 observations.")
    } else if (all(s > 0)) { # Positive variance = at least two points
      return(AK * s * n^(-p))
    } else {
      return(rep(1, d))
    }
  } else {
    stop("Discontinuous density not implemented yet.")
  }

}

#' Probability integral transform
#'
#' @param x A numeric vector of data points.
#' @param xout A numeric vector. If supplied, then the transformed function at the grid points different from \code{x} takes values
#' equidistant between themselves and the ends of the interval to which they belong.
#'
#' @return A numeric vector of values strictly between 0 and 1 of the same length as \code{xout} (or \code{x}, if \code{xout} is \code{NULL}).
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
pit <- function(x, xout = NULL) {
  x.transformed <- stats::ecdf(x)(x) - 0.5 / length(x)
  if (is.null(xout)) {
    return(x.transformed)
  } else {
    if (isTRUE(all.equal(x, xout, tolerance = .Machine$double.eps))) return(x.transformed)
    x.uniq.sorted <- sort(unique(x))
    ecdf.uniq.sorted <- sort(x.transformed[!duplicated(x)])
    xout.uniq.sorted <- sort(unique(xout))
    xout.cut <- cut(xout, breaks = c(-Inf, x.uniq.sorted, Inf))
    exact <- xout %in% x.uniq.sorted
    xout.cut[xout %in% x.uniq.sorted] <- NA # Exact matches will be dealt with at the last step
    xout.list <- split(xout, xout.cut)
    ecdf.uniq01 <- c(0, ecdf.uniq.sorted, 1)
    xout.list.uniq <- lapply(xout.list, function(x) sort(unique(x)))
    xout.spaced <- lapply(seq_len(length(xout.list.uniq)), function(i) {
      xg <- xout.list.uniq[[i]]
      l <- length(xg)
      if (l > 0) return(seq(ecdf.uniq01[i], ecdf.uniq01[i+1], length.out = l+2)[c(-1, -l-2)]) else return(NULL)
    })
    xout.spaced <- unlist(xout.spaced)
    xout.uniq.sorted.noorig <- xout.uniq.sorted[!(xout.uniq.sorted %in% x)]
    ret <- xout.spaced[match(xout, xout.uniq.sorted.noorig)]
    ret[exact] <- x.transformed[match(xout[exact], x)]
    return(ret)
  }
}

#' @importFrom data.table := .GRP .N
.deduplicate <- function(x) {
  x <- data.table::as.data.table(as.data.frame(x))
  dup.rate <- mean(duplicated(x))
  x[, `:=`(`(id)` = .GRP, `(count)` = .N), by = names(x)]
  # The penultimate column, (id), contains the id of the unique combination;
  # The last one, `(count)`, contains the counts
  id <- x[["(id)"]]
  x <- unique(x, by = "(id)")
  count <- x[["(count)"]]
  x[, c("(id)", "(count)") := NULL]
  x <- as.matrix(x)
  if (ncol(x) == 1) colnames(x) <- NULL
  return(list(x = x, id = id, count = count, dup.rate = dup.rate))
}

#' Check the data for kernel estimation
#'
#' @param x A numeric vector, matrix, or data frame containing observations. For density, the
#'   points used to compute the density. For kernel regression, the points corresponding to
#'   explanatory variables.
#' @param y Optional: a vector of dependent variable values.
#' @param xout A vector or a matrix of data points with \code{ncol(xout) = ncol(x)}
#'   at which the user desires to compute the weights, density, or predictions.
#'   In other words, this is the requested evaluation grid.
#'   If \code{NULL}, then \code{x} itself is used as the grid.
#' @param weights A numeric vector of observation weights (typically counts) to
#'   perform weighted operations. If null, \code{rep(1, NROW(x))} is used. In
#'   all calculations, the total number of observations is assumed to be the
#'   sum of \code{weights}.
#' @param bw Bandwidth for the kernel: a scalar or a vector of the same length as \code{ncol(x)}.
#'   Since it is the crucial parameter in many applications, a warning is thrown if the bandwidth
#'   is not supplied, and then, Silverman's rule of thumb (via \code{bw.row()}) is applied
#'   to *every dimension* of \code{x}.
#' @param kernel Character describing the desired kernel type. NB: due to limited machine precision, even Gaussian has finite support.
#' @param order An integer: 2, 4, or 6. Order-2 kernels are the standard kernels that
#'   are positive everywhere. Orders 4 and 6 produce some negative values, which reduces bias but may hamper density estimation.
#' @param convolution Logical: if FALSE, returns the usual kernel. If TRUE, returns
#'   the convolution kernel that is used in density cross-validation.
#' @param sparse Logical: TODO (ignored)
#' @param deduplicate.x Logical: if TRUE, full duplicates in the input \code{x}
#'   and \code{y} are counted and transformed into weights; subsetting indices
#'   to reconstruct the duplicated data set from the unique one are also returned.
#' @param deduplicate.xout Logical: if TRUE, full duplicates in the input \code{xout}
#'   are counted; subsetting indices to reconstruct the duplicated data set from
#'   the unique one are returned.
#' @param no.dedup Logical: if TRUE, sets \code{deduplicate.x} and \code{deduplicate.xout}
#'   to FALSE (shorthand).
#' @param PIT If TRUE, the Probability Integral Transform (PIT) is applied to all columns
#'   of \code{x} via \code{ecdf} in order to map all values into the [0, 1] range. May
#'   be an integer vector of indices of columns to which the PIT should be applied.
#'
#' @description
#' Checks if the order is 2, 4, or 6, transforms the objects into matrices,
#' checks the dimensions, provides the bandwidth, creates default arguments
#' to pass to the C++ functions, carries out de-duplication for speed-up etc.
#'
#' @return A list of arguments that are taken by [kernelDensity()] and [kernelSmooth()].
#' @export
#'
#' @examples
#' # De-duplication facilities
#' set.seed(1)  # Creating a data set with many duplicates
#' n.uniq <- 10000
#' n <- 60000
#' inds <- ceiling(runif(n, 0, n.uniq))
#' x.uniq <- matrix(rnorm(n.uniq*10), ncol = 10)
#' x <- x.uniq[inds, ]
#' y <- runif(n.uniq)[inds]
#' xout <- x.uniq[ceiling(runif(n.uniq*3, 0, n.uniq)), ]
#' w <- runif(n)
#' print(system.time(a1 <- prepareKernel(x, y, xout, w, bw = 0.5)))
#' print(system.time(a2 <- prepareKernel(x, y, xout, w, bw = 0.5,
#'                   deduplicate.x = FALSE, deduplicate.xout = FALSE)))
#' print(c(object.size(a1), object.size(a2)) / 1024) # Kilobytes used
#' # Speed-memory trade-off: 4 times smaller, takes 0.2 s, but reduces the
#' # number of matrix operations by a factor of
#' 1 - prod(1 - a1$duplicate.stats[1:2])    # 95% fewer operations
#' sum(a1$weights) - sum(a2$weights)  # Should be 0 or near machine epsilon
prepareKernel <- function(x,
                           y = NULL,
                           xout = NULL,
                           weights = NULL,
                           bw = NULL,
                           kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                           order = 2,
                           convolution = FALSE,
                           sparse = FALSE,
                           deduplicate.x = TRUE,
                           deduplicate.xout = TRUE,
                           no.dedup = FALSE,
                           PIT = FALSE
) {
  kernel <- kernel[1]
  if (!(order %in% c(2, 4, 6))) stop("The kernel order muse be 2, 4, or 6.")
  if (convolution && order > 2) stop("At this moment, convolution kernels have been implemented for kernel order 2 only.")
  if (no.dedup) deduplicate.x <- deduplicate.xout <- FALSE

  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.vector(x)) x <- matrix(x, ncol = 1) # The C++ code is equally fast for vectors and matrices
  if (is.null(weights)) weights <- rep(1, NROW(x))
  if (length(weights) != NROW(x)) stop("The length of 'weights' must be equal to the number of obervations (NROW(x)).")
  if (!is.null(y)) {
    if (NCOL(y) > 1) stop("y must be a numeric vector.")
    if (!is.null(dim(y))) y <- drop(y)
    if (length(y) != nrow(x)) stop("The length of 'y' must be equal to the number of obervations (NROW(x)).")
  }
  n <- NROW(x)

  x.matches <- xout.matches <- NULL
  duplicate.stats <- c(dup.rate.x = NA, dup.rate.xout = NA, seconds.x = NA, seconds.xout = NA)

  xout.not.given <- FALSE
  if (is.null(xout)) {
    xout.not.given <- TRUE
    xout <- x # Creating a placeholder that undergoes only class changes, but not calculations
    deduplicate.xout <- deduplicate.x
  }
  if (is.data.frame(xout)) xout <- as.matrix(xout)
  if (is.vector(xout)) xout <- matrix(xout, ncol = 1)


  d <- ncol(x) # Dimension of the problem
  if (d != ncol(xout)) stop("x and xout must be have the same number of columns (i.e. the same dimension).")

  # If xout is not given, de-duplicate x first, and copy everything to xout
  if (deduplicate.x) {
    xy <- x
    if (!is.null(y)) xy <- cbind(x, `(y)` = y)
    tic0 <- Sys.time()
    xy <- .deduplicate(xy)
    duplicate.stats[1] <- xy$dup.rate
    duplicate.stats[3] <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
    # If there are any duplicates, add up their weights and remove duplicates
    if (duplicate.stats[1] > 0) {
      x.matches <- xy$id
      if (!is.null(y)) {
        y <- xy$x[, ncol(xy$x)]
        x <- xy$x[, -ncol(xy$x), drop = FALSE]
      } else {
        x <- xy$x
      }
      weights <- unname(vapply(split(weights, x.matches), sum, FUN.VALUE = numeric(1))) # Adding up weights
    } else {
      deduplicate.x <- FALSE # No duplicates found = no action needed
    }
  }

  if (xout.not.given) {
    deduplicate.xout <- deduplicate.x # Copy the attributes
    xout.matches <- x.matches
    xout <- x
    duplicate.stats[2] <- duplicate.stats[1]
    duplicate.stats[3:4] <- duplicate.stats[3] * 0.5 # The time is considered to be evenly split between tasks
  } else { # If xout is given
    if (deduplicate.xout) {
      # If there are any grid duplicates, carry out similar procedures
      tic0 <- Sys.time()
      xy <- .deduplicate(xout)
      duplicate.stats[2] <- xy$dup.rate
      duplicate.stats[4] <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
      if (duplicate.stats[2] > 0) {
        xout.matches <- xy$id
        xout <- xy$x
      } else {
        deduplicate.xout <- FALSE
      }
    }
  }
  nout <- NROW(xout)

  if (PIT) {
    for (i in 1:d) {
      xout[, i] <- pit(x = x[, i], xout = xout[, i])
      x[, i] <- pit(x[, i])
    }
  }

  if (is.null(bw)) {
    bw <- bw.rot(x, kernel = kernel)
    warning("No bandwidth supplied, using Silverman's multi-dimensional rule of thumb: bw = (",
            paste(sprintf("%1.2e", bw), collapse = ", "), ").")
  }
  if (length(bw) == 1) {
    bw <- matrix(bw, nrow = nout, ncol = d)
  } else if (length(bw) == nout) {
    bw <- matrix(rep(bw, d), nrow = nout, ncol = d)
  } else if (length(bw) == d) {
    bw <- t(replicate(nout, bw))
  } else if (length(bw) != d*nout) {
    stop("The bandwidths must be a vector of length 1, or ncol(xout), or nrow(xout), or a matrix with the same dimensions as xout.")
  }

  return(list(x = x, y = y, xout = xout, weights = weights,
              order = order, bw = bw, kernel = kernel,
              deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout,
              x.matches = x.matches, xout.matches = xout.matches,
              duplicate.stats = duplicate.stats))
}


#' Kernel-based weights
#'
#' @inheritParams prepareKernel
#' @param sparse Logical: TODO (should be ignored?)
#'
#' Note that if \code{pit = TRUE}, then the kernel-based weights become nearest-neighbour weights (i.e. not much different from the ones used
#' internally in the built-in \code{loess} function) since the distances now depend on the ordering of data, not the values per se.
#'
#' Technical remark: if the kernel is Gaussian, then, the ratio of the tail density
#' to the maximum value (at 0) is less than mach.eps/2 when abs(x) > 2*sqrt(106*log(2)) ~ 8.572.
#' This has implications the relative error of the calculation: even the
#' kernel with full support (theoretically) may fail to produce numerically distinct
#' values if the argument values are more than ~8.5 standard deviations away from the mean.
#'
#' @return A matrix of weights of dimensions nrow(xout) x nrow(x).
#' @export
#'
#' @importClassesFrom Matrix dgCMatrix
#' @examples
#' set.seed(1)
#' x   <- sort(rnorm(1000)) # Observed values
#' g   <- seq(-10, 10, 0.1) # Grid for evaluation
#' w   <- kernelWeights(x, g, bw = 2, kernel = "triangular")
#' wsp <- kernelWeights(x, g, bw = 2, kernel = "triangular", sparse = TRUE)
#' print(c(object.size(w), object.size(wsp)) / 1024) # Kilobytes used
#' image(g, x, w)
#' all.equal(w[, 1],  # Internal calculation for one column
#'             kernelFun((g - x[1])/2, "triangular", 2, FALSE))
#'
#' # Bare-bones interface to the C++ functions
#' # Example: 4th-order convolution kernels
#' x <- seq(-3, 5, length.out = 301)
#' ks <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian")
#' kmat <- sapply(ks, function(k) kernelFun(x, k, 4, TRUE))
#' matplot(x, kmat, type = "l", lty = 1, bty = "n", lwd = 2)
#' legend("topright", ks, col = 1:5, lwd = 2)
kernelWeights <- function(x,
                          xout = NULL,
                          bw = NULL,
                          kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                          order = 2,
                          convolution = FALSE,
                          sparse = FALSE,
                          PIT = FALSE,
                          deduplicate.x = FALSE,
                          deduplicate.xout = FALSE,
                          no.dedup = FALSE
) {
  arg <- prepareKernel(x = x, xout = xout, weights = NULL, bw = bw, kernel = kernel, PIT = PIT, order = order, convolution = convolution,
                        deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout, no.dedup = no.dedup)
  if (!sparse)
    result <- kernelWeightsCPP(x = arg$x, xout = arg$xout, bw = arg$bw, kernel = arg$kernel, order = arg$order, convolution = convolution) else
      result <- sparseKernelWeightsCPP(x = arg$x, xout = arg$xout, bw = arg$bw, kernel = arg$kernel, order = arg$order, convolution = convolution)
    # Checking if de-duplication was done and re-populating the full matrix if necessary
    if (arg$deduplicate.x) result <- result[, arg$x.matches]
    if (arg$deduplicate.xout) result <- result[arg$xout.matches, ]
    return(result)
}


checkDupStat <- function(ds, thresh = 0.1) {
  # Here, isTRUE is needed because without de-duplication, the values of ds are NA
  if (isTRUE(sum(ds[1:2]) == 0) & sum(ds[3:4]) > 0.1) {  # Only substantial slow-downs are reported
    message("No duplicates found in input 'x' and 'xout'. Consider setting 'no.dedup = TRUE' to eliminate the overhead (> 100 ms).")
  } else if (isTRUE(ds[1] == 0) & ds[3] > 0.1) {
    message("No duplicates found in input 'x'. Consider setting 'deduplicate.x = FALSE' to eliminate the overhead (> 100 ms).")
  } else if (isTRUE(ds[2] == 0) & ds[4] > 0.1) {
    message("No duplicates found in input 'xout'. Consider setting 'deduplicate.xout = FALSE' to eliminate the overhead (> 100 ms).")
  }
  return(NULL)
}


#' Kernel density estimation
#'
#' @inheritParams kernelWeights
#' @inheritParams prepareKernel
#' @param chunks Integer: the number of chunks to split the task into (limits
#'   RAM usage but increases overhead). \code{0} = auto-select (making sure that
#'   no matrix has more than 2^27 elements).
#' @param return.grid Logical: if \code{TRUE}, returns \code{xout} and appends the estimated density as the last column.
#'
#' The number of chunks for kernel density and regression estimation is chosen
#' in such a manner that the number of elements in the internal weight matrix
#' should not exceed \eqn{2^{27} = 1.3\cdot 10^8}{2^{27} = 1.3e+8}, which caps
#' RAM use (64 bits = 8 bytes per element) at 1 GB. Larger matrices are processed
#' in parallel in chunks of size at most \eqn{2^{26} = 6.7\cdot 10^7}{2^{26} = 6.7e+7}
#' elements. The number of threads is 4 by default, which can be changed by
#' \code{RcppParallel::setThreadOptions(numThreads = 8)} or something similar.
#'
#' @return A vector of density estimates evaluated at the grid points or, if \code{return.grid}, a matrix with the density in the last column.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(rt(10000, df = 5)) # Observed values
#' g <- seq(-6, 6, 0.05) # Grid for evaluation
#' d2 <- kernelDensity(x, g, bw = 0.3, kernel = "epanechnikov", no.dedup = TRUE)
#' d4 <- kernelDensity(x, g, bw = 0.4, kernel = "quartic", order = 4, no.dedup = TRUE)
#' plot(g, d2, ylim = range(0, d2, d4), type = "l"); lines(g, d4, col = 2)
#'
#' # De-duplication facilities for faster operations
#' set.seed(1)  # Creating a data set with many duplicates
#' n.uniq <- 1000
#' n <- 4000
#' inds <- ceiling(runif(n, 0, n.uniq))
#' x.uniq <- matrix(rnorm(n.uniq*10), ncol = 10)
#' x <- x.uniq[inds, ]
#' xout <- x.uniq[ceiling(runif(n.uniq*3, 0, n.uniq)), ]
#' w <- runif(n)
#' data.table::setDTthreads(1) # For measuring the pure gains and overhead
#' RcppParallel::setThreadOptions(numThreads = 1)
#' kd1 <- kernelDensity(x, xout, w, bw = 0.5)
#' kd2 <- kernelDensity(x, xout, w, bw = 0.5, no.dedup = TRUE)
#' stat1 <- attr(kd1, "duplicate.stats")
#' stat2 <- attr(kd2, "duplicate.stats")
#' print(stat1[3:5]) # De-duplication time -- worth it
#' print(stat2[3:5]) # Without de-duplication, slower
#' unname(prod((1 - stat1[1:2])) / (stat1[5] / stat2[5])) # > 1 = better time
#' # savings than expected, < 1 = worse time savings than expected
#' all.equal(as.numeric(kd1), as.numeric(kd2))
#' max(abs(kd1 - kd2)) # Should be around machine epsilon or less
kernelDensity <- function(x,
                          xout = NULL,
                          weights = NULL,
                          bw = NULL,
                          kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                          order = 2,
                          convolution = FALSE,
                          chunks = 0,
                          PIT = FALSE,
                          deduplicate.x = TRUE,
                          deduplicate.xout = TRUE,
                          no.dedup = FALSE,
                          return.grid = FALSE
) {
  arg <- prepareKernel(x = x, xout = xout, weights = weights, bw = bw,
                        kernel = kernel, PIT = PIT, order = order, convolution = convolution,
                        deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout, no.dedup = no.dedup)
  tic0 <- Sys.time()
  result <- kernelDensityCPP(x = arg$x, xout = arg$xout, weights = arg$weights, bw = arg$bw,
                             kernel = arg$kernel, order = arg$order, convolution = convolution,
                             chunks = as.integer(chunks))
  if (return.grid) result <- cbind(arg$xout, density = result)
  if (arg$deduplicate.xout) { # If de-duplication was done
    if (return.grid) result <- result[arg$xout.matches, ] else result <- result[arg$xout.matches]
  }
  ds <- c(arg$duplicate.stats, seconds.density = as.numeric(difftime(Sys.time(), tic0, units = "secs")))
  checkDupStat(ds)
  attr(result, "duplicate.stats") <- ds
  return(result)
}


#' Local kernel smoother
#'
#' @inheritParams kernelDensity
#' @param y A numeric vector of dependent variable values.
#' @param LOO Logical: If \code{TRUE}, the leave-one-out estimator is returned.
#' @param degree Integer: 0 for locally constant estimator (Nadaraya--Watson), 1 for
#'   locally linear (Cleveland's LOESS), 2 for locally quadratic (use with care, less stable, requires larger bandwidths)
#' @param trim Trimming function for small weights to speed up locally weighted regression (if \code{degree} is 1 or 2).
#' @param robust.iterations The number of robustifying iterations (due to Cleveland, 1979). If greater than 0, \code{xout} is ignored.
#' @param robust Character: "huber" for Huber's local regression weights, "bisquare" for more robust bi-square ones
#' @param return.grid If \code{TRUE}, prepends \code{xout} to the return results.
#'
#' Standardisation is recommended for the purposes of numerical stability (sometimes
#'   \code{lm()} might choke when the dependent variable takes very large absolute
#'   values and its square is used).
#'
#' The robust iterations are carried out, if requested, according to @cleveland1979robust.
#' Huber weights are never zero; bisquare weights create a more robust re-descending
#' estimator.
#'
#' Note: if \code{x} and \code{xout} are different but robust iterations were requested,
#'   the robustification can take longer. TODO: do not estimate on (x, grid),
#'   do the calculation with K.full straight away.
#'
#' Note: if \code{LOO} is used, it makes sense to de-duplicate observations first.
#'   By default, this behaviour is not enforced in this function, but when it is
#'   called in cross-validation routines, the de-duplication is forced. It makes
#'   no sense to zero out once observation out of many repeated.
#'
#' @return A vector of predicted values or, if \code{return.grid} is \code{TRUE},
#'   a matrix with the predicted values in the last column.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 300
#' x <- sort(rt(n, df = 6)) # Observed values
#' g <- seq(-4, 5, 0.1) # Grid for evaluation
#' f <- function(x) 1 + x + sin(x) # True E(Y | X) = f(X)
#' y <- f(x) + rt(n, df = 4)
#' # 3 estimators: locally constant + 2nd-order kernel,
#' # locally constant + 4th-order kernel, locally linear robust
#' b2lc <- suppressWarnings(bw.CV(x, y = y, kernel = "quartic")
#'                          + 0.8)
#' b4lc <- suppressWarnings(bw.CV(x, y = y, kernel = "quartic", order = 4,
#'               try.grid = FALSE, start.bw = 3) + 1)
#' b2ll <- bw.CV(x, y = y, kernel = "quartic", degree = 1, robust.iterations = 1,
#'               try.grid = FALSE, start.bw = 3, verbose = TRUE)
#' m2lc <- kernelSmooth(x, y, g, bw = b2lc, kernel = "quartic", no.dedup = TRUE)
#' m4lc <- kernelSmooth(x, y, g, bw = b4lc, kernel = "quartic", order = 4, no.dedup = TRUE)
#' m2ll <- kernelSmooth(x, y, g, bw = b2ll, kernel = "quartic",
#'                      degree = 1, robust.iterations = 1, no.dedup = TRUE)
#' plot(x, y, xlim = c(-6, 7), col = "#00000088", bty = "n")
#' lines(g, f(g), col = "white", lwd = 5); lines(g, f(g))
#' lines(g, m2lc, col = 2); lines(g, m4lc, col = 3); lines(g, m2ll, col = 4)

#' # De-duplication facilities for faster operations
#' set.seed(1)  # Creating a data set with many duplicates
#' n.uniq <- 1000
#' n <- 4000
#' inds <- sort(ceiling(runif(n, 0, n.uniq)))
#' x.uniq <- sort(rnorm(n.uniq))
#' y.uniq <- 1 + x.uniq + sin(x.uniq*2) + rnorm(n.uniq)
#' x <- x.uniq[inds]
#' y <- y.uniq[inds]
#' xout <- x.uniq[sort(ceiling(runif(n.uniq*3, 0, n.uniq)))]
#' w <- runif(n)
#' data.table::setDTthreads(1) # For measuring the pure gains and overhead
#' RcppParallel::setThreadOptions(numThreads = 1)
#' kr1 <- kernelSmooth(x, y, xout, w, bw = 0.2)
#' kr2 <- kernelSmooth(x, y, xout, w, bw = 0.5, no.dedup = TRUE)
#' stat1 <- attr(kr1, "duplicate.stats")
#' stat2 <- attr(kr2, "duplicate.stats")
#' print(stat1[3:5]) # De-duplication time -- worth it
#' print(stat2[3:5]) # Without de-duplication, slower
#' unname(prod((1 - stat1[1:2])) / (stat1[5] / stat2[5])) # > 1 = better time
#' # savings than expected, < 1 = worse time savings than expected
#' all.equal(as.numeric(kr1), as.numeric(kr2))
#' max(abs(kr1 - kr2)) # Should be around machine epsilon or less
#'
#' # Example in 2 dimensions
#' # TODO
kernelSmooth <- function(x, y, xout = NULL, weights = NULL,
                         bw = NULL,
                         kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                         order = 2,
                         convolution = FALSE,
                         chunks = 0,
                         PIT = FALSE, LOO = FALSE,
                         degree = 0,
                         trim = function(x) 0.01 / length(x),
                         robust.iterations = 0, robust = c("bisquare", "huber"),
                         deduplicate.x = TRUE, deduplicate.xout = TRUE,
                         no.dedup = FALSE, return.grid = FALSE
) {
  bw0 <- bw  # To be used later if there are robust iterations
  if (!(degree %in% 0:2)) stop("kernelSmooth: degree must be 0, 1, or 2.")
  robust <- robust[1]
  if (!(robust %in% c("huber", "bisquare"))) stop("kernelSmooth: 'robust' must be Huber (less robust) or 'bisquare' (redescending, robust).")
  if (LOO && !is.null(xout)) {
    warning("kernelSmooth: Leave-One-Out estimation requested, but a custom xout passed! Ignoring it.")
    xout <- NULL
  }
  arg <- prepareKernel(x = x, y = y, xout = xout, weights = weights, bw = bw, kernel = kernel,
                        PIT = PIT, order = order, convolution = convolution,
                        deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout, no.dedup = no.dedup)

  tic0 <- Sys.time()
  if (degree == 0) {
    result <- kernelSmoothCPP(x = arg$x, y = arg$y, xout = arg$xout, weights = arg$weights,
                              bw = arg$bw, kernel = arg$kernel, order = as.integer(arg$order),
                              LOO = as.logical(LOO), convolution = as.logical(convolution),
                              chunks = as.integer(chunks))
    # Define the bad rows later
  } else {
    x <- arg$x
    y <- arg$y
    xout <- arg$xout
    weights <- arg$weights
    bw <- arg$bw
    kernel <- arg$kernel
    order <- arg$order

    robW <- switch(robust, huber = function(x) ifelse(abs(x) < 1.345, 1, 1.345 / abs(x)),
                   bisquare = function(x) ifelse(abs(x) < 4.685, (1 - (x/4.685)^2)^2, 0))

    is.exact <- isTRUE(all.equal(x, xout))
    # No PIT here because arg$x is already transformed
    K <- kernelWeights(x = x, xout = xout, bw = bw, kernel = kernel, order = order, convolution = convolution)
    K <- sweep(K, 2, weights, "*")
    if (LOO) diag(K) <- 0
    K <- K / rowSums(K)
    w.list <- apply(K, 1, sparseVectorToList, trim = trim)
    needs.full <- (!is.exact) & robust.iterations > 0
    if (needs.full) { # for final weights on x & xout because they are not equal
      # No PIT here because arg$x is already transformed
      # However, we require a bandwidth matrix for X, not xout, hence we can use only the initial v
      arg.full <- prepareKernel(x = x, y = y, xout = x, weights = weights, bw = bw0, kernel = kernel,
                           order = order, convolution = convolution,
                           deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout, no.dedup = no.dedup)
      K.full <- kernelWeights(x = arg.full$x, bw = arg.full$bw, kernel = arg.full$kernel, order = arg.full$order, convolution = convolution)
      K.full <- sweep(K.full, 2, weights, "*")
      if (LOO) diag(K.full) <- 0
      K.full <- K.full / rowSums(K.full)
      w.full.list <- apply(K.full, 1, sparseVectorToList, trim = trim)
    }

    m <- apply(x, 2, stats::median) # Standardising for LM fit stability
    s <- apply(x, 2, stats::IQR)
    xs <- sweep(sweep(x, 2, m, "-"), 2, s, "/")
    xouts <- sweep(sweep(xout, 2, m, "-"), 2, s, "/")
    if (degree == 2) {
      xs <- cbind(xs, xs^2)
      xouts <- cbind(xouts, xouts^2)
    }

    WOLS <- function(nonzw, robustw = NULL) {
      dimb <- ncol(xs)+1
      # If there are no non-zero neighbours or one point has too much influence, declare failure
      # Zero indices can happen in LOO cross-validation
      if (length(nonzw$idx) <= 1 || any(nonzw$ct > 0.999)) return(rep(NA, dimb))
      wreg <- sqrt(nonzw$ct)
      if (!is.null(robustw)) wreg <- wreg * sqrt(robustw[nonzw$idx])
      xw <- cbind(1, xs[nonzw$idx, , drop = FALSE]) * wreg
      yw <- y[nonzw$idx] * wreg
      # If there are fewer valid observations than necessary for a fit, return NA
      if (nrow(xw) < dimb) return(rep(NA, dimb))
      # If for any other reason the fit fails, safeguard
      b <- tryCatch(stats::.lm.fit(xw, yw)$coefficients, error = function(e) return(rep(NA, dimb)))
      return(b)
    }

    findBadRows <- function(coefhat) {
      bad.rows <- which(apply(coefhat, 1, function(x) any(is.na(x))))
      msg <- c("Local weighted fit numerically unstable: one point has a huge kernel weight > 0.999.",
               "Potentially no neighbours with weights > trimming value, terminating to avoid a singular fit.",
               paste0("Problematic 'xout' row indices: ", paste0(bad.rows, collapse = ", "), ")"),
               "Try increasing the bandwidth to get more data points in the vicinity.")
      if (length(bad.rows) > 0) warning(paste0(msg, collapse = "\n"))
      return(bad.rows)
    }

    if (!needs.full) {
      # First stage: if no robust iterations were requested, just do the job
      # If robust iterations were requested and x = xout, do the first iteration
      coefhat <- do.call(rbind, lapply(seq_along(w.list), function(i) WOLS(w.list[[i]])))
      result <- rowSums(coefhat * cbind(1, xouts))
      bad <- findBadRows(coefhat) # Print warnings
    } else {
      # Otherwise, when x != xout; a local regression for each point in x
      # must be estimated first
      deltak <- rep(1, length(w.full.list))
      for (i in 1:robust.iterations) {
        coefhat.full <- do.call(rbind, lapply(seq_along(w.full.list), function(i) WOLS(w.full.list[[i]], deltak)))
        result.full <- rowSums(coefhat.full * cbind(1, xs))
        bad <- findBadRows(coefhat.full)
        resid <- y - result.full
        ss <- stats::median(abs(resid), na.rm = TRUE)
        deltak <- robW(resid / ss / 6) # Robust LOESS weights from Cleveland 1979
      }
      # Now that robust observation weights have been computed, do OLS on a grid
      coefhat <- do.call(rbind, lapply(seq_len(length(w.list)), function(i) WOLS(w.list[[i]], deltak)))
      result  <- rowSums(coefhat * cbind(1, xouts))
    }
  }

  # De-duplication summary
  ds <- c(arg$duplicate.stats, seconds.smoother = as.numeric(difftime(Sys.time(), tic0, units = "secs")))
  checkDupStat(ds)

  if (return.grid) result <- cbind(xout, mu = result)
  if (arg$deduplicate.xout) { # If de-duplication was done
    result <- if (return.grid) result[arg$xout.matches, ] else result[arg$xout.matches]
  }
  if (degree == 0) bad <- which(!is.finite(result)) # This belongs after the restoration of duplicates
  if (any(is.nan(result)))
    warning(paste0("Some smoothed values are NaN, which happens (among other reasons) when ",
                   "the sum of smoothing weights is exactly zero (bandwidth ",
                   paste0(sprintf("%1.2f", arg$bw), collapse = "; "), " -- try increasing it) ",
                   "or when there are NA's in the input 'x' data (address the missing values)."))
  if ((!any(is.nan(result))) && any(!is.finite(result)))
    warning(paste0("Some smoothed values are NA or Inf, which is really strange. ",
                   "Possible reasons: NA's in the input 'x' data (address the missing values) ",
                   "or something else. Check the inputs via 'all(is.finite(x))' ",
                   "or run 'debugonce(kernelSmooth)' and retry the last call."))

  attr(result, "duplicate.stats") <- ds
  attr(result, "bad.indices") <- bad
  return(result)
}

#' Density and/or kernel regression estimator with conditioning on discrete variables
#'
#' @param x A vector or a matrix/data frame of discrete explanatory variables (exogenous).
#'   Non-integer values are fine because the data are split into bins defined by interactions of these variables.
#' @param y Optional: a vector of dependent variable values.
#' @param compact Logical: return unique values instead of full data with repeated observations?
#' @param fun A function that computes a statistic of \code{y} inside every category defined by \code{x}.
#'
#' @return A list with \code{x}, density estimator (\code{fhat}) and, if \code{y} was provided, regression estimate.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(rnorm(1000))
#' p <- 0.5*pnorm(x) + 0.25 # Propensity score
#' d <- as.numeric(runif(1000) < p)
#' # g = discrete version of x for binning
#' g <- as.numeric(as.character(cut(x, -4:4, labels = -4:3+0.5)))
#' dhat.x <- kernelSmooth(x = x, y = d, bw = 0.4, no.dedup = TRUE)
#' dhat.g <- kernelDiscreteDensitySmooth(x = g, y = d)
#' dhat.comp <- kernelDiscreteDensitySmooth(g, d, compact = TRUE)
#' plot(x, p, ylim = c(0, 1), bty = "n", type = "l", lty = 2)
#' points(x, dhat.x, col = "#00000044")
#' points(dhat.comp, col = 2, pch = 16, cex = 2)
#' lines(dhat.comp$x, dhat.comp$fhat, col = 4, pch = 16, lty = 3)
kernelDiscreteDensitySmooth <- function(x,
                                        y = NULL,
                                        compact = FALSE,
                                        fun = mean
) {
  if (is.matrix(x)) x <- as.data.frame(x)
  if (is.data.frame(x)) x <- as.integer(interaction(x, drop = TRUE, lex.order = TRUE))
  if (!is.vector(x)) stop("x must be a numeric vector, matrix, or a data frame!")
  n <- NROW(x)
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
#' @inheritParams kernelDensity
#' @param by A variable containing unique identifiers of discrete categories.
#' @param byout A variable containing unique identifiers of discrete categories
#'   for the output grid (same points as \code{xout})
#' @param parallel Logical: if \code{TRUE}, parallelises the calculation over
#'   the unique values of \code{by}. At this moment, supports only
#'   \code{parallel::mclapply} (therefore, will not work on Windows).
#' @param cores Integer: the number of CPU cores to use. High core count = high RAM usage.
#'   If the number of unique values of 'by' is less than the number of cores requested,
#'   then, only \code{length(unique(by))} cores are used.
#' @param preschedule Logical: passed as \code{mc.preschedule} to \code{mclapply}.
#' @param ... Passed to \code{kernelDensity}.
#'
#' @return A numeric vector of the density estimate of the same length as \code{nrow(xout)}.
#' @export
#'
#' @examples
#' # Estimating 3 densities on something like a panel
#' set.seed(1)
#' n <- 200
#' x <- c(rnorm(n), rchisq(n, 4)/4, rexp(n, 1))
#' by <- rep(1:3, each = n)
#' xgrid <- seq(-3, 6, 0.1)
#' out <- expand.grid(x = xgrid, by = 1:3)
#' fhat <- kernelMixedDensity(x = x, xout = out$x, by = by, byout = out$by)
#' plot(xgrid, dnorm(xgrid)/3, type = "l", bty = "n", lty = 2, ylim = c(0, 0.35),
#'      xlab = "", ylab = "Density")
#' lines(xgrid, dchisq(xgrid*4, 4)*4/3, lty = 2, col = 2)
#' lines(xgrid, dexp(xgrid, 1)/3, lty = 2, col = 3)
#' for (i in 1:3) {
#'   lines(xgrid, fhat[out$by == i], col = i, lwd = 2)
#'   rug(x[by == i], col = i)
#' }
#' legend("top", c("00", "10", "01", "11"), col = 2:5, lwd  = 2)
kernelMixedDensity <- function(x, by, xout = NULL, byout = NULL, weights = NULL, parallel = FALSE, cores = 1, preschedule = TRUE, ...) {
  .kernelMixed(x = x, by = by, xout = xout, byout = byout, weights = weights, type = "density", parallel = parallel, cores = cores, preschedule = preschedule, ...)
}


#' Smoothing with conditioning on discrete and continuous variables
#'
#' @inheritParams kernelSmooth
#' @inheritParams kernelMixedDensity
#' @param ... Passed to \code{kernelSmooth} (usually \code{bw}, \code{gaussian} for both; \code{degree} and \code{robust.iterations} for "smooth"),
#'
#' @return A numeric vector of the kernel estimate of the same length as \code{nrow(xout)}.
#' @export
#'
#' @examples
#' set.seed(1)
#' n <- 1000
#' z1 <- rbinom(n, 1, 0.5)
#' z2 <- rbinom(n, 1, 0.5)
#' x <- rnorm(n)
#' u <- rnorm(n)
#' y <- 1 + x^2 + z1 + 2*z2 + z1*z2 + u
#' by <- as.integer(interaction(list(z1, z2)))
#' out <- expand.grid(x = seq(-4, 4, 0.25), by = 1:4)
#' yhat <- kernelMixedSmooth(x = x, y = y, by = by, bw = 1, degree = 1,
#'                           xout = out$x, byout = out$by)
#' plot(x, y)
#' for (i in 1:4) lines(out$x[out$by == i], yhat[out$by == i], col = i+1, lwd = 2)
#' legend("top", c("00", "10", "01", "11"), col = 2:5, lwd  = 2)
#'
#' # The function works faster if there are duplicated values of the
#' # conditioning variables in the prediction grid and there are many
#' # observations; this is illustrated by the following example
#' # without a custom grid
#' # In this example, ignore the fact that the conditioning variable is rounded
#' # and therefore contains measurement error (ruining consistency)
#' x  <- rnorm(10000)
#' xout <- rnorm(5000)
#' xr <- round(x)
#' xrout <- round(xout)
#' w <- runif(10000, 1, 3)
#' y  <- 1 + x^2 + rnorm(10000)
#' by <- rep(1:4, each = 2500)
#' byout <- rep(1:4, each = 1250)
#' system.time(kernelMixedSmooth(x = x, y = y, by = by, weights = w,
#'                               xout = xout, byout = byout, bw = 1))
#' #  user  system elapsed
#' # 0.144   0.000   0.144
#' system.time(km1 <- kernelMixedSmooth(x = xr, y = y, by = by, weights = w,
#'                                      xout = xrout, byout = byout, bw = 1))
#' #  user  system elapsed
#' # 0.021   0.000   0.022
#' system.time(km2 <- kernelMixedSmooth(x = xr, y = y, by = by, weights = w,
#'                      xout = xrout, byout = byout, bw = 1, no.dedup = TRUE))
#' #  user  system elapsed
#' # 0.138   0.001   0.137
#' all.equal(km1, km2)
#'
#' # Parallel capabilities shine in large data sets
#' if (.Platform$OS.type != "windows") {
#' # A function to carry out the same estimation in multiple cores
#' pFun <- function(n) kernelMixedSmooth(x = rep(x, 2), y = rep(y, 2),
#'          weights = rep(w, 2), by = rep(by, 2),
#'          bw = 1, degree = 0, parallel = TRUE, cores = n)
#' system.time(pFun(1))  # 0.6--0.7 s
#' system.time(pFun(2))  # 0.4--0.5 s
#' }
kernelMixedSmooth <- function(x, y, by, xout = NULL, byout = NULL, weights = NULL, parallel = FALSE, cores = 1, preschedule = TRUE, ...) {
  .kernelMixed(x = x, y = y, by = by, xout = xout, byout = byout, weights = weights, type = "smooth", parallel = parallel, cores = cores, preschedule = preschedule, ...)
}

.kernelMixed <- function(x, y = NULL, by, xout = NULL, byout = NULL,
                         weights = NULL, PIT = FALSE,
                         type = c("density", "smooth"),
                         deduplicate.x = TRUE, deduplicate.xout = TRUE, no.dedup = FALSE,
                         parallel = FALSE, cores = 1L, preschedule = TRUE,
                         ...
) {
  type <- type[1]
  dot.args <- list(...)
  if (is.null(y) && type == "smooth") stop("Supply the mandatory 'y' argument to obtain a kernel regression smoother.")
  if (any(by != round(by))) stop("'by' must be an integer (consider using 'as.integer').")
  if (length(by) != NROW(x)) stop("The length of 'byout' must be the same as NROW(xout) (because they correspond to the same observations).")
  if (is.null(xout)) xout <- x # This is not a redundancy; prepareKernel must receive full data for de-duplication
  if (is.null(byout)) byout <- by
  if (length(byout) != NROW(xout)) stop("The length of 'byout' must be the same as NROW(xout) (because they correspond to the same observations).")
  if (!all(unique(byout) %in% unique(by))) stop("The unique 'byout' prediction categories have new values not present in the input training data.")
  by2 <- as.integer(factor(c(by, byout)))
  by <- by2[seq_along(by)]
  byout <- by2[(length(by)+1):length(by2)]

  if (any(table(by) < 2)) warning("Some categories have only 1 observation; the distribution is degenerate. At least 2 obs. per category are needed.")

  arg <- prepareKernel(x = cbind(x, by), y = y, xout = cbind(xout, byout),
                        weights = weights, bw = 1, PIT = PIT,
                        deduplicate.x = deduplicate.x, deduplicate.xout = deduplicate.xout, no.dedup = no.dedup)
  arg$by <- as.integer(arg$x[, ncol(arg$x)])
  arg$byout <- if (!isTRUE(dot.args$LOO)) as.integer(arg$xout[, ncol(arg$xout)]) else arg$by
  arg$x <- arg$x[, -ncol(arg$x), drop = FALSE]
  arg$xout <- if (!isTRUE(dot.args$LOO)) arg$xout[, -ncol(arg$xout), drop = FALSE] else arg$x
  if (isTRUE(dot.args$LOO) && arg$deduplicate.x) {
    arg$deduplicate.xout <- TRUE
    arg$xout.matches <- arg$x.matches
    arg$duplicate.stats["dup.rate.xout"] <- arg$duplicate.stats["dup.rate.x"]
    arg$duplicate.stats["seconds.xout"] <- 0
  }
  # This arg$xout be ignored later, but is necessary to create the de-duplicated output with the same number of observations
  n <- sum(arg$weights)

  res <- rep(NA_real_, nrow(arg$xout))
  k <- max(arg$by) # Number of partitions
  s.list <- lapply(1L:k, function(i) arg$by == i)
  sout.list <- lapply(1L:k, function(i) arg$byout == i)
  innerFun <- function(i) {
    s <- s.list[[i]]
    sout <- sout.list[[i]]
    x.sub <- arg$x[s, , drop = FALSE]
    w.sub <- arg$weights[s]
    y.sub <- arg$y[s]
    xout.sub <- arg$xout[sout, , drop = FALSE]
    if (NROW(xout.sub) < 1) return(NULL) # If the new data do not have these values, skip
    res.sub <- switch(type,
                      density = kernelDensity(x = x.sub,            weights = w.sub, xout = xout.sub, no.dedup = TRUE, ...) * (sum(w.sub) / n),
                      smooth  = kernelSmooth(x  = x.sub, y = y.sub, weights = w.sub, xout = if (!isTRUE(dot.args$LOO)) xout.sub else NULL, no.dedup = TRUE, ...))
    # Initially, xout is set to x in case LOO is requested; this produces de-duplicated xout and byout of correct size
    # NULL is passed to kernelSmooth in the case of LOO CV to avoid warnings (does not affect the result)
    return(res.sub)
  }
  if (cores > k) cores <- k
  res.list <- if (parallel && cores > 1) parallel::mclapply(X = 1L:k, FUN = innerFun, mc.preschedule = preschedule, mc.cores = cores) else lapply(1L:k, innerFun)
  for (i in 1L:k) res[sout.list[[i]]] <- res.list[[i]] # Same order as before

  if (arg$deduplicate.xout) res <- res[arg$xout.matches]
  return(res)
}

#' Density cross-validation
#'
#' @inheritParams kernelDensity
#' @inheritParams kernelWeights
#' @param bw Candidate bandwidth values: scalar, vector, or a matrix (with columns corresponding to columns of \code{x}).
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#'
#' Note: since DCV requires computing the leave-one-out estimator,
#'   repeated observations are combined first; the de-duplication is therefore
#'   forced in cross-validation. The only situation where de-duplication can be
#'   skipped is passing de-duplicated data sets from outside (e.g. inside
#'   optimisers).
#'
#' @return A numeric vector of the same length as \code{bw} or \code{nrow(bw)}.
#'
#' @export
#' @examples
#' set.seed(1)
#' x <- rlnorm(100); x <- c(x[1], x)  # x with 1 duplicate
#' bws <- exp(seq(-3, 0.5, 0.1))
#' plot(bws, DCV(x, bws), log = "x", bty = "n", main = "Density CV")
DCV <- function(x, bw, weights = NULL, same = FALSE, kernel = "gaussian", order = 2,
                PIT = FALSE, chunks = 0, no.dedup = FALSE) {
  arg <- prepareKernel(x = x, weights = weights, bw = bw[1],
                        # bw[1] skips the length check in case multiple bw values are given as a vector not equal to NCOL(x)
                        kernel = kernel, PIT = PIT, order = order,
                        convolution = FALSE, deduplicate.x = !no.dedup)
  one.dim <- ncol(arg$x) == 1 # Are our data one-dimensional?
  n <- sum(arg$weights)

  if (one.dim) {
    many.bw <- (length(bw) > 1)
  } else {
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    # If the bandwidth is a matrix, each row = bandwidth; if the same bandwidth is enforces, then, a vector with > 1 element = many
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else {
      bw <- list(bw) # If there is only one bw, make it a list
    }
  }
  CV <- function(b) {
    # A sub-function to compute the CV for one BW, parallelisable
    if (any(b <= 0)) return(Inf)
    if (!one.dim && length(b) == 1) b <- rep(b, ncol(x))
    # No PIT here because arg$x is already transformed
    # Term 1: int ^f(x)^2 dx
    KK <- kernelWeights(x = arg$x, bw = b, kernel = arg$kernel, order = arg$order, convolution = TRUE, no.dedup = TRUE)
    KK <- sweep(KK, 2, arg$weights, "*")
    pb <- prod(b)
    term1 <- sum(KK) / (n^2 * pb)
    # Computing the LOO estimator efficiently:
    # fhat_i(X[i]) = 1 / (n-1) / prod(b) sum_{j != i} K(X[j] - X[i], b)
    # (n-1) fhat_i(X[i]) = 1 / prod(b) [-K(0) + sum_{j} K(X[j] - X[i], b)]
    # (n-1) fhat_i(X[i]) = n fhat(X[i]) - K(0) / prod(b)
    # fhat_i(X[i]) = n/(n-1) fhat(X[i]) - K(0) / prod(b) / (n-1)
    # n-1 gets replaced with n - w[i] in case of weights
    K0   <- as.numeric(kernelWeights(x = matrix(0, ncol = length(b)), bw = b, kernel = arg$kernel, order = arg$order, no.dedup = TRUE))
    fhat <- kernelDensity(x = arg$x, weights = arg$weights, bw = b, kernel = arg$kernel, order = arg$order, chunks = chunks, no.dedup = TRUE)
    fhat.LOO <- (n*fhat - K0/pb) / (n - arg$weights)
    term2 <- -2 * mean(fhat.LOO)
    return(term1 + term2)
  }
  CV.values <- vapply(bw, CV, FUN.VALUE = numeric(1))
  attr(CV.values, "duplicate.stats") <- arg$duplicate.stats
  return(CV.values)
}

#' Least-squares cross-validation function for the Nadaraya-Watson estimator
#'
#' @inheritParams kernelSmooth
#' @param bw Candidate bandwidth values: scalar, vector, or a matrix (with columns corresponding to columns of \code{x}).
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#' @param cores Integer: the number of CPU cores to use. High core count = high RAM usage.
#'
#' Note: since LSCV requires zeroing out the diagonals of the weight matrix,
#'   repeated observations are combined first; the de-duplication is therefore
#'   forced in cross-validation. The only situation where de-duplication can be
#'   skipped is passing de-duplicated data sets from outside (e.g. inside
#'   optimisers).
#'
#' @return A numeric vector of the same length as \code{bw} or \code{nrow(bw)}.
#' @export
#'
#' @examples
#' set.seed(1)  # Creating a data set with many duplicates
#' n.uniq <- 1000
#' n <- 4000
#' inds <- sort(ceiling(runif(n, 0, n.uniq)))
#' x.uniq <- sort(rnorm(n.uniq))
#' y.uniq <- 1 + 0.2*x.uniq + 0.3*sin(x.uniq) + rnorm(n.uniq)
#' x <- x.uniq[inds]
#' y <- y.uniq[inds]
#' w <- 1 + runif(n, 0, 2) # Relative importance
#' data.table::setDTthreads(1) # For measuring pure gains and overhead
#' RcppParallel::setThreadOptions(numThreads = 1)
#' bw.grid <- seq(0.1, 1.2, 0.05)
#' ncores <- if (.Platform$OS.type == "windows") 1 else 2
#' CV <- LSCV(x, y, bw.grid, weights = w, cores = ncores)  # Parallel capabilities
#' bw.opt <- bw.grid[which.min(CV)]
#' g <- seq(-3.5, 3.5, 0.05)
#' yhat <- kernelSmooth(x, y, xout = g, weights = w,
#'                      bw = bw.opt, deduplicate.xout = FALSE)
#' oldpar <- par(mfrow = c(2, 1), mar = c(2, 2, 2, 0)+.1)
#' plot(bw.grid, CV, bty = "n", xlab = "", ylab = "", main = "Cross-validation")
#' plot(x.uniq, y.uniq, bty = "n", xlab = "", ylab = "", main = "Optimal fit")
#' points(g, yhat, pch = 16, col = 2, cex = 0.5)
#' par(oldpar)
LSCV <- function(x, y, bw, weights = NULL, same = FALSE, degree = 0, kernel = "gaussian",
                 order = 2, PIT = FALSE, chunks = 0, robust.iterations = 0, cores = 1) {
  arg <- prepareKernel(x = x, y = y, weights = weights, bw = bw[1],
                       # bw[1] skips the length check in case multiple bw values are given as a vector not equal to NCOL(x)
                       kernel = kernel, PIT = PIT, order = order,
                       convolution = FALSE, deduplicate.x = TRUE)
  one.dim <- ncol(arg$x) == 1 # Are the data one-dimensional? Determines how vector bw is handled.
  n <- sum(arg$weights)

  if (is.data.frame(bw)) bw <- as.matrix(bw)
  if (one.dim) {
    many.bw <- (length(bw) > 1)
  } else {
    many.bw <- (!is.null(dim(bw))) | (is.vector(bw) & length(bw) > 1 & same)
    # If the bandwidth is a matrix, each row = bandwidth; if the same bandwidth is enforces, then, a vector with > 1 element = many
    if (many.bw) {
      if (!is.vector(bw)) bw <- lapply(seq_len(nrow(bw)), function(i) bw[i, ]) # If the input is a matrix, split it into a list
      # Otherwise, [mc]lapply will happily eat a vector
    } else {
      bw <- list(bw) # If there is only one bw, make it a list
    }
  }
  ASE <- function(b) { # Accepts the already-deduplicated data
    if (any(b <= 0)) return(Inf)
    muhat.i <- kernelSmooth(x = arg$x, y = arg$y, bw = b, weights = arg$weights, kernel = arg$kernel,
                            order = arg$order, degree = degree, LOO = TRUE, chunks = chunks,
                            robust.iterations = robust.iterations, no.dedup = TRUE)
    m <- sum((arg$y - muhat.i)^2 * arg$weights) / n
    if (!is.finite(m)) m <- Inf
    return(m)
  }
  ASE.values <- unlist(parallel::mclapply(bw, ASE, mc.cores = cores))
  attr(ASE.values, "duplicate.stats") <- arg$duplicate.stats
  return(ASE.values)
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
#' If \code{y} is NULL and only \code{x} is supplied, returns the density-cross-validated bandwidth (DCV).
#' If \code{y} is supplied, then, returns the least-squares-cross-validated bandwidth (LSCV).
#'
#' @inheritParams kernelDensity
#' @inheritParams kernelSmooth
#' @inheritParams kernelWeights
#' @inheritParams prepareKernel
#' @param y A numeric vector of responses (dependent variable) if the user wants least-squares cross-validation.
#' @param robust.iterations Passed to \code{kernelSmooth} if \code{y} is not \code{NULL} (for least-squares CV).
#' @param degree Passed to \code{kernelSmooth} if \code{y} is not \code{NULL} (for least-squares CV).
#' @param start.bw Numeric vector: initial value for bandwidth search.
#' @param same Logical: use the same bandwidth for all columns of \code{x}?
#' @param tol Relative tolerance used by the optimiser as the stopping criterion.
#' @param try.grid Logical: if true, 10 different bandwidths around the rule-of-thumb
#'   one are tried with multiplier \code{1.2^(-3:6)}
#' @param ndeps Numerical-difference epsilon. Puts a lower bound on the result: the estimated optimal bw
#'   cannot be less than this value.
#' @param verbose Logical: print out the optimiser return code for diagnostics?
#' @param attach.attributes Logical: if TRUE, returns the output of `optim()` for diagnostics.
#' @param control List: extra arguments to pass to the control-argument list of `optim`.
#'
#' @return Numeric vector or scalar of the optimal bandwidth.
#' @export
#'
#' @examples
#' set.seed(1)  # Creating a data set with many duplicates
#' n.uniq <- 200
#' n <- 500
#' inds <- sort(ceiling(runif(n, 0, n.uniq)))
#' x.uniq <- sort(rnorm(n.uniq))
#' y.uniq <- 1 + 0.1*x.uniq + sin(x.uniq) + rnorm(n.uniq)
#' x <- x.uniq[inds]
#' y <- y.uniq[inds]
#' w <- 1 + runif(n, 0, 2) # Relative importance
#' data.table::setDTthreads(1) # For measuring the pure gains and overhead
#' RcppParallel::setThreadOptions(numThreads = 1)
#' bw.grid <- seq(0.1, 1.3, 0.2)
#' CV <- LSCV(x, y, bw.grid, weights = w)
#' bw.init <- bw.grid[which.min(CV)]
#' bw.opt <- bw.CV(x, y, w) # 0.49, very close
#' g <- seq(-3.5, 3.5, 0.05)
#' yhat <- kernelSmooth(x, y, g, w, bw.opt, deduplicate.xout = FALSE)
#' oldpar <- par(mfrow = c(2, 1), mar = c(2, 2, 2, 0)+.1)
#' plot(bw.grid, CV, bty = "n", xlab = "", ylab = "", main = "Cross-validation")
#' points(bw.opt, LSCV(x, y, bw.opt, w), col = 2, pch = 15)
#' plot(x.uniq, y.uniq, bty = "n", xlab = "", ylab = "", main = "Optimal fit")
#' points(g, yhat, pch = 16, col = 2, cex = 0.5)
#' par(oldpar)
bw.CV <- function(x, y = NULL, weights = NULL,
                  kernel = "gaussian", order = 2, PIT = FALSE,
                  chunks = 0,
                  robust.iterations = 0, degree = 0,
                  start.bw = NULL, same = FALSE,
                  tol = 1e-4, try.grid = TRUE, ndeps = 1e-5,
                  verbose = FALSE, attach.attributes = FALSE,
                  control = list()) {
  CV <- if (!is.null(y)) "LSCV" else "DCV"
  arg <- prepareKernel(x = x, y = y, weights = weights, bw = 1,
                        # The value '1' skips the bw length check
                        kernel = kernel, PIT = PIT, order = order,
                        convolution = FALSE, deduplicate.x = TRUE)
  arg.list <- list(x = arg$x, weights = arg$weights, kernel = arg$kernel, order = arg$order,
                   chunks = chunks, same = same) # Processed data to pass further
  if (CV == "LSCV") arg.list <- c(arg.list, list(y = arg$y), degree = degree, robust.iterations = robust.iterations)

  f.to.min <- function(b) {
    if (verbose) print(b)
    if (isTRUE(any(b <= ndeps))) return(Inf)
    ret <- if (CV == "DCV") do.call(DCV, c(arg.list, list(b = b))) else do.call(LSCV, c(arg.list, list(b = b)))
    return(ret)
  }

  if (is.null(start.bw)) start.bw <- bw.rot(arg$x, kernel = arg$kernel, na.rm = TRUE) * 1.5
  xgaps <- apply(arg$x, 2, function(a) max(diff(sort(a))))
  if (verbose) cat("Rule-of-thumb bw: (", paste0(start.bw, collapse = ", "), "); max. gap between obs.: (",
                   paste0(xgaps, collapse = ", "), "). Taking the maximum as the starting point.\n", sep = "")
  start.bw <- pmax(start.bw, xgaps * (1 + 1e-3))
  if (same && ncol(arg$x) > 1) start.bw <- stats::quantile(start.bw, 0.75) # To the over-smoothing side

  optim.control <- list(reltol = tol, REPORT = 1, trace = if (verbose) 2 else 0, ndeps = rep(ndeps, length(start.bw)))
  start.bw <- unname(start.bw) # For proper handling of named arguments inside do.call
  f0 <- suppressWarnings(f.to.min(start.bw))
  # A quick grid search to avoid multi-modality
  if (try.grid) {
    bgrid <- lapply(-3:6, function(p) start.bw * 1.25^p)
    CV.grid <- suppressWarnings(sapply(bgrid, f.to.min))
    if (any(is.finite(CV.grid))) {
      start.bw <-  bgrid[[which.min(CV.grid)]]
      f0 <- CV.grid[which.min(CV.grid)]
    } else {
      start.bw <- bgrid[[10]] # The largest one
      f0 <- CV.grid[[10]]
    }
  }
  if (verbose) cat("After a quick grid evaluation, starting bandwidth search at", start.bw, "...\n")
  if (!is.finite(f0)) {
    for (i in 1:25) {
      if (verbose) cat("Bandwidth (", paste0(start.bw, collapse = ", "), ") too small, increasing by 20%...\n", sep = "")
      start.bw <- start.bw * 1.25
      f0 <- suppressWarnings(f.to.min(start.bw))
      if (is.finite(f0)) break
    }
  }
  opt.result <- tryCatch(stats::optim(par = start.bw, fn = f.to.min, method = "BFGS", control = optim.control),
                         error = function(e) return(e))
  if (inherits(opt.result, "error")) {
    if (grepl("non-finite finite-difference", opt.result$message)) {
      warning("CV gradient could not be computed (initial bw at the boundary? the optimiser tries strange values?). Retrying optimisation with the Nelder-Mead method.")
    } else {
      warning("Generic optimiser error. Retrying optimisation with the Nelder-Mead method.")
    }
    opt.result <- tryCatch(stats::optim(par = start.bw, fn = f.to.min, method = "Nelder-Mead",
                                        control = list(reltol = tol, REPORT = 1, trace = if (verbose) 2 else 0), warn.1d.NelderMead = FALSE),
                           error = function(e) return(e))
  }
  if (inherits(opt.result, "error")) {
    warning(paste0("'optim' failed to optimise the bandwidth. Returning a very rough rule-of-thumb value. Reason: ",
                   as.character(opt.result)))
    return(start.bw)
  }

  if (verbose) message(paste0("optim exit code ", opt.result$convergence, ", done in (", paste0(opt.result$counts, collapse = ", "), ") iterations."))
  bw <- opt.result$par
  if (attach.attributes) attr(bw, "optim") <- opt.result[names(opt.result) != "par"]
  return(bw)
}


#' Basic univatiate kernel functions
#'
#' Computes 5 most popular kernel functions of orders 2 and 4 with the potential of returning
#' an analytical convolution kernel for density cross-validation. These kernels appear
#' in \insertCite{silverman1986density}{smoothemplik}.
#'
#' @param x A numeric vector of values at which to compute the kernel function.
#' @param kernel Kernel type: uniform, Epanechnikov, triangular, quartic, or Gaussian.
#' @param order Kernel order. 2nd-order kernels are always non-negative.
#'   *k*-th-order kernels have all moments from 1 to (k-1) equal to zero, which is
#'   achieved by having some negative values.
#' \eqn{\int_{-\infty}^{+\infty} x^2 k(x) = \sigma^2_k = 1}.
#' This is useful because in this case, the constant \code{k_2} in formulæ 3.12 and 3.21
#' from \insertCite{silverman1986density;textual}{smoothemplik} is equal to 1.
#' @param convolution Logical: return the convolution kernel? (Useful for density cross-validation.)
#'
#' @details
#' The kernel functions take non-zero values on \eqn{[-1, 1]}{[-1, 1]}, except for
#' the Gaussian one, which is supposed to have full support, but due to the rapid
#' decay, is indistinguishable from machine epsilon outside
#' \eqn{[-8.2924, 8.2924]}{[-8.2924, 8.2924]}.
#'
#'
#' @return A numeric vector of the same length as input.
#' @importFrom Rdpack reprompt
#' @export
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ks <- c("uniform", "triangular", "epanechnikov", "quartic", "gaussian"); names(ks) <- ks
#' os <- c(2, 4); names(os) <- paste0("o", os)
#' cols <- c("#000000CC", "#0000CCCC", "#CC0000CC", "#00AA00CC", "#BB8800CC")
#' put.legend <- function() legend("topright", legend = ks, lty = 1, col = cols, bty = "n")
#' xout <- seq(-4, 4, length.out = 301)
#' plot(NULL, NULL, xlim = range(xout), ylim = c(0, 1.1),
#'   xlab = "", ylab = "", main = "Unscaled kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xout, kernelFun(xout, kernel = ks[i]), col = cols[i])
#' oldpar <- par(mfrow = c(1, 2))
#' plot(NULL, NULL, xlim = range(xout), ylim = c(-0.1, 0.8), xlab = "", ylab = "",
#'   main = "4th-order kernels", bty = "n"); put.legend()
#' for (i in 1:5) lines(xout, kernelFun(xout, kernel = ks[i], order = 4), col = cols[i])
#' par(mfrow = c(1, 1))
#' plot(NULL, NULL, xlim = range(xout), ylim = c(-0.25, 1.4), xlab = "", ylab = "",
#'   main = "Convolution kernels", bty = "n"); put.legend()
#' for (i in 1:5) {
#'   for (j in 1:2) lines(xout, kernelFun(xout, kernel = ks[i], order = os[j],
#'   convolution = TRUE), col = cols[i], lty = j)
#' }; legend("topleft", c("2nd order", "4th order"), lty = 1:2, bty = "n")
#' par(oldpar)
#'
#' # All kernels integrate to correct values; we compute the moments
#' mom <- Vectorize(function(k, o, m, c) integrate(function(x) x^m * kernelFun(x, k, o,
#'   convolution = c), lower = -Inf, upper = Inf)$value)
#' for (m in 0:4) {
#'   cat("\nComputing integrals of x^", m, " * f(x). \nSimple unscaled kernel:\n", sep = "")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = FALSE)), 4))
#'   cat("Convolution kernel:\n")
#'   print(round(outer(os, ks, function(o, k) mom(k, o, m = m, c = TRUE)), 4))
#' }
#'
kernelFun <- function(x,
                      kernel = c("gaussian", "uniform", "triangular", "epanechnikov", "quartic"),
                      order = c(2, 4),
                      convolution = FALSE
) {
  is.arr <- !is.null(d <- dim(x))
  kernel <- kernel[1]
  order <- order[1]
  ret <- .Call(`_smoothemplik_kernelFunCPP`, x, kernel, order, convolution)
  if (is.arr) {
    ret <- array(ret, dim = d, dimnames = dimnames(x))
  } else {
    ret <- drop(ret)
    names(ret) <- names(x)
  }
  return(ret)
}
