#' Smoothed Empirical Likelihood function value
#'
#' Evaluates SEL function for a given moment function at a certain parameter value.
#'
#' @param rho The moment function depending on parameters and data (and potentially other parameters). Must return a numeric vector.
#' @param theta A parameter at which the moment function is evaluated.
#' @param data A data object on which the moment function is computed.
#' @param type Character: \code{"EL"} for empirical likelihood, \code{"EuL"} for Euclidean likelihood, \code{"EL0"} for one-dimensional
#'   empirical likelihood. \code{"EL0"} is *strongly* recommended for 1-dimensional moment functions because it is
#'   faster and more robust: it searches for the Lagrange multiplier directly and has nice fail-safe options
#'   for convex hull failure.
#' @param sel.weights Either a matrix with valid kernel smoothing weights with rows adding up to 1,
#'   or a function that computes the kernel weights based on the \code{data} argument passed to \code{...}.
#' @param EL.args A list of arguments passed to \code{EL()}, \code{EL0()}, or \code{EuL()}.
#' @param kernel.args A list of arguments passed to \code{kernelWeights()} if
#'   \code{sel.weights} is a function.
#' @param minus If TRUE, returns SEL times -1 (for optimisation via minimisation).
#' @param parallel If TRUE, uses \code{parallel::mclapply} to speed up the computation.
#' @param cores The number of cores used by \code{parallel::mclapply}.
#' @param chunks The number of chunks into which the weight matrix is split for memory saving.
#' One chunk is good for sample sizes 2000 and below. If equal to the number of observations, then,
#' the smoothed likelihoods are computed in series, which saves memory but computes kernel weights at
#' every step of a loop, increasing CPU time.
#' If \code{parallel} is \code{TRUE}, parallelisation occurs within each chunk.
#' @param sparse Logical: convert the weight matrix to a sparse one?
#' @param verbose If \code{TRUE}, a progress bar is made to display the evaluation progress in case partial or full memory saving is in place.
#' @param bad.value Replace non-finite individual SEL values with this value.
#'   May be useful if the optimiser does not allow specific non-finite values (like L-BFGS-B).
#' @param attach.attributes If \code{"none"}, returns just the sum of expected likelihoods;
#' otherwise, attaches certain attributes for diagnostics:
#' \code{"ELRs"} for expected likelihoods,
#' \code{"residuals"} for the residuals (moment function values),
#' \code{"lam"} for the Lagrange multipliers lambda in the EL problems,
#' \code{"nabla"} for d/d(lambda)EL (should be close to zero because this must be true for any \code{theta}),
#' \code{"converged"} for the convergence of #' individual EL problems,
#' \code{"exitcode"} for the \code{EL} exit codes (0 for success),
#' \code{"probabilities"} for the matrix of weights (very large, not recommended for sample sizes larger than 2000).
#' @param ... Passed to \code{rho}.
#'
#' @return A scalar with the SEL value and, if requested, attributes containing the diagnostic information attached to it.
#' @export
#'
#' @examples
#' set.seed(1)
#' x <- sort(rlnorm(50))
#' # Heteroskedastic DGP
#' y <- abs(1 + 1*x + rnorm(50) * (1 + x + sin(x)))
#' mod.OLS <- lm(y ~ x)
#' rho <- function(theta, ...) y - theta[1] - theta[2]*x  # Moment fn
#' w <- kernelWeights(x, PIT = TRUE, bw = 0.25, kernel = "epanechnikov")
#' w <- w / rowSums(w)
#' image(x, x, w, log = "xy")
#' theta.vals <- list(c(1, 1), coef(mod.OLS))
#' SEL <- function(b, ...) smoothEmplik(rho = rho, theta = b, sel.weights = w, ...)
#' sapply(theta.vals, SEL) # Smoothed empirical likelihood
#' # SEL maximisation
#' ctl <- list(fnscale = -1, reltol = 1e-6, ndeps = rep(1e-5, 2),
#'             trace = 1, REPORT = 5)
#' b.init <- coef(mod.OLS)
#' b.init <- c(1.790207, 1.007491)  # Only to speed up estimation
#' b.SEL <- optim(b.init, SEL, method = "BFGS", control = ctl)
#' print(b.SEL$par) # Closer to the true value (1, 1) than OLS
#' plot(x, y)
#' abline(1, 1, lty = 2)
#' abline(mod.OLS, col = 2)
#' abline(b.SEL$par, col = 4)
#'
#' # Euclidean likelihood
#' SEuL <- function(b, ...)  smoothEmplik(rho = rho, theta = b,
#'                                        type = "EuL", sel.weights = w, ...)
#' b.SEuL <- optim(coef(mod.OLS), SEuL, method = "BFGS", control = ctl)
#' abline(b.SEuL$par, col = 3)
#' cbind(SEL = b.SEL$par, SEuL = b.SEuL$par)
#'
#' # Now we start from (0, 0), for which the Taylor expansion is necessary
#' # because all residuals at this starting value are positive and the
#' # unmodified EL ratio for the test of equality to 0 is -Inf
#' smoothEmplik(rho=rho, theta=c(0, 0), sel.weights = w, EL.args = list(chull.fail = "none"))
#' smoothEmplik(rho=rho, theta=c(0, 0), sel.weights = w)
#'
#' # The next example is very slow; approx. 1 minute
#' \donttest{
#' # Experiment: a small bandwidth so that the spanning condition should fail often
#' # It yields an appalling estimator
#' w <- kernelWeights(x, PIT = TRUE, bw = 0.15, kernel = "epanechnikov")
#' w <- w / rowSums(w)
#' # The first option is faster but it may sometimes fails
#' b.SELt <- optim(c(0, 0), SEL, EL.args = list(chull.fail = "taylor"),
#'                 method = "BFGS", control = ctl)
#' b.SELw <- optim(c(0, 0), SEL, EL.args = list(chull.fail = "wald"),
#'                 method = "BFGS", control = ctl)
#' # In this sense, Euclidean likelihood is robust to convex hull violations
#' b.SELu <- optim(c(0, 0), SEuL, method = "BFGS", control = ctl)
#' b0grid <- seq(-1.5, 7, length.out = 51)
#' b1grid <- seq(-1.5, 4.5, length.out = 51)
#' bgrid <- as.matrix(expand.grid(b0grid, b1grid))
#' fi <- function(i) smoothEmplik(rho, bgrid[i, ], sel.weights = w, type = "EL0",
#'                   EL.args = list(chull.fail = "taylor"))
#' ncores <- max(floor(parallel::detectCores()/2 - 1), 1)
#' chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")  # Limit to 2 cores for CRAN checks
#' if (nzchar(chk) && chk == "TRUE") ncores <- min(ncores, 2L)
#' selgrid <- unlist(parallel::mclapply(1:nrow(bgrid), fi, mc.cores = ncores))
#' selgrid <- matrix(selgrid, nrow = length(b0grid))
#' probs <- c(0.25, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99, 1-10^seq(-4, -16, -2))
#' levs <- qchisq(probs, df = 2)
#' # levs <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000)
#' labs <- round(levs, 1)
#' cols <- rainbow(length(levs), end = 0.7, v = 0.7)
#' oldpar <- par(mar = c(4, 4, 2, 0) + .1)
#' selgrid2 <- -2*(selgrid - max(selgrid, na.rm = TRUE))
#' contour(b0grid, b1grid, selgrid2, levels = levs,
#'         labels = labs, col = cols, lwd = 1.5, bty = "n",
#'         main = "'Safe' likelihood contours", asp = 1)
#' image(b0grid, b1grid, log1p(selgrid2))
#' # The narrow lines are caused by the fact that if two observations are close together
#' # at the edge, the curvature at that point is extreme
#'
#' # The same with Euclidean likelihood
#' seulgrid <- unlist(parallel::mclapply(1:nrow(bgrid), function(i)
#'   smoothEmplik(rho, bgrid[i, ], sel.weights = w, type = "EuL"),
#'     mc.cores = ncores))
#' seulgrid <- matrix(seulgrid, nrow = length(b0grid))
#' seulgrid2 <- -50*(seulgrid - max(seulgrid, na.rm = TRUE))
#' par(mar = c(4, 4, 2, 0) + .1)
#' contour(b0grid, b1grid, seulgrid2, levels = levs,
#'         labels = labs, col = cols, lwd = 1.5, bty = "n",
#'         main = "'Safe' likelihood contours", asp = 1)
#' image(b0grid, b1grid, log1p(seulgrid2))
#' par(oldpar)
#' }
smoothEmplik <- function(rho, theta, data, sel.weights = NULL,
                         type = c("EL", "EuL", "EL0"),
                         kernel.args = list(bw = NULL, kernel = "epanechnikov", order = 2, PIT = TRUE, sparse = TRUE),
                         EL.args = list(chull.fail = "taylor", weight.tolerance = NULL),
                         minus = FALSE,
                         parallel = FALSE, cores = 1,
                         chunks = NULL, sparse = FALSE, verbose = FALSE,
                         bad.value = -Inf,
                         attach.attributes = c("none", "all", "ELRs", "residuals", "lam", "nabla", "converged", "exitcode", "probabilities"),
                         ...
) {
  type       <- match.arg(type)
  # Constructing residuals
  rho.series <- rho(theta, data, ...)
  n <- NROW(rho.series)  # length or nrow
  if (is.null(chunks)) chunks <- ceiling(n / 2000)
  if (any("none" %in% attach.attributes)) attach.attributes <- "none"
  if (isTRUE(attach.attributes)) attach.attributes <- "all"
  # Special treatment for attach.probs because this is passed to the weighted(EuL|EL0) call to same memory
  attach.probs <- ("probabilities" %in% attach.attributes) | isTRUE(attach.attributes == "all")

  # Since SEL is a non-parametric method and relies on smoothing with kernel weights,
  # using large matrices for large problems can be inefficient
  # Instead, we allow chunking and gradual weight matrix build-up
  if (chunks < 1 || chunks > n)
    stop("smoothEmplik: 'chunks' must be an integer between 1 and the sample size n.")
  if (!is.function(sel.weights) && chunks > 1) {
    stop("smoothEmplik: When chunks > 1, 'sel.weights' must be a function returning weights for given indices (to avoid storing all weights at once).")
  }
  if (!is.null(sel.weights) && !is.function(sel.weights) && !(is.list(sel.weights) ||
        (is.matrix(sel.weights)) || class(sel.weights)[1] %in% c("dgCMatrix", "dgeMatrix")))
    stop("smoothEmplik: 'sel.weights' must be a function, a list, or a matrix of weights.")

  if (chunks == 1) {
    chunk.list <- list(1:n)
  } else if (chunks == n) {
    chunk.list <- as.list(1:n)
  } else {
    group <- cut(1:n, breaks = chunks, labels = FALSE)
    chunk.list <- split(1:n, group)
  }

  empliklist <- vector("list", n)
  if (verbose) pb <- utils::txtProgressBar(min = 0, max = length(chunk.list), style = 3)

  calcOne <- function(i) { # Call the appropriate weighted likelihood function based on `type`
    if (type == "EL") {
      return(EL(z = rho.series, ct = w[i, ], mu = 0, SEL = TRUE,
                        weight.tolerance = EL.args$weight.tolerance, return.weights = attach.probs))
    } else if (type == "EuL") {
      return(EuL(z = rho.series, ct = w[i, ], mu = 0, SEL = TRUE,
                         weight.tolerance = EL.args$weight.tolerance, return.weights = attach.probs))
    } else if (type == "EL0") {
      return(EL0(z = rho.series, ct = w[i, ], mu = 0, SEL = TRUE,
                         chull.fail = EL.args$chull.fail, weight.tolerance = EL.args$weight.tolerance,
                         return.weights = attach.probs))
    } else {
      stop("The 'type' argument must be 'EL' or 'EuL'.")
    }
  }

  if ((is.matrix(sel.weights) || inherits(sel.weights, "dgeMatrix")) && sparse) {
    sel.weights <- Matrix::Matrix(sel.weights, sparse = TRUE)
  }

  for (k in seq_along(chunk.list)) {
    inds <- chunk.list[[k]]
    # This is the memory-consuming operation
    w <- if (!is.function(sel.weights)) sel.weights else suppressWarnings(sel.weights(inds, data))
    if (is.null(dim(w))) w <- matrix(w, nrow = 1)
    if (is.matrix(w) || class(w)[1] %in% c("dgCMatrix", "dgeMatrix")) {
      # w is a weight matrix: row = observation in this batch, column = data point.
      if (nrow(w) != length(inds) || ncol(w) != n)
        stop("smoothEmplik: sel.weights returned a matrix with incompatible dimensions.")
      # This non-pure function relies on the existence of 'w' in the memory
      if (parallel && cores > 1) {
        empliklist[inds] <- parallel::mclapply(X = seq_along(inds), FUN = calcOne, mc.cores = cores)
      } else {
        empliklist[inds] <- lapply(seq_along(inds), calcOne)
      }
    } else {
      stop("smoothEmplik: sel.weights returned an unsupported type (must be function, list, or matrix).")
    }
    if (verbose) utils::setTxtProgressBar(pb, k)
  }
  if (verbose) close(pb)

  log.ELR.values <- unlist(lapply(empliklist, "[[", "logelr"))
  if (any(bad <- !is.finite(log.ELR.values))) log.ELR.values[bad] <- bad.value
  if (minus) log.ELR.values <- -log.ELR.values
  log.SELR <- sum(log.ELR.values)
  ret <- log.SELR
  if (isTRUE(attach.attributes == "none")) return(ret)
  aa <- isTRUE(attach.attributes == "all")
  if ("ELRs" %in% attach.attributes || aa) attr(ret, "ELRs") <- log.ELR.values
  if ("residuals" %in% attach.attributes || aa) attr(ret, "residuals") <- rho.series
  if ("lam" %in% attach.attributes || aa) attr(ret, "lam") <- unlist(lapply(empliklist, "[[", "lam"))
  if ("nabla" %in% attach.attributes || aa)
    attr(ret, "nabla") <- if (type == "EL0") unlist(lapply(empliklist, "[[", "f.root")) else unlist(lapply(empliklist, "[[", "gradnorm"))
  if ("converged" %in% attach.attributes || aa)  # uniroot has many exit codes, but the custom optimiser has only 2 defined by the termination criterion
    attr(ret, "converged") <- if (type == "EL0") unlist(lapply(empliklist, "[[", "converged")) else unlist(lapply(empliklist, "[[", "exitcode")) == 0
  if ("exitcode" %in% attach.attributes || aa) attr(ret, "exitcode") <- unlist(lapply(empliklist, "[[", "exitcode"))
  if (attach.probs) attr(ret, "probabilities") <- lapply(empliklist, "[[", "wts")
  return(ret)
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
