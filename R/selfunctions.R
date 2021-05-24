#' Silverman's rule-of-thumb bandwidth
#'
#' A fail-safe function that would return a nice Silverman-like bandwidth suggestion for data for which
#' the standard deviation might be NA or 0.
#'
#' It is obtained under the assumption that the true density is multivariate normal with zero covariances
#' (i.e. a diagonal variance-covariance matrix
#' \eqn{\Sigma = \mathrm{\mathop{diag}}(\sigma^2_k)}{\Sigma = diag(\sigma^2_k)} with
#' \eqn{\det \Sigma = \prod_k \sigma^2_k}{det \Sigma = prod(\sigma^2_k)} and \eqn{\Sigma^{-1} = diag(\sigma^{-2}_k)}{\Sigma = diag(1/\sigma^2_k)}).
#' Then, the formula 4.12 in Silverman (1986) depends only on \eqn{\alpha}{\alpha}, \eqn{\beta}{\beta}
#' (which depend only on the kernel and are fixed for a multivariate normal), and on the L2-norm of the
#' second derivative of the density. The (i, i)th element of the Hessian of multi-variate normal
#' (\eqn{\phi(x_1, \ldots, x_d) = \phi(X)}{\phi(x_1, ..., x_d) = \phi(X)}) is
#' \eqn{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}{\phi(X)(x_i^2 - \sigma^2_i)/\sigma_i^4}.
#'
#' @param x A numeric vector without non-finite values.
#' @param na.rm Logical: should missing values be removed? Setting it to TRUE may cause issued because variable-wise removal of NAs may return a bandwidth that is inappropriate for the final data set to which it is applied.
#' @return A bandwidth that might be optimal for non-parametric density estimation of \code{x}.
#' @examples
#' set.seed(1); bw.rot(stats::rnorm(100)) # Should be 0.3787568 in R version 4.0.4
#' set.seed(1); bw.rot(matrix(stats::rnorm(500), ncol = 10)) # 0.4737872 ... 0.7089850
#' @export
bw.rot <- function(x, na.rm = FALSE) {
  if (any(is.na(x))) {
    if (na.rm) warning("bw.rot: There are missing values in the data, and you should do something about it because proper analysis is impossible with NA, and your results might be unreliable with these bandwidths.") else
      stop("bw.rot: There are missing values in the data, but non-parametric methods rely on data with finite numeric values only.")
  }
  one.dim <- is.vector(x) # Are our data one-dimensional?
  if (one.dim) x <- matrix(x, ncol = 1)
  d <- ncol(x)
  n <- nrow(x)
  s <- apply(x, 2, function(x) stats::sd(x, na.rm = na.rm))
  AK <- (4 / (2*d + 1))^(1 / (d + 4)) # (4.15) from Silverman (1986)
  if (any(!is.finite(s))) {
    stop("bw.rot: Could not compute the bandwidth, check your data, most likely it has length 1.")
  } else if (all(s > 0)) {
    return(AK * s * length(x)^(-1/(d+4)))
  } else {
    return(rep(1, d))
  }
}

#' Empirical likelihood for one-dimensional vectors
#'
#' Empirical likelihood with counts to solve one-dimensional problems efficiently with Brent's root search algorithm. Conducts an empirical likelihood ratio test of the hypothesis that the mean of \code{z} is \code{mu}
#' The names of the elements in the returned list are consistent with Art B. Owen's original R code.
#' @param z The data vector.
#' @param ct The count variable that indicates the multiplicity ob observations. Can be fractional. Very small counts below the tolerance threshold are zeroed.
#' @param shift The value to add in the denominator (useful in case there are extra Lagrange multipliers): 1 + lambda'Z + shift.
#' @param mu Hypothesized mean of \code{z} in the moment condition.
#' @param SEL If FALSE, then the boundaries for the lambda search are based on the total sum of counts, like in vanilla empirical likelihood,
#' due to Owen (2001, formula 2.9), otherwise according to Cosma et al. (2019, p. 170, the topmost formula).
#' @param wttol Weight tolerance for counts for numerical stability (similar to the ones in Art B. Owen's 2017 code, but adapting to the sample size).
#' @param truncto Counts under \code{wttol} will be set to this value. In most cases, setting this to \code{0} or \code{wttol} is a viable solution of the zero-denominator problem.
#' @param uniroot.control A list passed to the \code{uniroot}.
#' @param return.weights Logical: if TRUE, individual EL weights are computed and returned. Setting this to FALSE gives huge memory savings in large data sets, especially when smoothing is used.
#' @return A list with the following elements:
#' \itemize{
#' \item{logelr }{Logarithm of the empirical likelihood ratio.}
#' \item{lam }{The Lagrange multiplier.}
#' \item{wts }{Observation weights/probabilities (of the same length as \code{z}).}
#' \item{converged }{TRUE if the algorithm converged, FALSE otherwise (usually means that \code{mu} is not in the convex hull of the data, that is, within the range of \code{z}).}
#' \item{iter }{The number of iterations used (from \code{uniroot}).}
#' \item{bracket }{The admissible interval for lambda (that is, yielding weights between 0 and 1).}
#' \item{estim.prec }{The approximate estimated precision of lambda (from \code{uniroot}).}
#' \item{f.root }{The value of the function at the root (from \code{uniroot}); values \code{> sqrt(.Machine$double.eps)} indicate covergence problems.}
#' \item{exitcode }{An integer indicating the reason of termination.}
#' \describe{
#' \item{0:}{success, an interior solution for lambda found.}
#' \item{1:}{the value of the derivative is \code{> sqrt(.Machine$double.eps)} (the value of the function at the returned root it not exactly zero, and this is not an issue with the tolerance).}
#' \item{2:}{the root was found and the function value seems to be zero, but the root is very close (\code{< sqrt(.Machine$double.eps)}) to the boundary.}
#' \item{3:}{like \code{1} and \code{2} simultaneously: the value of the target function is > 1e-8, and the result is close to the boundary (< 1e-8).}
#' \item{4:}{an error occurred while calling \code{uniroot}.}
#' \item{5:}{\code{mu} is not strictly in the convex hull of \code{z} (spanning condition not met).}
#' }
#' }
#' @references
#' Owen, A. B. (2001). *Empirical likelihood*. CRC press.
#'
#' Cosma, A., Kostyrka, A. V., & Tripathi, G. (2019). Inference in conditional moment restriction models when there is selection due to stratification. In *The Econometrics of Complex Survey Data*. Emerald Publishing Limited.
#' @export
weightedEL <- function(z,
                       ct = NULL,
                       shift = NULL,
                       mu = 0,
                       SEL = FALSE,
                       wttol = 0.01 / length(z),
                       truncto = 0,
                       uniroot.control = list(),
                       return.weights = FALSE
) {
  norig <- length(z)
  if (is.null(ct)) ct <- rep(1, norig)
  if (is.null(shift)) shift <- rep(0, norig)
  if (min(ct) < 0) stop("Negative weights are not welcome.")
  # If originally the weights were too small, too many points would be truncated
  # Warn if any non-zero weights are smaller than wttol
  if (any(0 < ct & ct < wttol)) {
    warning(paste0("Positive counts below ", wttol, " have been replaced with ", truncto, "."))
    ct[ct < wttol] <- truncto
  }
  nonz <- which(ct > 0)
  z <- z[nonz]
  ct <- ct[nonz]
  shift <- shift[nonz]
  if (SEL) ct <- ct / sum(ct) # We might have truncated some weights, so re-normalisation is needed!
  # The denominator for EL with counts is the sum of total counts, and for SEL, it is the number of observations!
  N <- if (SEL) norig else sum(ct)
  if (N <= 0) stop("Total weights must be positive.")
  if (is.matrix(z)) {
    if (ncol(z) > 1) stop("Only one-dimensional vectors or matrices are supported.") else z <- as.numeric(z)
  }
  z <- z - mu

  # Simple derivation:
  # p_i = c_i/N * 1/(1 + lambda*z_i + shift)
  # We know that p_i<1 \forall i, so
  # c_i/N * 1/(1 + lambda*z_i + shift) < 1
  # c_i/N < 1 + lambda*z_i + shift
  # c_i/N - 1 - shift < lambda*z_i
  # Two cases are possible: either z_i<0, or z_i>0 (we do not want to divide by z_i=0).
  # Denote negz those i for which z_i<0. Then
  # {(c_i/N - 1 - shift)/z_i > lambda, \forall i: z_i<0
  # {(c_i/N - 1 - shift)/z_i < lambda, \forall i: z_i>0,
  # or
  # {lambda < (c_i/N - 1 - shift)/z_i, \forall i: z_i<0
  # {lambda > (c_i/N - 1 - shift)/z_i, \forall i: z_i>0,
  # which implies
  # {lambda < min_{i: z_i<0} (c_i/N - 1 - shift)/z_i
  # {lambda > max_{i: z_i>0} (c_i/N - 1 - shift)/z_i
  if (return.weights) wts <- rep(0, norig)
  # Checking the spanning condition; 0 must be in the convex hull of z, that is, min(z) < 0 < max(z)
  if (min(z) < 0 & max(z) > 0) {
    negz <- z < 0
    comp <- (ct / N - 1 - shift) / z
    minlambda <- max(comp[!negz])
    maxlambda <- min(comp[negz])
    con <- list(tol = 1e-20, maxiter = 200, trace = 0)
    con[names(uniroot.control)] <- uniroot.control

    dllik <- function(lambda) return(sum(ct * z / (1 + lambda * z + shift)))
    lambda <- tryCatch(stats::uniroot(dllik, interval = c(minlambda, maxlambda), tol = con$tol, maxiter = con$maxiter, trace = con$trace), # If there is a warning, we still return the object
      warning = function(w) return(list(stats::uniroot(dllik, interval = c(minlambda, maxlambda), tol = con$tol, maxiter = con$maxiter, trace = con$trace), w)),
      error = function(e) return(NULL)
    )
    if (!is.null(lambda)) { # Some result with or without a warning as the second element of the list
      if (class(lambda[[2]]) == "simpleWarning") {
        print(lambda[[2]]) # The user should know that went wrong in uniroot()
        lambda <- lambda[[1]]
      }
      lam <- lambda$root
      if (return.weights) wts[nonz] <- (ct / N) / (1 + z * lam + shift)
      logelr <- -sum(ct * log(1 + z * lam + shift))
      converged <- if (class(lambda[[2]]) == "simpleWarning") FALSE else TRUE
      iter <- lambda$iter
      estim.prec <- lambda$estim.prec
      f.root <- lambda$f.root
      exitcode <- 0
      if (abs(f.root) > sqrt(.Machine$double.eps)) {
        warning("FOC not met: the value of d/dlambda(weighted EL) is different from 0 by more than sqrt(Mach.eps)!")
        exitcode <- 1
      }
      if (abs(lam - maxlambda) < sqrt(.Machine$double.eps) | abs(lam - minlambda) < sqrt(.Machine$double.eps)) {
        warning("Lambda is very close to the boundary (closer than sqrt(Mach.eps))!")
        exitcode <- exitcode + 2
      }
    } else {
      lam <- Inf
      logelr <- -Inf
      converged <- FALSE
      iter <- NA
      estim.prec <- NA
      f.root <- NA
      exitcode <- 4
      warning("Root finder returned an error!")
    }
  } else {
    lam <- Inf
    logelr <- -Inf
    converged <- FALSE
    iter <- NA
    minlambda <- -Inf
    maxlambda <- Inf
    estim.prec <- NA
    f.root <- NA
    exitcode <- 5
    warning("mu is not strictly in the convex hull of z!")
  }

  return(list(logelr = logelr, lam = lam, wts = if (return.weights) wts else NULL, converged = converged, iter = iter, bracket = c(minlambda, maxlambda), estim.prec = estim.prec, f.root = f.root, exitcode = exitcode))
}

#' SEL for every observation
#'
#' @param rho The moment function depending on parameters and data (and potentially other parameters). Must return a numeric vector.
#' @param sel.weights A matrix with valid kernel smoothing weights with rows adding up to 1, or a list of kernel weights for smoothing where the sum of each element is 1 (must be returned by \code{sparseVectorToList}).
#' @param parallel If TRUE, uses \code{parallel::mclapply} to speed up the computation.
#' @param cores The number of cores used by \code{parallel::mclapply}.
#' @param grad Currently not used.
#' @param ... Passed to \code{rho}; usually \code{theta} and \code{data}.
#'
#' @return A list with two values:
#' \item{rho}{One-dimensional vector of the moment function values.}
#' \item{empliklist}{A list with the results of one-dimension optimisation conducted by \code{\link{weightedEL}}.}
#' @export
smoothEmplikList <- function(rho,
                             sel.weights = NULL,
                             parallel = FALSE, cores = 1,
                             grad = FALSE,
                             ...
) {
  if (is.null(sel.weights)) {
    data <- list(...)$data
    bw <- bw.rot(data$X)
    sel.weights <- kernelWeights(data$X, bw = bw)
    sel.weights <- sel.weights / rowSums(sel.weights)
    warning(paste0("smoothEmplikList: you forgot to provide SEL weights, assuming normal kernel weights with RoT bandwidth ", round(bw, 5), "!"))
  }

  # Constructing residuals
  rho.series <- rho(...)

  # Attempting to save memory
  if (is.list(sel.weights)) {
    if (parallel) { # Returns a list, one item for each conditioning vector point
      empliklist <- parallel::mclapply(1:nrow(data), function(i) suppressWarnings(weightedEL(rho.series[sel.weights[[i]]$idx], ct = sel.weights[[i]]$ct, SEL = TRUE)), mc.cores = cores)
    } else {
      empliklist <- lapply(sel.weights, function(x) suppressWarnings(weightedEL(rho.series[x$idx], ct = x$ct, SEL = TRUE)))
    }
  } else { # If it is a matrix
    if (parallel) {
      empliklist <- parallel::mclapply(1:nrow(data), function(i) suppressWarnings(weightedEL(rho.series, ct = sel.weights[i, ], SEL = TRUE)), mc.cores = cores)
    } else {
      empliklist <- apply(sel.weights, MARGIN = 1, function(w) suppressWarnings(weightedEL(rho.series, ct = w, SEL = TRUE)))
    }
  }

  return(list(rho = rho.series, empliklist = empliklist))
}

# A function that takes parameters, data, smoothing weights, bandwidth, and returns the smoothed empirical likelihood

#' Smoothed Empirical Likelihood function value
#'
#' Evaluates SEL function for a given moment function at a certain parameter value.
#'
#' @param rho Passed to \code{\link{smoothEmplikList}}.
#' @param sel.weights Passed to \code{\link{smoothEmplikList}}.
#' @param trim A vector of trimming function values to multiply the output of \code{rho(...)} with. If NULL, no trimming is done.
#' @param minus If TRUE, returns SEL times -1 (for optimisation via minimisation).
#' @param parallel Passed to \code{\link{smoothEmplikList}}.
#' @param cores Passed to \code{\link{smoothEmplikList}}.
#' @param bad.value Replace non-finite individual SEL values with this value. May be useful if the optimiser does not allow specific non-finite values.
#' @param ... Passed to \code{\link{smoothEmplikList}}.
#'
#' @return A scalar with the SEL value.
#' @export
#'
#' @examples
smoothEmplik <- function(rho,
                         sel.weights = NULL,
                         trim = NULL,
                         minus = FALSE,
                         parallel = FALSE, cores = 1,
                         bad.value = -Inf,
                         ...
) {
  sellist <- smoothEmplikList(rho = rho, sel.weights = sel.weights, parallel = parallel, cores = cores, ...)
  if (is.null(trim)) trim <- rep(1, length(sellist$rho))
  logsemplik <- sum(trim * unlist(lapply(sellist$empliklist, "[[", "logelr"))) * (1 - 2 * as.numeric(minus))
  if (any(bad <- !is.finite(logsemplik))) logsemplik[bad] <- bad.value
  return(logsemplik)
}

#' Constrained smoothed empirical likelihood
#'
#' This function takes two vectors of parameters: free and fixed ones, so that SEL could be optimised with equality constraints.
#' The fixed parameter vector must have full length, and have NA's in places where the free parameters should go
#; Then, the free parameter values will be inserted in places with NA, and the entire vector passed to smoothEmplik
#' Title
#'
#' @param rho Passed to \code{\link{smoothEmplik}}.
#' @param par.free A numeric vector of *only the free parameters* passed to \code{rho}. In case of constrained optimisation, must be shorter than
#' @param par.fixed A numeric vector of the same length as \code{par.free}: numeric values should be in the places where the values are to be kept fixed, and NA where the free parameters should be. Use NULL or a vector of NA's for unrestricted optimisation.
#' @param sel.weights Passed to \code{\link{smoothEmplik}}.
#' @param trim Passed to \code{\link{smoothEmplik}}.
#' @param minus Passed to \code{\link{smoothEmplik}}.
#' @param parallel Passed to \code{\link{smoothEmplik}}.
#' @param cores Passed to \code{\link{smoothEmplik}}.
#' @param bad.value Passed to \code{\link{smoothEmplik}}.
#' @param ... Passed to \code{\link{smoothEmplik}}.
#'
#' @return A scalar with the constrained (or, if \code{par.fixed = NULL}, unconstrained SEL.)
#' @export
#'
#' @examples
constrSmoothEmplik <- function(rho,
                               par.free = NULL,
                               par.fixed = NULL,
                               sel.weights = NULL, # Must be valid SEL weights with rows adding up to 1!
                               trim = NULL,
                               minus = FALSE, # Return minus log-likelihood for optimisers?
                               parallel = FALSE, cores = 1, # If true and cores>1 (not on Windows), parallel::mclapply is used
                               bad.value = -Inf, # May be useful if an optimiser does not allow specific non-finite values
                               ... # Passed to rho
) {
  if (is.null(par.free)) { # No free parameters supplied
    theta <- par.fixed
  } else if (isTRUE(all(is.na(par.fixed))) | is.null(par.fixed)) { # All free parameters
    theta <- par.free
  } else { # Some free, some fixed parameters
    theta <- par.fixed
    theta[is.na(theta)] <- par.free
  }
  SEL <- smoothEmplik(rho,
    sel.weights = sel.weights, minus = minus, parallel = parallel, cores = cores, ...,
    theta = theta, bad.value = bad.value
  )
  return(SEL)
}

# A function that takes parameters, data, smoothing weights, bandwidth, and returns the gradient of the smoothed empirical likelihood
# smoothEmplikGrad <- function(par,
#                              data, # Must always contain the X variable for construction of weights!
#                              sel.weights = NULL, # Must be valid SEL weights with rows adding up to 1!
#                              trim = NULL,
#                              model = "missing", # "simple" or "missing"
#                              miss.bw = NULL, # Bandwidth for non-parametric imputation in models with missing endogenous variables
#                              minus = FALSE, # Return minus log-likelihood for optimisers?
#                              parallel = FALSE, cores = 1, # If true and cores>1 (not on Windows), mclapply is used
#                              pihat = NULL,
#                              ystarhat = NULL) {
#   sellist <- smoothEmplikList(par = par, data = data, sel.weights = sel.weights, model = model, miss.bw = miss.bw, minus = minus, parallel = parallel, cores = cores, pihat = pihat, ystarhat = ystarhat)
#   if (is.null(trim)) trim <- rep(1, nrow(data))
#   p <- length(par)
#   ans <- rep(0, p)
#
#   rho.j <- sellist$rho
#   rho.j.by.theta <- cbind(1, -data$Z)
#   lambda.i <- unlist(lapply(sellist$empliklist, "[[", "lam"))
#   lambda.i.by.theta <- matrix(0, ncol = p, nrow = nrow(data))
#   sel.denom <- 1 + outer(lambda.i, rho.j)
#   for (th in 1:p) {
#     lambda.i.by.theta[, th] <- rowSums(sel.weights * rho.j.by.theta[, th] / sel.denom^2) / rowSums(sel.weights * rho.j^2 / sel.denom^2)
#     ans[th] <- -sum(outer(trim, trim) * sel.weights * (outer(lambda.i.by.theta[, th], rho.j) + outer(lambda.i, rho.j.by.theta[, th])) / sel.denom)
#   }
#   ans <- ans * (1 - 2 * as.numeric(minus))
#   return(ans)
# }

# smoothEmplik(c(1, 1), data = data, sel.weights = sel.weights, trim = NULL, model = "missing", pihat = pihat, ystarhat = ystarhat)
# smoothEmplikGrad(par, data = data, sel.weights = sel.weights, trim = NULL, model = "missing", pihat = pihat, ystarhat = ystarhat)
# library(numDeriv)
# grad(function(x) smoothEmplik(x, data = data, sel.weights = sel.weights, trim = NULL, model = "missing", pihat = pihat, ystarhat = ystarhat), par)

# Function: optimise SEL with or without constraint, starting with smart initial values for SEL using GMM
# In order to speed up convergence and reduce maximisation time, we check a grid of points not far away from GMM initval and pick the best one
# This function covers a wide grid of possible points
# Receives the dataset from the external function and does the computaion
#' Title
#'
#' @param rho Passed to \code{\link{smoothEmplik}} and \code{\link{constrSmoothEmplik}}.
#' @param start.values Initial values for unrestricted parameters.
#' @param restricted.params Fixed parameter values used for constrained optimisation, hypothesis testing, power and size calculations etc.
#' @param sel.weights Passed to \code{\link{smoothEmplik}} and \code{\link{constrSmoothEmplik}}
#' @param verbose If TRUE, reports optimisation progress, otherwise remains silent.
#' @param parallel Passed to \code{\link{smoothEmplik}} and \code{\link{constrSmoothEmplik}}.
#' @param cores Passed to \code{\link{smoothEmplik}} and \code{\link{constrSmoothEmplik}}.
#' @param optmethod A string indicating the optimisation method: "nlm" (the default, invoking \code{stats::nlm}) is slightly quicker, and if it fails, "BFGS" is used
#' @param nlm.step.max Passed to \code{nlm} if \code{optmethod == "nlm"}; in case of convergence issues, can be reduced, but if it is too small and 5 \code{nlm} iterations are done with the maximum step size, optimisation is re-done via BFGS.
#' @param maxit Maximum number of numerical optimiser steps. If it has not converged in this number of steps, fail gracefully with a meaningfull return.
#' @param ... Passed to \code{\link{smoothEmplik}} and \code{\link{constrSmoothEmplik}}.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{par }{The argmax of the SEL that was found by the optimiser.}
#' \item{value }{The SEL value corresponding to \code{par}.}
#' \item{restricted }{A logical vector of the same length as \code{par} indicating if the respective element was kept fixed during the optimisation.}
#' \item{xtimes }{A vector with two components indicating time (in seconds) it took to find the candidate initial value and the final optimum respectively.}
#' }
#' @export
#'
#' @examples
maximiseSEL <- function(rho,
                        start.values = NULL,
                        restricted.params = NULL,
                        sel.weights = NULL,
                        verbose = FALSE,
                        parallel = FALSE,
                        cores = 1,
                        optmethod = c("nlm", "BFGS", "Nelder-Mead"),
                        nlm.step.max = NULL,
                        maxit = 50,
                        ...) {
  optmethod <- optmethod[1]
  tic0 <- Sys.time()

  # Case 1: If all parameters are restricted, just evaluate the SEL at the given point
  if (!is.null(restricted.params) & all(is.finite(restricted.params))) {
    SEL <- smoothEmplik(rho = rho, sel.weights = sel.weights, parallel = parallel, cores = cores, theta = restricted.params, ...)
    diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
    return(list(par = restricted.params, value = SEL, restricted = rep(TRUE, length(restricted.params)), code = 0, xtimes = c(initial = 0, opt = diff.opt)))
  }

  # Case 2: # If the model is estimated without restrictions; safeguarging against logical(0)
  if (isTRUE(all(is.na(restricted.params))) | is.null(restricted.params)) {
    restricted <- rep(FALSE, length(start.values))
  } else {
    restricted <- !is.na(restricted.params)
  }

  if (verbose) print("Maximising SEL...")
  SELToOptim <- function(theta, ...) constrSmoothEmplik(rho, par.free = theta, par.fixed = restricted.params, sel.weights = sel.weights, minus = TRUE, parallel = parallel, cores = cores, ...)
  # constrSmoothEmplik(rho = rho.complete.case, par.free = start.values, par.fixed = restricted.params, sel.weights = sel.weights, minus = TRUE, parallel = parallel, cores = cores, data = data)
  if (optmethod == "nlm") {
    if (is.null(nlm.step.max)) nlm.step.max <- sqrt(sum(start.values^2)) / 3
    optim.SEL <- tryCatch(stats::nlm(p = start.values, f = SELToOptim, print.level = verbose * 2, stepmax = nlm.step.max, ...), error = opterror)
    if (optim.SEL$code %in% c(4, 5)) { # If there are optimisation issues
      warning(paste0("nlm exit code: ", optim.SEL$code, ", restarting with BFGS."))
      optmethod <- "BFGS"
      optim.SEL <- tryCatch(stats::optim(par = start.values, fn = SELToOptim, control = list(trace = as.numeric(verbose), REPORT = 1, maxit = maxit), method = optmethod, ...), error = opterror)
    }
  } else {
    optim.SEL <- tryCatch(stats::optim(par = start.values, fn = SELToOptim, control = list(trace = as.numeric(verbose), REPORT = 1, maxit = maxit), method = optmethod, ...), error = opterror)
  }

  diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
  if (any(restricted)) {
    thetahat <- restricted.params
    thetahat[is.na(thetahat)] <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  } else {
    thetahat <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  }
  SEL <- if (optmethod != "nlm") -optim.SEL$value else -optim.SEL$minimum
  if (!isTRUE(abs(SEL) < 1e8)) {
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

#' Fail gracefully
#'
#' Return parseable output in case the optimiser fails
#'
#' @return A list of the same format as the output of \code{\link{maximiseSEL}}.
fail <- function() {
  return(list(par = c(NA, NA), value = NA, restricted = c(FALSE, FALSE), code = NA, iterations = NA, xtimes = c(0, 0)))
}


#' Convert a weight vector to list
#'
#' This function saves memory (which is crucial in large samples) and allows one to speed up the code by minimising the number of
#' time-consuming subsetting operations and memory-consuming matrix multiplications. We do not want to rely on extra packages for
#' sparse matrix manipulation since the EL smoothing weights are usually fixed at the beginning, and need not be recomputed dynamically,
#' so we recommend applying this function to the rows of a matrix. In order to avoid numerical instability, the weights are trimmed
#' at \code{0.01 / length(x)}.
#'
#' @param x A numeric vector (with many close-to-zero elements).
#' @param trim A trimming function that returns a threshold value below which the weights are ignored. In common applications, this function should tend to 0 as the length of \code{x} increases.
#'
#' @return A list with indices of large enough elements.
#' @export
#'
#' @examples
sparseVectorToList <- function(x, trim = function(x) 0.01 / length(x)) {
  idx <- which(x >= trim(x))
  x <- x[idx]
  return(list(idx = idx, ct = x))
}

jitterText <- function(x, y, labels, times = 16, radius = 0.1, bgcol = "#FFFFFF88", col = "black", ...) {
  dirx <- cos(seq(0, 2*pi, length.out = times + 1)[-1])
  diry <- sin(seq(0, 2*pi, length.out = times + 1)[-1])
  for (i in 1:times) graphics::text(x + dirx[i] * radius, y + diry[i] * radius, labels, col = bgcol, ...)
  graphics::text(x, y, labels, col = col, ...)
}


#' Newey's optimal IV estimator (1993) for linear models
#'
#' @param y The vector containing dependent variable values
#' @param incl A vector or a matrix of included instruments (can be NULL)
#' @param endog A vector or a matrix of endogenous variables (can be NULL)
#' @param excl A vector or a matrix of excluded instruments (can be NULL)
#' @param bw A scalar of a vector or bandwidths to predict sigma2(incl, excl) (can be NULL)
#' @param iterations Iterations or Newey's weighting procedure; offers some finite-sample (but not asymptotic) improvement if greater than 1
#' @param coef.names A character vector to assign as the names of the estimated coefficients (useful if some variables are endogenous and their projections inherit uninformative names)
#'
#' @return A list with optimal IV estimates, their variance-covariance matrix, and model residuals.
#' @export
#'
#' @examples
#' d <- generateData(seed = 1, n = 500)
#' d <- d[as.logical(d$D), ]
#' incl <- NULL
#' endog <- d$Z
#' excl <- d$X
#' y <- d$Y
#' lmEff(y, incl, endog, excl)
#' m1 <- AER::ivreg(y ~ endog | excl)
#' m2 <- lmEff(y, incl, endog, excl, iterations = 0)
lmEff <- function(y, incl = NULL, endog = NULL, excl = NULL, bw = NULL, iterations = 2, coef.names = NULL) {
  if (is.null(incl) & is.null(endog) & is.null(excl)) stop("There is no model; no variables on the RHS!")
  IV <- cbind(1, incl, excl)
  IV0 <- cbind(incl, excl)
  if (!is.null(endog)) {
    if (is.vector(endog)) endog <- as.matrix(endog)
    if (is.null(excl)) stop("There must be excluded instruments!") else {
      if (is.vector(excl)) excl <- as.matrix(excl)
      if (ncol(endog) > ncol(excl)) stop("There must be more excluded instruments than endogenous variables!")
    }
    if (is.vector(excl)) excl <- as.matrix(excl)
    Xhat <- stats::lm.fit(x = IV, y = endog)$fitted.values
    Xhat <- cbind(incl, Xhat)
    colnames(Xhat) <- c(names(incl), names(endog))
  } else Xhat <- incl
  mod <- stats::lm(y ~ Xhat)
  if (iterations < 1) {
    warning("No iterations for optimal instrument estimation; returning the vanilla inefficient IV estimator with the wrong VCOV!")
  } else {
    bw0 <- bw
    for (i in 1:iterations) {
      e2 <- as.numeric(y - cbind(1, incl, endog) %*% mod$coefficients)^2
      if (is.null(bw0)) bw <- bw.CV(x = IV0, y = e2, CV = "LSCV", degree = 0, robust.iterations = 0)
      mod.var <- kernelSmooth(x = IV0, y = e2, bw = bw, degree = 0)
      # plot(IV0[, 1], e2)
      # points(IV0[, 1], mod.var, cex = 0.7, col = "red", pch = 16)
      mod <- stats::lm(y ~ Xhat, weights = 1 / mod.var)
    }
  }
  names(mod$coefficients) <- coef.names
  v <- stats::vcov(mod)
  rownames(v) <- colnames(v) <- coef.names
  U <- as.numeric(y - cbind(1, incl, endog) %*% mod$coefficients)
  return(list(coefficients = mod$coefficients, vcov = v, residuals = U))
}

getSELWeights <- function(X, bw) {
  sel.weights <- kernelWeights(X, bw = bw)
  sel.weights <- sel.weights / rowSums(sel.weights)
  sel.weights <- apply(sel.weights, 1, sparseVectorToList)
}

opterror <- function(e = NULL) return(list(code = 5))
