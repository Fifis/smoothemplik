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
#' @param n.orig An optional scalar to denote the original sample size (useful in the rare cases re-normalisation is needed).
#' @param weight.tolerance Weight tolerance for counts for numerical stability (similar to the ones in Art B. Owen's 2017 code, but adapting to the sample size).
#' @param truncto Counts under \code{weight.tolerance} will be set to this value. In most cases, setting this to \code{0} or \code{weight.tolerance} is a viable solution of the zero-denominator problem.
#' @param uniroot.control A list passed to the \code{uniroot}.
#' @param return.weights Logical: if TRUE, individual EL weights are computed and returned. Setting this to FALSE gives huge memory savings in large data sets, especially when smoothing is used.
#' @param verbose Logical: if \code{TRUE}, prints warnings.
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
                       n.orig = NULL,
                       weight.tolerance = 0.01 / length(z),
                       truncto = 0,
                       uniroot.control = list(),
                       return.weights = FALSE,
                       verbose = FALSE
) {
  if (is.null(n.orig)) n.orig <- length(z)
  if (is.null(ct)) ct <- rep(1, length(z))
  if (is.null(shift)) shift <- rep(0, length(z))
  if (min(ct) < 0) stop("Negative weights are not welcome.")
  # If originally the weights were too small, too many points would be truncated
  # Warn if any non-zero weights are smaller than weight.tolerance
  if (any(0 < ct & ct < weight.tolerance)) {
    if (verbose) warning(paste0("Positive counts below ", weight.tolerance, " have been replaced with ", truncto, "."))
    ct[ct < weight.tolerance] <- truncto
  }
  nonz <- which(ct > 0)
  z <- z[nonz]
  ct <- ct[nonz]
  shift <- shift[nonz]
  if (SEL) ct <- ct / sum(ct) # We might have truncated some weights, so re-normalisation is needed!
  # The denominator for EL with counts is the sum of total counts, and for SEL, it is the number of observations!
  N <- if (SEL) n.orig else sum(ct)
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
  wts <- if (return.weights) rep(0, n.orig) else NULL
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
        if (verbose) print(lambda[[2]]) # The user should know that went wrong in uniroot()
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
        if (verbose) warning("FOC not met: the value of d/dlambda(weighted EL) is different from 0 by more than sqrt(Mach.eps)!")
        exitcode <- 1
      }
      if (abs(lam - maxlambda) < sqrt(.Machine$double.eps) | abs(lam - minlambda) < sqrt(.Machine$double.eps)) {
        if (verbose) warning("Lambda is very close to the boundary (closer than sqrt(Mach.eps))!")
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
      if (verbose) warning("Root finder returned an error!")
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
    if (verbose) warning("mu is not strictly in the convex hull of z!")
  }

  return(list(logelr = logelr, lam = lam, wts = wts, converged = converged, iter = iter, bracket = c(minlambda, maxlambda), estim.prec = estim.prec, f.root = f.root, exitcode = exitcode))
}

#' Smoothed Empirical Likelihood function value
#'
#' Evaluates SEL function for a given moment function at a certain parameter value.
#'
#' @param rho The moment function depending on parameters and data (and potentially other parameters). Must return a numeric vector.
#' @param theta A parameter at which the moment function is evaluated.
#' @param data A data object on which the moment function is computed.
#' @param sel.weights Either a matrix with valid kernel smoothing weights with rows adding up to 1, or a list of kernel weights for smoothing where the sum of each element is 1 (must be returned by \code{sparseVectorToList}), or a function that computes the kernel weights based on the \code{data} argument passed to \code{...}. If \code{memory.saving} is \code{"partial"} or \code{"full"}, then it must be a function that computes the kernel weights for the data set.
#' @param trim A vector of trimming function values to multiply the output of \code{rho(...)} with. If NULL, no trimming is done.
#' @param weight.tolerance Passed to \code{weightedEL} (uses the same default value).
#' @param minus If TRUE, returns SEL times -1 (for optimisation via minimisation).
#' @param parallel If TRUE, uses \code{parallel::mclapply} to speed up the computation.
#' @param cores The number of cores used by \code{parallel::mclapply}.
#' @param memory.saving A string. \code{"none"} implies no memory-saving tricks, and the entire problem is processed in the computer memory at once (good for sample sizes 2000 and below; if \code{sel.weights} is not provided or is a function, the weight matrix / list is computed at once.). If \code{"full"}, then the smoothed likelihoods are computed in series, which saves memory but computes kernel weights at every step of a loop, increasing CPU time; the SEL weights, normally found in the rows of the \code{sel.weights} matrix, are computed on the fly. If \code{"partial"}, then, the problem is split into \code{chunks} sub-problems with smaller weight matrices / lists. If \code{parallel} is \code{TRUE}, parallelisation occurs within each chunk.
#' @param chunks The number of chunks into which the weight matrix is split. Only used if \code{memory.saving} is \code{"partial"}. If there are too many chunks (resulting in fewer than 2 observations per chunk), then it is treated as if \code{memory.saving} were \code{"full"}.
#' @param print.progress If \code{TRUE}, a progress bar is made to display the evaluation progress in case partial or full memory saving is in place.
#' @param bad.value Replace non-finite individual SEL values with this value. May be useful if the optimiser does not allow specific non-finite values (like L-BFGS-B).
#' @param attach.attributes If \code{"none"}, returns just the sum of expected likelihoods; otherwise, attaches certain attributes for diagnostics: \code{"ELRs"} for expected likelihoods, \code{"residuals"} for the residuals (moment function values), \code{"lam"} for the Lagrange multipliers in the EL problems, \code{"converged"} for the convergence of individual EL problems, \code{"exitcode"} for the \code{weightedEL} exit codes (0 for success), \code{"probabilities"} for the matrix of weights (very large, not recommended for sample sizes larger than 2000).
#' @param ... Passed to \code{rho}.
#'
#' @return A scalar with the SEL value and, if requested, attributes containing the diagnostic information attached to it.
#' @export
#'
#' @examples
smoothEmplik <- function(rho,
                         theta,
                         data,
                         sel.weights = NULL,
                         trim = NULL, weight.tolerance = NULL,
                         minus = FALSE,
                         parallel = FALSE, cores = 1,
                         memory.saving = c("none", "full", "partial"), chunks = 10, print.progress = FALSE,
                         bad.value = -Inf,
                         attach.attributes = c("none", "all", "ELRs", "residuals", "lam", "converged", "exitcode", "probabilities"),
                         ...
) {
  # Constructing residuals
  rho.series <- rho(theta, data, ...)
  n <- length(rho.series)
  if (is.null(weight.tolerance)) weight.tolerance <- 0.01 / n
  if (any("none" %in% attach.attributes)) attach.attributes <- "none"
  attach.probs <- ("probabilities" %in% attach.attributes) | isTRUE(attach.attributes == "all")

  # Since SEL is a non-parametric method and relies on smoothing with kernel weights, using large matrices for large problems
  # can be terribly inefficient. Instead, we split it into chunks.
  memory.saving <- memory.saving[1]
  empliklist <- vector("list", n)
  if (memory.saving == "full") {
      if (!is.function(sel.weights)) stop("smoothEmplik: sel.weights must be a function that accepts an index and a data frame and returns a numeric vector of kernel weights.")
    if (print.progress) pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
      for (i in 1:n) {
        w <- suppressWarnings(as.numeric(sel.weights(i, data)))
        empliklist[[i]] <- weightedEL(rho.series, ct = w, SEL = TRUE, weight.tolerance = weight.tolerance, return.weights = attach.probs)
        if (print.progress) {if (i %% floor(n / 50) == 0) utils::setTxtProgressBar(pb, i)}
      }
    if (print.progress) close(pb)
    } else {
      if (memory.saving == "none") {
        chunks <- 1
        chunk.list <- list(1:n)
      } else if (memory.saving == "partial") {
        if (chunks > n/2) {
          chunks <- max(ceiling(n/1000), 2)
          warning(paste0("Too many chunks for memory saving (fewer than 2 observations per chunk); to reduce overhead, using chunks = ", ceiling(n/1000), "."))
        }
        chunk.labels <- cut(1:n, breaks = chunks, labels = 1:chunks)
        chunk.list <- split(1:n, chunk.labels)
      } else stop("smoothEmplik: Wrong value of the memory.saving argument.")
      if (print.progress) pb <- utils::txtProgressBar(min = 0, max = chunks, style = 3)
      for (k in 1:chunks) {
        inds <- chunk.list[[k]]
        if (k == 1 & memory.saving == "none" & (is.matrix(sel.weights) | is.list(sel.weights))) {
          w <- sel.weights
        } else w <- suppressWarnings(sel.weights(chunk.list[[k]], data))

        # If w is a sparse matrix, it can be converted into a list to save memory (helps in most cases)
        ELiFunc <- if (is.list(w)) function(i) weightedEL(rho.series[w[[i]]$idx], ct = w[[i]]$ct, SEL = TRUE, weight.tolerance = weight.tolerance, return.weights = attach.probs) else if (is.matrix(w)) function(i) weightedEL(rho.series, ct = w[i, ], SEL = TRUE, weight.tolerance = weight.tolerance, return.weights = attach.probs) else {
          if (memory.saving == "none") stop("smoothEmplik: sel.weights must be a list or a matrix of weights for every observation of the data set, or a function returning a matrix or a list (favouring CPU over memory).") else stop("smoothEmplik: sel.weights must be a function returning a matrix or a list (favouring memory over CPU).")
        }
        empliklist[inds] <- if (parallel & cores > 1) parallel::mclapply(X = seq_along(inds), FUN = ELiFunc, mc.cores = cores) else lapply(seq_along(inds), ELiFunc)
        if (print.progress) utils::setTxtProgressBar(pb, k)
      }
      if (print.progress) close(pb)
    }

  if (is.null(trim)) trim <- rep(1, n)
  log.ELR.values <- unlist(lapply(empliklist, "[[", "logelr"))
  if (any(bad <- !is.finite(log.ELR.values))) log.ELR.values[bad] <- bad.value
  if (minus) log.ELR.values <- -log.ELR.values
  log.SELR <- sum(trim * log.ELR.values)
  ret <- log.SELR
  if (isTRUE(attach.attributes == "none")) return(ret)
  if ("ELRs" %in% attach.attributes | isTRUE(attach.attributes == "all")) attr(ret, "ELRs") <- log.ELR.values
  if ("residuals" %in% attach.attributes | isTRUE(attach.attributes == "all")) attr(ret, "residuals") <- rho.series
  if ("lam" %in% attach.attributes | isTRUE(attach.attributes == "all")) attr(ret, "lam") <- unlist(lapply(empliklist, "[[", "lam"))
  if ("converged" %in% attach.attributes | isTRUE(attach.attributes == "all")) attr(ret, "converged") <- unlist(lapply(empliklist, "[[", "converged"))
  if ("exitcode" %in% attach.attributes | isTRUE(attach.attributes == "all")) attr(ret, "exitcode") <- unlist(lapply(empliklist, "[[", "exitcode"))
  if (attach.probs) attr(ret, "probabilities") <- lapply(empliklist, "[[", "wts")
  return(ret)
}

#' Constrained smoothed empirical likelihood
#'
#' This function takes two vectors of parameters: free and fixed ones, so that SEL could be optimised with equality constraints.
#' The fixed parameter vector must have full length, and have NA's in places where the free parameters should go.
#' Then, the free parameter values will be inserted in places with NA, and the entire vector passed to smoothEmplik.
#'
#' @param par.free A numeric vector of *only the free parameters* passed to \code{rho}. In case of constrained optimisation, must be shorter than the number of parameters in the model.
#' @param par.fixed A numeric vector of the same length as \code{par.free}: numeric values should be in the places where the values are to be kept fixed, and NA where the free parameters should be. Use NULL or a vector of NA's for unrestricted optimisation.
#' @param ... Passed to \code{\link{smoothEmplik}}.
#'
#' @return A scalar with the constrained (or, if \code{par.fixed = NULL}, unconstrained SEL.)
#' @export
#'
#' @examples
constrSmoothEmplik <- function(par.free = NULL,
                               par.fixed = NULL,
                               ...
) {
  if (is.null(par.free)) { # No free parameters supplied
    theta <- par.fixed
  } else if (isTRUE(all(is.na(par.fixed))) | is.null(par.fixed)) { # All free parameters
    theta <- par.free
  } else { # Some free, some fixed parameters
    theta <- par.fixed
    theta[is.na(theta)] <- par.free
  }
  SEL <- smoothEmplik(theta = theta, ...)
  return(SEL)
}

#' Optimise SEL with or without constraints
#'
#' @param start.values Initial values for unrestricted parameters.
#' @param restricted.params Fixed parameter values used for constrained optimisation, hypothesis testing, power and size calculations etc.
#' @param verbose If TRUE, reports optimisation progress, otherwise remains silent.
#' @param optmethod A string indicating the optimisation method: "nlm" (the default, invoking \code{stats::nlm}) is slightly quicker, and if it fails, "BFGS" is used
#' @param nlm.step.max Passed to \code{nlm} if \code{optmethod == "nlm"}; in case of convergence issues, can be reduced, but if it is too small and 5 \code{nlm} iterations are done with the maximum step size, optimisation is re-done via BFGS.
#' @param maxit Maximum number of numerical optimiser steps. If it has not converged in this number of steps, fail gracefully with a meaningfull return.
#' @param ... Passed to \code{\link{constrSmoothEmplik}} or, in case all parameters are fixed, \code{\link{smoothEmplik}}.
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
maximiseSEL <- function(start.values = NULL,
                        restricted.params = NULL,
                        verbose = FALSE,
                        optmethod = c("nlm", "BFGS", "Nelder-Mead"),
                        nlm.step.max = NULL,
                        maxit = 50,
                        ...) {
  optmethod <- optmethod[1]
  tic0 <- Sys.time()

  # Case 1: If all parameters are restricted, just evaluate the SEL at the given point
  if (!is.null(restricted.params) & all(is.finite(restricted.params))) {
    SEL <- smoothEmplik(theta = restricted.params, ...)
    diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
    return(list(par = restricted.params, value = SEL, restricted = rep(TRUE, length(restricted.params)), code = 0, xtimes = c(initial = 0, opt = diff.opt)))
  }

  # Case 2: # If the model is estimated without restrictions; safeguarging against logical(0)
  if (isTRUE(all(is.na(restricted.params))) | is.null(restricted.params)) {
    restricted <- rep(FALSE, length(start.values))
  } else {
    restricted <- !is.na(restricted.params)
  }

  if (verbose) cat("Maximising SEL.\n")
  opt.controls <- list(trace = as.numeric(verbose), REPORT = 1, maxit = maxit, fnscale = -1)
  if (!is.null(list(...)$minus)) stop("Do not pass 'minus' to 'maximiseSEL'.")
  SELToOptim <- function(theta, ...) constrSmoothEmplik(par.free = theta, par.fixed = restricted.params, minus = FALSE, ...)
  if (optmethod == "nlm") {
    if (is.null(nlm.step.max)) nlm.step.max <- sqrt(sum(start.values^2)) / 3
    optim.SEL <- tryCatch(stats::nlm(p = start.values, f = SELToOptim, fscale = -1, print.level = verbose * 2, stepmax = nlm.step.max, ...), error = optError)
    if (optim.SEL$code %in% c(4, 5)) { # If there are optimisation issues
      warning(paste0("nlm exit code: ", optim.SEL$code, ", restarting with BFGS."))
      optmethod <- "BFGS"
      optim.SEL <- tryCatch(stats::optim(par = start.values, fn = SELToOptim, control = opt.controls, method = optmethod, ...), error = optError)
    }
  } else {
    optim.SEL <- tryCatch(stats::optim(par = start.values, fn = SELToOptim, control = opt.controls, method = optmethod, ...), error = optError)
  }

  diff.opt <- as.numeric(difftime(Sys.time(), tic0, units = "secs"))
  if (any(restricted)) {
    thetahat <- restricted.params
    thetahat[is.na(thetahat)] <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  } else {
    thetahat <- if (optmethod != "nlm") optim.SEL$par else optim.SEL$estimate
  }
  SEL <- if (optmethod != "nlm") optim.SEL$value else optim.SEL$minimum
  if (!isTRUE(abs(SEL) < 1e20)) {
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
#' at \code{0.01 / length(x)}. Using too much trimming may cause the spanning condition to fail (the moment function values can have the same sign in some neighbourhoods).
#'
#' @param x A numeric vector (with many close-to-zero elements).
#' @param trim A trimming function that returns a threshold value below which the weights are ignored. In common applications, this function should tend to 0 as the length of \code{x} increases.
#'
#' @return A list with indices of large enough elements.
#' @export
#'
#' @examples
sparseVectorToList <- function(x, trim = NULL) {
  if (is.null(trim)) trim <- function(x) 0.01 / length(x)
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

getSELWeights <- function(X, bw, trim = NULL) {
  sel.weights <- kernelWeights(X, bw = bw)
  sel.weights <- sel.weights / rowSums(sel.weights)
  sel.weights <- apply(sel.weights, 1, sparseVectorToList, trim = trim)
}

optError <- function(e = NULL) return(list(code = 5))

checkRes <- function(res, split, min.opp.sign = 3) {
  r <- split(res, split)
  ret <- unlist(lapply(r, function(x) {
    if (all(is.na(x))) return("All NA") else
      if (sum(x > 0, na.rm = TRUE) >= min.opp.sign & sum(x < 0, na.rm = TRUE) >= min.opp.sign) return("OK") else
        if (max(x, na.rm = TRUE) > 0 & min(x, na.rm = TRUE) < 0) return("Weak") else return("Fail")
  }))
  ret <- factor(ret, levels = c("All NA", "OK", "Weak", "Fail"))
  return(ret)
}

# A function to check the spanning condition for every distinct value defining the support
# Returns indices defining the bins in which the data fall according to conditioning on discrete variables
# Empirical criterion: grab the closest until there are at least `min.opp.sign` positive
# and `min.opp.sign` negative residuals.
# Closest = by Mahalanobis distance
defineValidSupport <- function(res, data, var.names, type = c("spanning", "cellsize"), min.obs = 3, verbose = TRUE) {
  type <- type[1]
  if (!(type %in% c("spanning", "cellsize", "none"))) stop("defineValidSupport: type must be 'spanning' to check the spanning condition or 'cellsize' to simply ensure large enough cells.")
  spl <- as.integer(interaction(data[, var.names], drop = TRUE, sep = ";"))
  ch <- checkRes(res = res, split = spl, min.opp.sign = min.obs)
  group.counts <- table(spl)
  bad <- switch(type, spanning = ch %in% c("Fail", "Weak"), cellsize = group.counts < min.obs)
  na  <- switch(type, spanning = ch %in% "All NA", cellsize = rep(FALSE, length(group.counts)))

  if (verbose) {
    if (type == "spanning" & any(bad)) {
      cat("The spanning condition necessary for SEL does not hold for ", sum(ch == "Fail"), " groups (", sum(group.counts[ch == "Fail"]), " observations) in the data set.\n", sep = "")
      cat("Additionally, ", sum(ch == "Weak"), " groups (", sum(group.counts[ch == "Weak"]), " observations) are in the groups with fewer than ",  min.obs, " obs. of opposite sign, creating potential numerical instabilities.\n", sep = "")
    } else if (type == "cellsize" & any(bad)) {
      cat("The minimum cell size (", min.obs, ") condition does not hold for ", sum(bad), " groups (", sum(group.counts[bad]), " observations) in the data set.\n", sep = "")
    } else if (min.obs <= 1) cat("No restrictions requested, skipping the check and returning the group indices as they are!\n")
  }
  if ((type == "spanning" & all(ch == "OK")) | (type == "cellsize" & all(!bad)) | min.obs <= 1) {
    return(spl)
  } else {
    if (any(bad)) {
      bad.cats <- which(bad)
      na.cats <- which(na)
      n.uniq <- as.integer(max(spl))
      if (!all(1:n.uniq == 1:max(spl))) stop("defineValidSupport: internal error while generating split.")

      mah.var <- stats::var(data[, var.names])
      if (any(is.na(mah.var))) stop("Missingness is not allowed in the variables defining the cells.")
      mah.v_1 <- solve(mah.var)

      d.list <- split(data[, var.names], f = spl)
      group.means <- as.matrix(do.call(rbind, lapply(d.list, colMeans)))
      group.distances <- lapply(1:n.uniq, function(i) {
        di <- sweep(group.means, 2, group.means[i, ], "-")
        return(rowSums(di %*% mah.v_1 * di))
      })
      group.distances <- do.call(rbind, group.distances)
      diag(group.distances) <- Inf
      # image(sqrt(group.distances), asp = 1) # Look at this for debugging to see if there is any internal structure
      # Now we declare that the distances between bad groups are infinite
      # so that the matches were sought only in the good ones
      group.distances[bad.cats, bad.cats] <- Inf
      # Also, the ones where all dependent variable observations are missing are not good for bin enlargement
      group.distances[bad.cats, na.cats] <- Inf
      group.distances[na.cats, bad.cats] <- Inf

      # image(sqrt(group.distances), asp = 1)
      bad.cats.closest <- apply(group.distances[bad.cats, , drop = FALSE], 1, which.min)

      spl.new <- spl
      for (i in seq_along(bad.cats)) spl.new[spl == bad.cats[i]] <- bad.cats.closest[i]
      spl.new <- as.integer(factor(spl.new)) # Making numeration contiguous
      if (verbose) {
        cat("Before: ", n.uniq, " groups, after: ", length(unique(spl.new)), " groups after moving ", sum(group.counts[bad]), " observations (", round(sum(group.counts[bad])/nrow(data) * 100), "%) to larger groups.\n", sep = "")
      }
      return(spl.new)
    }
  }
}

smoothEmplikDiscrete <- function(rho,
                                 theta, data,
                                 by = NULL, parallel = FALSE, cores = 1,
                                 trim = NULL,
                                 minus = FALSE,
                                 bad.value = -Inf,
                                 weight.tolerance = NULL,
                                 attach.attributes = FALSE,
                                 ...) {
  if (is.null(by)) stop("You forgot to supply the factor or integer variable 'by' indicating the unique values of the conditioning set.")

  rho.series <- rho(theta, data, ...)
  n <- length(rho.series)
  rho.list <- split(rho.series, f = by)
  if (is.null(weight.tolerance)) weight.tolerance <- 0.01 / n

  SELi <- function(x) {
    if (all(is.na(x))) return(list(logelr = 0, lam = 0, wts = rep(1 / length(x), length(x)), converged = TRUE, iter = 0, bracket = c(0, 0), estim.prec = NA, f.root = NA, exitcode = 0))
    x <- x[!is.na(x)]
    return(suppressWarnings(weightedEL(x, SEL = TRUE, weight.tolerance = weight.tolerance)))
  }

  # Attempting to save memory
  if (parallel) { # Returns a list, one item for each conditioning vector point
    empliklist <- parallel::mclapply(rho.list, SELi, mc.cores = cores)
  } else {
    empliklist <- lapply(rho.list, SELi)
  }

  SELs <- unlist(lapply(empliklist, "[[", "logelr"))
  SELs.all <- unname(SELs[by])
  if (is.null(trim)) trim <- rep(1, n)
  SEL <- sum(trim * SELs.all) * (1 - 2 * as.numeric(minus))
  if (!is.finite(SEL)) SEL <- bad.value
  if (attach.attributes) attr(SEL, "SELs") <- SELs.all

  return(SEL)
}

smoothEmplikMixed <- function(rho, theta, data,
                              by = NULL,
                              sel.weights = NULL,
                              parallel.in = FALSE,
                              parallel.out = FALSE,
                              cores = 1,
                              trim = NULL,
                              bad.value = -Inf,
                              minus = FALSE,
                              weight.tolerance = NULL,
                              ...) {
  if (is.null(by)) stop("You forgot to supply the factor or integer variable 'by' indicating the unique values of the conditioning set.")
  if (parallel.in & parallel.out) {
    warning("Cannot parallelise at both levels, setting inner parallelisation to FALSE.")
    parallel.in <- FALSE
  }


  rho.series <- rho(theta, data, ...)
  n <- length(rho.series)
  if (is.null(weight.tolerance)) weight.tolerance <- 0.01 / n
  rho.list <- split(rho.series, f = by)
  if (is.null(trim)) trim <- rep(1, n)
  trim.list <- split(trim, f = by)
  if (!is.list(sel.weights)) stop("The current implementation only allows lists of weights.")

  SEL.block <- function(i) {
    x <- rho.list[[i]]
    w <- sel.weights[[i]]
    SEL.b.j <- function(j) {
      xj <- x[w[[j]]$idx]
      wj <- w[[j]]$ct
      if (all(is.na(xj))) return(list(logelr = 0, lam = 0, wts = rep(1 / length(xj), length(xj)), converged = TRUE, iter = 0, bracket = c(0, 0), estim.prec = NA, f.root = NA, exitcode = 0))
      suppressWarnings(weightedEL(z = xj, ct = wj, SEL = TRUE, weight.tolerance = weight.tolerance))
    }
    empliklist <- if (parallel.in) parallel::mclapply(1:length(x), SEL.b.j, mc.cores = cores) else lapply(1:length(x), SEL.b.j)
    logsemplik <- sum(trim.list[[i]] * unlist(lapply(empliklist, "[[", "logelr")))
    if (!is.finite(logsemplik)) logsemplik <- bad.value
    attr(logsemplik, "SELR") <- unlist(lapply(empliklist, "[[", "logelr"))
    return(logsemplik)
  }

  SEL <- if (parallel.out) parallel::mclapply(1:length(rho.list), SEL.block, mc.cores = cores) else lapply(1:length(rho.list), SEL.block)
  SSEL <- sum(unlist(SEL)) * (1 - 2 * as.numeric(minus))
  attr(SSEL, "SELs") <- unlist(lapply(SEL, attr, which = "SELR"))

  return(SSEL)
}



