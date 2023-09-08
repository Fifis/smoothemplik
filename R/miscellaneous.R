#' Adding leading or trailing zeros
#'
#' @param x A number or a numeric vector to format.
#' @param d The desired length of the string.
#'
#' @return A character vector.
.lead0 <- function(x, d) formatC(x, width = d, format = "d", flag = "0")

#' @rdname lead0
.trail0 <- function(x, d) sprintf(paste0("%.", d, "f"), x)

#' Converting matrices to LaTeX tables
#'
#' @param x A number vector or matrix
#' @param d Passed to trail0.
#'
#' @return No return, just prints the output.
#' @export
matrix2LaTeX <- function(x, d = 3) {
  if (is.null(dim(x))) x <- t(as.matrix(x))
  for (i in 1:nrow(x)) cat(gsub("-", "$-$", paste(.trail0(x[i, ], d), collapse = " & ")), "\\\\\n")
  return(invisible(NULL))
}


#' Bandwidth decrement rule.
#'
#' @param pow Convergence power (2 for regular estimators, 5 for non-parametric ones).
#'
#' @return A numeric vector of length 5.
#' @export
decFac <- function(pow) {
  r <- c(500, 1000, 2000, 4000, 8000)^(-1/pow)
  r <- r / r[1]
  return(r)
}

#' @param labels Passed to \code{text}.
#' @rdname linesWithHalo
#' @export
textWithHalo <- function(x, y, labels, nhalo = 32, col.halo = "#FFFFFF44", col.main = "#000000", hscale = 0.04, vscale = 0.04, ...) {
  offsets <- cbind(cos(seq(0, 2*pi, length.out = nhalo + 1)) * hscale, sin(seq(0, 2*pi, length.out = nhalo + 1)) * vscale)[-(nhalo + 1), ]
  for (i in 1:nhalo) graphics::text(x + offsets[i, 1], y + offsets[i, 2], labels = labels, col = col.halo, ...)
  graphics::text(x, y, labels = labels, col = col.main, ...)
}

#' Adding text or lines with an outline for better legibility
#'
#' @param x Passed to the plotting command.
#' @param y Passed to the plotting command.
#' @param nhalo The number of text duplicates to arrange in a circle to create the white glow effect.
#' @param col.halo Halo colour.
#' @param col.main Text colour.
#' @param hscale Halo horizontal radius.
#' @param vscale Halo vertical radius.
#' @param ... Passed to the plotting command.
#'
#' @export
linesWithHalo <- function(x, y, nhalo = 32, hscale = 0.01, vscale = 0.01, col.halo = "#FFFFFF88", col.main = "#000000", ...) {
  offsets <- cbind(cos(seq(0, 2*pi, length.out = nhalo + 1)) * hscale, sin(seq(0, 2*pi, length.out = nhalo + 1)) * vscale)[-(nhalo + 1), ]
  for (i in 1:nhalo) graphics::lines(x + offsets[i, 1], y + offsets[i, 2], col = col.halo, ...)
  graphics::lines(x, y, col = col.main, ...)
}

#' Fail gracefully
#'
#' Return parseable output in case the SEL optimiser fails
#'
#' @return A list of the same format as the output of \code{\link{maximiseSEL}}.
.fail <- function() return(list(par = c(NA, NA), value = NA, restricted = c(FALSE, FALSE), code = NA, iterations = NA, xtimes = c(0, 0)))


#' Finite-difference coefficients for arbitrary grids
#'
#' This function computes the coefficients that yield an numerical approximation
#' or arbitrary order to the true first derivative of any function, or, for a given
#' grid (stencil) \emph{a}, the weights yielding the derivative approximation using the
#' function evaluated on the grid:
#' \deqn{\frac{df}{dx} \approx h^{-1} \sum_{i=1}^n a_i f(x + s_i\cdot h)}{f'(x) ~ sum_i a[i] f(x + s[i]*h)}
#'
#' @param d Order of the derivative
#' @param order Order of accuracy in terms of the step size (usually
#' \eqn{h^o/o! \cdot f^{(o)}(x)}{h^o/o! * d^o/dx^o f(x)}).
#' @param side Centred or one-sided differences. Unless the function is computationally prohibitively expensive, two-sided differences are strongly recommended.
#' @param stencil In case the user desires a non-uniform grid or a custom grid, a vector of points at which the function is to be evaluated.
#'
#' @details
#' The finite-difference coefficients for any given stencil are given as a solution of a linear equation
#' system. There derivation of the system is due to \insertCite{taylor2016finite}{smoothemplik}, although a similar
#' approach is described in \insertCite{fornberg1988generation}{smoothemplik}. This function reproduces the tables
#' from the latter paper exactly.

#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' fdCoef() # Simple two-sided derivative
#' fdCoef(4)$weights # Should be (1/12, -2/3, 2/3, -1/12)
#' # Replicating Table 1 from Fornberg (1988) "Generation of Finite Difference Formulas..."
#' t1 <- t(sapply((1:4)*2, function(o) c(rep(0, 4-o/2), fdCoef(o)$weights, rep(0, 4-o/2))))
#' colnames(t1) <- (setdiff(-4:4, 0)); rownames(t1) <- (1:4)*2
#' print(t1, digits = 4)
#' # Replicating Table 2 (ibid) for a halfway stencil
#' s2 <- c(-7, -5, -3, -1, 1, 3, 5, 7)/2
#' t2 <- t(sapply(1:4, function(o) c(rep(0, 4-o), fdCoef(stencil = s2[(5-o):(o+4)])$weights,
#'                                   rep(0, 4-o))))
#' colnames(t2) <- s2; rownames(t2) <- (1:4)*2
#' print(t2, digits = 4)
#' # Replicating Table 3 (ibid)
#' t3 <- t(sapply(1:8, function(o) c(fdCoef(stencil = 0:o)$weights, rep(0, 8-o))))
#' colnames(t3) <- 0:8; rownames(t3) <- 1:8
#' print(t3, digits = 4)
#' @export
fdCoef <- function(d = 1, order = 2, side = c("central", "forward", "backward"), stencil = NULL) {
  side <- side[1]
  if (length(order) != 1) stop("The 'order' argument must have length 1.")
  if (order < 1 | (side == "central" & order/2 != round(order/2)) | (order != round(order)))
    stop("The order of accuracy must be a positive integer. For 2-sided derivatives, it must be even.")
  if (is.null(stencil)) stencil <- switch(side, backward = (-order):0,
                                          central = -(order/2):(order/2),
                                          forward = 0:order)
  l <- length(stencil)
  A <- t(sapply(1:l, function(i) stencil^(i-1)))
  b <- numeric(l); b[1+d] <- factorial(d)
  weights <- solve(A, b)
  nzw <- abs(weights) > .Machine$double.eps
  stencil <- stencil[nzw]
  weights <- weights[nzw]
  return(list(stencil = stencil, weights = weights))
}

#' Parallelised gradient computation
#'
#' Computes a two- or one-sided numerical derivative that approximates the gradient | Jacobian using the indicated number of cores for maximum efficiency.
#'
#' @param func A function that returns a numeric scalar or a vector. If the function is vector-valued, the, the result is the Jacobian.
#' @param x A point at which the gradient or Jacobian needs to be estimated.
#' @param h The numerical difference step size. Too large = the slope of the secant is a bad estimator of the gradient, too small = ill conditioning (0/0).
#' @param order Desired order of accuracy. The error is usually O(h^order). To achieve this order of accuracy, the function needs to be evaluated \code{order*length(x)} times.
#' @param side Passed to \code{fdCoef()}. Centred or one-sided differences. Unless the function is computationally prohibitively expensive, two-sided differences are strongly recommended.
#' @param parallel If TRUE, estimates the gradient via finite differences where the function is evaluated in parallel.
#' @param use.mclapply If TRUE, uses a much quicker and simpler \code{parallel::mclapply} an process forking; if FALSE, requires a properly set up cluster.
#' @param cores Number of forked processes.
#' @param cluster A cluster on which the computations are done.
#' @param load.balance If TRUE, disables pre-scheduling for \code{mclapply} or enables load balancing via \code{parLapplyLB}.
#' @param ... Passed to `func`.
#'
#' Note that for one-sided problems, the step size that make the formula error
#' equal to the truncation error is of the order Mach.eps^(1/2) and for two-sided, Mach.eps^(1/3).
#' However, the optimal step size depends on the value of the higher-order derivatives
#' that is not available in general (or required extra computation that is, in turn, prone to numerical error).
#'
#' @return If \code{func} returns a scalar, a vector of the same length as \code{x}.
#' If \code{func} returns a vector, then, a matrix of dimensions \code{length(f(x)) length(x)}
#'
#' @examples
#' \dontrun{
#' slowFunScalar <- function(x) {Sys.sleep(0.04); print(x, digits = 12); sum(sin(x))}
#' slowFunVector <- function(x) {Sys.sleep(0.04); print(x, digits = 12); c(sum(sin(x)), sum(exp(x)))}
#' true.g <- cos(1:4) # Analytical gradient
#' true.j <- rbind(cos(1:4), exp(1:4)) # Analytical Jacobian
#' system.time(g.slow <- numDeriv::grad(slowFunScalar, x = 1:4) - true.g)
#' system.time(j.slow <- numDeriv::jacobian(slowFunVector, x = 1:4) - true.j)
#' system.time(g.fast <- gradParallel(slowFunScalar, x = 1:4,
#'                                    parallel = TRUE, use.mclapply = TRUE, cores = 4) - true.g)
#' system.time(j.fast <- gradParallel(slowFunVector, x = 1:4,
#'                                    parallel = TRUE, use.mclapply = TRUE, cores = 4) - true.j)
#' system.time(j.fast4 <- gradParallel(slowFunVector, x = 1:4, order = 4,
#'                                     parallel = TRUE, use.mclapply = TRUE, cores = 4) - true.j)
#' rownames(j.slow) <- c("numDeriv.jacobian", "")
#' rownames(j.fast) <- c("fast.jacobian.order2", "")
#' rownames(j.fast4) <- c("fast.jacobian.order4", "")
#' # Discrepancy
#' rbind(numDeriv.grad = g.slow, fast.grad = g.fast, j.slow, j.fast, j.fast4)
#' # The order-4 derivative is more accurate for functions with large high-order derivatives
#'}
#'
#' @export
gradParallel <- function(func, x, side = c("central", "forward", "backward"), order = 2,
                         h = if (side[1] == "central") .Machine$double.eps^(1/3) else sqrt(.Machine$double.eps),
                         parallel = FALSE, use.mclapply = FALSE, cores = 2, cluster = NULL, load.balance = TRUE,
                         ...) {
  side <- side[1]
  n <- length(x)
  s <- fdCoef(order = order, side = side)

  # Parallelising the task in the most efficient as possible: we need as many evaluations
  # as length(x) * order
  dx <- s$stencil*h
  par.grid <- expand.grid(arg = 1:n, grid = 1:length(dx))
  par.grid$step <- dx[par.grid$grid]
  par.grid$weights <- s$weights[par.grid$grid]
  xx <- lapply(1:nrow(par.grid), function(i) {
    newx <- x
    newx[par.grid$arg[i]] <- newx[par.grid$arg[i]] + par.grid$step[i]
    newx
  })
  ff <- if (parallel & use.mclapply & cores > 1) {
    parallel::mclapply(X = xx, FUN = function(x) func(x, ...), mc.cores = cores, mc.preschedule = !load.balance)
  } else if (parallel & !use.mclapply & !is.null(cluster)) {
    if (load.balance) parallel::parLapplyLB(cl = cluster, X = xx, fun = function(x) func(x, ...)) else
      parallel::parLapply(cl = cluster, X = xx, fun = function(x) func(x, ...))
  } else {
    lapply(xx, function(x) func(x, ...))
  }
  # The output can be vector-valued, which is why we work carefully with rows
  ff <- do.call(rbind, ff)
  ffw <- ff * par.grid$weights
  ffl <- split(as.data.frame(ffw), f = par.grid$arg)
  if (ncol(ff) == 1) ffl <- lapply(ffl, as.matrix) # The dimensions must be preserved
  ffs <- lapply(ffl, colSums) # A list of vectors or scalars
  jac <- unname(as.matrix(do.call(cbind, ffs)))
  if (nrow(jac) == 1) {
    jac <- as.numeric(jac)
    if (!is.null(names(x))) names(jac) <- names(x)
  } else {
    if (!is.null(names(x))) colnames(jac) <- names(x)
  }
  return(jac / h)
}

