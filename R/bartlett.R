#' Bartlett correction factor for empirical likelihood with estimating equations
#'
#' Compute the Bartlett correction factor \eqn{b} for empirical likelihood based on
#' the moment conditions \eqn{E\{g(X;\theta_0)\}=0}. The function implements the rotation
#' in \insertCite{liu2010adjusted}{smoothemplik} and evaluates \eqn{b} either from raw moments (unadjusted) or from
#' the bias-reduced moment estimators recommended in their paper.
#'
#' @param x Numeric vector or matrix of estimating functions. If a matrix,
#'   rows are observations and columns are the components of \eqn{g}.
#' @param centre Logical. If `TRUE` (default), centre each column of `x` by its
#'   sample mean before computing the correction (this corresponds to plugging in
#'   a consistent \eqn{\hat\theta} so that
#'   \eqn{n^{-1}\sum g_i(\hat\theta)\approx 0}{1/n sum_i g_i(^theta) ~ 0}).
#' @param bias.adj Logical. If `TRUE` (default), use the bias-reduced
#'   moment estimators. When \eqn{n \le 4}, the adjustment is disabled automatically.
#'
#' @returns Numeric scalar: the estimated Bartlett correction factor \eqn{b}.
#'   For multivariate inputs, the value has an attribute \code{"components"} equal to
#'   \code{c(b1, b2)} where \code{b = b_1 - b_2}. If \code{bias.adj = TRUE}, attributes
#'   \code{"unadjusted"} and \code{"unadjusted.components"} store the corresponding
#'   unadjusted estimates.
#'
#' @details
#' Let \eqn{V(\theta) = \mathrm{Var}\{g(X, \theta)\}}{V(theta) = Var(g(X, theta))},
#' and let \eqn{P}{P} be the orthogonal matrix of eigenvectors of \eqn{V(\hat\theta)}{V(^theta)}.
#' Define the rotated variables \eqn{Y_i = g_i(\hat\theta) P} (observations in rows), and write
#' \eqn{\alpha^{rs\cdots t} = E(Y^r Y^s \cdots Y^t)}{ars..t = E(Y^r Y^s ... Y^t)} with
#' \eqn{\alpha^{rr}=E(Y_r^2)}{arr = E(Y[, r]^2)}.
#'
#' The Bartlett factor (Theorem 1 of \insertCite{liu2010adjusted}{smoothemplik}) can be written compactly as
#' \deqn{
#' b = \frac{1}{q}\Bigl\{
#'   \frac{1}{2} \sum_{r,s} \frac{\alpha^{rrss}}{\alpha^{rr} \alpha^{ss}}
#'   -
#'   \frac{1}{3} \sum_{r,s,t} \frac{(\alpha^{rst})^2}{\alpha^{rr} \alpha^{ss} \alpha^{tt}}
#' \Bigr\},
#' }
#' where \eqn{q} is the dimension of \eqn{g}. The first double sum is over all
#' pairs \code{(r, s)}, and the triple sum is over all triples \code{(r, s, t)}.
#'
#' For adjusted-EL applications, the implementation also uses the equivalent
#' decomposition \code{b = b_1 - b_2}.
#'
#' When \code{bias.adj = TRUE}, all moments are replaced by the
#' bias-reduced estimators given in Eq. (10) and the table beneath it in
#' \insertCite{liu2010adjusted}{smoothemplik}.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @examples
#' set.seed(1)
#'
#' # One-dimensional: Bartlett factor for the mean
#' x <- rchisq(50, df = 4)
#' bartlettFactor(x)  # Bias-adjusted
#' bartlettFactor(x, bias.adj=FALSE)
#'
#' # Multi-variate g(X; theta): columns are components of g
#' n <- 100
#' g <- cbind(rchisq(n, 4)-4, rchisq(n, 3)-3, rchisq(n, 6)-6, rnorm(n))
#' bartlettFactor(g)  # Bias-adjusted, centred
#' bartlettFactor(g, centre = FALSE)  # The true average was used in g
bartlettFactor <- function(x, centre = TRUE, bias.adj = TRUE) {
  n <- NROW(x)
  q <- NCOL(x)
  if (is.null(dim(x))) x <- as.matrix(x)
  if (centre) x <- sweep(x, 2, colMeans(x), "-")

  if (bias.adj && n <= 4) {
    warning("Bias adjustment of the Bartlett correction requires n > 4. Switching to unadjusted.")
    bias.adj <- FALSE
  }

  v <- stats::var(x) * (n-1) / n
  e <- eigen(v)
  P <- e$vectors
  xi <- e$values
  # v - P %*% diag(xi) %*% t(P)  # The decomposition is right
  # P %*% t(P)  # Also right
  Y <- x %*% P
  # colMeans(Y^2) / xi  # Must be 1!

  if (q == 1) {
    if (!is.null(dim(Y))) Y <- as.numeric(Y[, 1])
    a2 <- mean(Y^2)
    a3 <- mean(Y^3)
    a4 <- mean(Y^4)
    b0 <- a4 / a2^2 / 2 - a3^2 / a2^3 / 3
    if (bias.adj) {
      a6 <- mean(Y^6)
      a2t <- a2 * n / (n-1)
      a4t <- (n*a4 - 6*a2t^2) / (n-4)
      a2_2t <- a2t^2 - a4t/n
      a3t <- n*a3 / (n-3)
      a3_2t <- a3t^2 - (a6 - a3t^2)/n
      a2_3t <- a2t^3
      b1 <- a4t / a2_2t / 2 - a3_2t / a2_3t / 3
      attr(b1, "unadjusted") <- b0
      b0 <- b1
    }
    return(b0)
  } else {
    a1 <- a2 <- a3 <- a4 <- a5 <- a6 <- 0
    a1t <- a2t <- a3t <- a4t <- a5t <- a6t <- 0
    arr.vec <- sapply(1:q, function(r) mean(Y[, r]^2))
    if (bias.adj) arr.t.vec <- arr.vec * n / (n - 1)

    # Computing a1, a2
    for (r in 1:q) {  # Single-index loop
      arr <- arr.vec[r]
      arrr <- mean(Y[, r]^3)
      arrrr <- mean(Y[, r]^4)
      a1 <- a1 + arrrr / (arr^2)
      a2 <- a2 + arrr^2 / (arr^3)
      if (bias.adj) {
        arr.t <- arr.t.vec[r]
        arrr.t <- arrr * n / (n-3)
        arrrr.t <- (n*arrrr - 6*arr.t^2) / (n-4)
        a1t <- a1t + arrrr.t / (arr.t^2)
        a2t <- a2t + arrr.t^2 / (arr.t^3)
      }
    }

    # Double loop for indices r != s
    # Computing a3, a4, a5 = a4
    for (r in 1:q) {
      arr <- arr.vec[r]
      if (bias.adj) arr.t <- arr.t.vec[r]
      for (s in setdiff(1:q, r)) {
        ass <- arr.vec[s]
        arrss <- mean(Y[, r]^2 * Y[, s]^2)
        arss <- mean(Y[, r] * Y[, s]^2)
        arrs <- mean(Y[, r]^2 * Y[, s])
        arrssss <- mean(Y[, r]^2 * Y[, s]^4)
        arrrrss <- mean(Y[, r]^4 * Y[, s]^2)
        a3 <- a3 + arrss / (arr*ass)
        a4 <- a4 + arss^2 / (arr * ass^2)
        a5 <- a5 + arrs^2 / (arr^2 * ass)
        if (bias.adj) {
          ass.t <- arr.t.vec[s]
          arrss.t <- (n*arrss - 2*arr.t*ass.t) / (n-4)
          arrass.t <- arr.t*ass.t - arrss.t/n
          arss.t <- arss * n / (n-3)
          arrs.t <- arrs * n / (n-3)
          arss_2.t <- arss.t^2 - (arrssss - arss.t^2)/n
          arrs_2.t <- arrs.t^2 - (arrrrss - arrs.t^2)/n
          a3t <- a3t + arrss.t / (arrass.t)
          a4t <- a4t + arss_2.t / (arr.t * ass.t^2)
          a5t <- a5t + arrs_2.t / (arr.t^2 * ass.t)
        }
      }
    }

    if (q >= 3) {
      for (r in 1:(q-2)) {
        arr <- arr.vec[r]
        if (bias.adj) arr.t <- arr.t.vec[r]
        for (s in (r+1):(q-1)) {
          ass <- arr.vec[s]
          if (bias.adj) ass.t <- arr.t.vec[s]
          for (t in (s+1):q) {
            att <- arr.vec[t]
            arst <- mean(Y[, r] * Y[, s] * Y[, t])
            a6 <- a6 + arst^2 / (arr * ass * att)
            if (bias.adj) {
              att.t <- arr.t.vec[t]
              arst.t <- arst * n / (n-3)
              arrsstt <- mean(Y[, r]^2 * Y[, s]^2 * Y[, t]^2)
              arst_2.t <- arst.t^2 - (arrsstt - arst.t^2) / n
              a6t <- a6t + arst_2.t / (arr.t * ass.t * att.t)
            }
          }
        }
      }

      # # Debug method: computing unadjusted b in 2 components over all indices
      # # balt MUST coincide with bd without bias adjustments
      # s1 <- s2 <- 0
      # for (r in 1:q) {
      #   arr <- arr.vec[r]
      #   for (s in 1:q) {
      #     ass <- arr.vec[s]
      #     arrss <- mean(Y[, r]^2 * Y[, s]^2)
      #     s1 <- s1 + arrss / (arr*ass)
      #   }
      # }
      # for (r in 1:q) {
      #   arr <- arr.vec[r]
      #   for (s in 1:q) {
      #     ass <- arr.vec[s]
      #     for (t in 1:q) {
      #       att <- arr.vec[t]
      #       arst <- mean(Y[, r]*Y[, s]*Y[, t])
      #       s2 <- s2 + arst^2 / (arr * ass * att)
      #     }
      #   }
      # }
      # balt <- 1/q * (s1/2 - s2/3)

    }  # End if q>=3

    b1 <- (a1/2 - a2/3 + a3/2 - a4/2) / q
    b2 <- (a5/2 + 2*a6) / q
    if (!is.finite(b1) || !is.finite(b2))
      stop("Bartlett factor cannot be computed for rank-deficient data sets with linearly dependent columns.")
    if (b1 < 0 || b2 < 0) stop("Bug in bartlettFactor: the sums b1 and b2 cannot be negative. Please report it to GitHub.")
    bd <- b1 - b2
    attr(bd, "components") <- c(b1, b2)
    if (bias.adj) {
      b1t <- (a1t/2 - a2t/3 + a3t/2 - a4t/2) / q
      b2t <- (a5t/2 + 2*a6t) / q
      if (b1t < 0 || b2t < 0) {
        warning("The Bartlett sums b1 and b2 cannot be negative; removing the adjustment.")
        b1t <- b1
        b2t <- b2
      }
      bdt <- b1t - b2t
      attr(bdt, "components") <- c(b1t, b2t)
      attr(bdt, "unadjusted") <- as.numeric(bd)
      attr(bdt, "unadjusted.components") <- c(b1, b2)
      bd <- bdt
    }
    return(bd)
    # return(list(bd, balt = balt))  # For debugging with the debug block
  }  # End if d > 1
}
