% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scelcount.R
\name{cemplik}
\alias{cemplik}
\title{Self-concordant empirical likelihood for a vector mean with counts}
\source{
This original code was written for \insertCite{owen2013self}{smoothemplik}
and [published online](https://artowen.su.domains/empirical/) by Art B. Owen
(February 2014) and slightly reworked to contain fewer inner functions and loops.
}
\usage{
cemplik(
  z,
  ct = NULL,
  mu = NULL,
  lam = NULL,
  eps = NULL,
  M = Inf,
  order = 4,
  wttol = 0.001,
  thresh = 1e-30,
  itermax = 100,
  ALPHA = 0.3,
  BETA = 0.8,
  BACKEPS = 0,
  verbose = FALSE
)
}
\arguments{
\item{z}{A numeric vector or a matrix with one data vector per column.}

\item{ct}{A numeric vector of non-negative counts.}

\item{mu}{Hypothesised mean, default (0 ... 0) in R^ncol(z)}

\item{lam}{Starting lambda, default (0 ... 0)}

\item{eps}{Lower cut-off for \code{mllog()}, default 1/nrow(z)}

\item{M}{Upper cutoff for \code{mllog()}, default Inf}

\item{order}{Positive integer such that the Taylor approximation of this order to log(x) is self-concordant; usually 4 or higher. Passed to \code{mllog}.}

\item{wttol}{Weight tolerance for counts to improve numerical stability}

\item{thresh}{Convergence threshold for log-likelihood (the default is aggressive)}

\item{itermax}{Upper bound on number of Newton steps (seems ample)}

\item{ALPHA}{Backtracking line search parameter: acceptance of a decrease in function value by ALPHA*f of the prediction
based on the linear extrapolation.}

\item{BETA}{Backtracking line search reduction factor. 0.1 corresponds to a very crude search, 0.8 corresponds
to a less crude search.}

\item{BACKEPS}{Backtrack threshold: the search can miss by this much. Consider setting it to 1e-10
if backtracking seems to be failing due to round-off.}

\item{verbose}{Logical: print output diagnostics?}
}
\value{
A list with the following values:
\describe{
    \item{logelr}{Log of empirical likelihood ratio (equal to 0 if the hypothesised mean is equal to the sample mean)}
    \item{lam}{Vector of Lagrange multipliers}
    \item{wts}{Observation weights/probabilities (vector of length n)}
    \item{converged}{\code{TRUE} if algorithm converged. \code{FALSE} usually means that mu is not in the convex hull of the data. Then, a very small likelihood is returned (instead of zero).}
    \item{iter}{Number of iterations taken.}
    \item{ndec}{Newton decrement (see Boyd & Vandenberghe).}
    \item{gradnorm}{Norm of the gradient of log empirical likelihood.}
}
}
\description{
Implements the empirical-likelihood-ratio test for the mean of the coordinates of \code{z}
(with the hypothesised value \code{mu}). The counts need not be integer;
in the context of local likelihoods, they can be kernel observation weights.
}
\details{
Negative weights are not allowed. They could be useful in some applications, but they can destroy
convexity or even boundedness. They also make the Newton step fail to be of least squares type.

This function relies on the improved computational strategy for the empirical likelihood.
The search of the lambda multipliers is carried out via a dampened Newton method with guaranteed
convergence owing to the fact that the log-likelihood is replaced by its Taylor approximation
of any desired order (default: 4, the minimum value that ensures self-concordance).

Tweak \code{ALPHA} and \code{BETA} with extreme caution. See \insertCite{boyd2004convex}{smoothemplik},
pp. 464--466 for details. It is necessary that \code{0 < ALPHA < 1/2} and \code{0 < BETA < 1}.
\code{ALPHA = 0.03} seems better than 0.01 on some 2-dimensional test data (sometimes fewer iterations).

The argument names are matching the original names in Art B. Owen's implementation.
The highly optimised one-dimensional counterpart, \code{weightedEL}
}
\examples{
earth <- c(
  5.5, 5.61, 4.88, 5.07, 5.26, 5.55, 5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29,
  5.44, 5.34, 5.79, 5.1, 5.27, 5.39, 5.42, 5.47, 5.63, 5.34, 5.46, 5.3, 5.75, 5.68, 5.85
)
cemplik(earth, mu = 5.517, verbose = TRUE) # 5.517 is the modern accepted value

# Linear regression through empirical likelihood
coef.lm <- coef(lm(mpg ~ hp + am, data = mtcars))
xmat <- cbind(1, as.matrix(mtcars[, c("hp", "am")]))
yvec <- mtcars$mpg
foc.lm <- function(par, x, y) {  # The sample average of this
  resid <- y - drop(x \%*\% par)   # must be 0
  resid * x
}
minusEL <- function(par) -cemplik(foc.lm(par, xmat, yvec), itermax = 10)$logelr
coef.el <- optim(c(mean(yvec), 0, 0), minusEL)$par
abs(coef.el - coef.lm) / coef.lm  # Relative difference

# Likelihood ratio testing without any variance estimation
# Define the profile empirical likelihood for the coefficient on am
minusPEL <- function(par.free, par.am)
  -cemplik(foc.lm(c(par.free, par.am), xmat, yvec), itermax = 20)$logelr
# Constrained maximisation assuming that the coef on par.am is 3.14
coef.el.constr <- optim(coef.el[1:2], minusPEL, par.am = 3.14)$par
print(-2 * cemplik(foc.lm(c(coef.el.constr, 3.14), xmat, yvec))$logelr)
# Exceeds the critical value qchisq(0.95, df = 1)
}
\references{
\insertAllCited{}
}
