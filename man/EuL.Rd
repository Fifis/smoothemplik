% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elfunctions.R
\name{EuL}
\alias{EuL}
\title{Multi-variate Euclidean likelihood with analytical solution}
\usage{
EuL(
  z,
  mu = NULL,
  ct = NULL,
  vt = NULL,
  shift = NULL,
  SEL = TRUE,
  weight.tolerance = NULL,
  trunc.to = 0,
  return.weights = FALSE,
  verbose = FALSE,
  chull.diag = FALSE
)
}
\arguments{
\item{z}{Numeric data vector.}

\item{mu}{Hypothesized mean of \code{z} in the moment condition.}

\item{ct}{Numeric count variable with non-negative values that indicates the multiplicity of observations.
Can be fractional. Very small counts below the threshold \code{weight.tolerance} are zeroed.}

\item{vt}{Numeric vector: non-negative variance weights for estimating the conditional
variance of \code{z}. Probabilities are returned only for the observations where \code{vt > 0}.}

\item{shift}{The value to add in the denominator (useful in case there are extra Lagrange multipliers): 1 + lambda'Z + shift.}

\item{SEL}{If \code{FALSE}, then the boundaries for the lambda search are based on the total sum of counts, like in vanilla empirical likelihood,
due to formula (2.9) in \insertCite{owen2001empirical}{smoothemplik}, otherwise according to Cosma et al. (2019, p. 170, the topmost formula).}

\item{weight.tolerance}{Weight tolerance for counts to improve numerical stability
(similar to the ones in Art B. Owen's 2017 code, but adapting to the sample size).}

\item{trunc.to}{Counts under \code{weight.tolerance} will be set to this value.
In most cases, setting this to \code{0} or \code{weight.tolerance} is a viable solution of the zero-denominator problem.}

\item{return.weights}{Logical: if TRUE, individual EL weights are computed and returned.
Setting this to FALSE gives huge memory savings in large data sets, especially when smoothing is used.}

\item{verbose}{Logical: if \code{TRUE}, prints warnings.}

\item{chull.diag}{Logical: if \code{TRUE}, checks if there is a definite convex hull failure
  in at least one dimension (\code{mu} being smaller than the smallest or larger
  than the largest element). Note that it does not check if \code{mu} is strictly in the
  convex hull because this procedure is much slower and is probably unnecessary.

The arguments \code{ct} and \code{vt} are responsible for smoothing of the moment function
and conditional variance, respectively. The objective function is
\deqn{\min_{p_{ij}} \frac1n \sum_{i=1}^n \sum_{j=1}^n \mathbb{I}_{ij} \frac{(p_{ij} -
  c_{ij})^2}{2v_{ij}}}{min_(p_ij) 1/n * sum_i sum_j I_ij (p_ij - c_ij)^2 / (2v_{ij})},
where \eqn{\mathbb{I}_{ij}}{I_ij} is 1 if \eqn{v_{ij} \ne 0}{v_ij != 0}.

This estimator is numerically equivalent to the Sieve Minimum Distance estimator
of \insertCite{ai2003efficient}{smoothemplik} with kernel sieves, but this interface
provides more flexibility through the two sets of weights. If \code{ct} and
\code{vt} are not provided, their default value is set to 1, and the resulting
estimator is the CUE-GMM estimator: a quadratic form in which the unconditional
mean vector is weighted by the inverse of the unconditional variance.}
}
\value{
A list with the same structure as that in [EL()].
}
\description{
Multi-variate Euclidean likelihood with analytical solution
}
\examples{
set.seed(1)
z <- cbind(rnorm(10), runif(10))
colMeans(z)
a <- EuL(z, return.weights = TRUE)
a$wts
sum(a$wts)  # Unity
colSums(a$wts * z)  # Zero
}
\references{
\insertAllCited{}
}
\seealso{
[EL()]
}
