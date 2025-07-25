% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{tlog}
\alias{tlog}
\title{d-th derivative of the k-th-order Taylor expansion of log(x)}
\usage{
tlog(x, a = as.numeric(c(1)), k = 4L, d = 0L)
}
\arguments{
\item{x}{Numeric: a vector of points for which the logarithm is to be evaluated}

\item{a}{Scalar: the point at which the polynomial approximation is computed}

\item{k}{Non-negative integer: maximum polynomial order in the Taylor expansion
of the original function. \code{k = 0} returns a constant.}

\item{d}{Non-negative integer: derivative order

Note that this function returns the d-th derivative of the k-th-order Taylor expansion, not the
k-th-order approximation of the d-th derivative. Therefore, the degree of the resulting polynomial
is \eqn{d-k}{d-k}.}
}
\value{
The approximating Taylor polynomial around \code{a} of the order \code{d-k} evaluated at \code{x}.
}
\description{
d-th derivative of the k-th-order Taylor expansion of log(x)
}
\examples{
cl <- rainbow(9, end = 0.8, v = 0.8, alpha = 0.8)
a <- 1.5
x <- seq(a*2, a/2, length.out = 101)
f <- function(x, d = 0)  if (d == 0) log(x) else ((d\%\%2 == 1)*2-1) * 1/x^d * gamma(d)
oldpar <- par(mfrow = c(2, 3), mar = c(2, 2, 2.5, 0.2))
for (d in 0:5) {
y <- f(x, d = d)
plot(x, y, type = "l", lwd = 7, bty = "n", ylim = range(0, y),
       main = paste0("d^", d, "/dx^", d, " Taylor(Log(x))"))
  for (k in 0:8) lines(x, tlog(x, a = a, k = k, d = d), col = cl[k+1], lwd = 1.5)
  points(a, f(a, d = d), pch = 16, cex = 1.5, col = "white")
}
legend("topright", as.character(0:8), title = "Order", col = cl, lwd = 1)
par(oldpar)
}
