% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smoothers.R
\name{kernelDiscreteDensitySmooth}
\alias{kernelDiscreteDensitySmooth}
\title{Density and/or kernel regression estimator with conditioning on discrete variables}
\usage{
kernelDiscreteDensitySmooth(x, y = NULL, compact = FALSE, fun = mean)
}
\arguments{
\item{x}{A vector or a matrix/data frame of discrete explanatory variables (exogenous).
Non-integer values are fine because the data are split into bins defined by interactions of these variables.}

\item{y}{Optional: a vector of dependent variable values.}

\item{compact}{Logical: return unique values instead of full data with repeated observations?}

\item{fun}{A function that computes a statistic of \code{y} inside every category defined by \code{x}.}
}
\value{
A list with \code{x}, density estimator (\code{fhat}) and, if \code{y} was provided, regression estimate.
}
\description{
Density and/or kernel regression estimator with conditioning on discrete variables
}
\examples{
set.seed(1)
x <- sort(rnorm(1000))
p <- 0.5*pnorm(x) + 0.25 # Propensity score
d <- as.numeric(runif(1000) < p)
# g = discrete version of x for binning
g <- as.numeric(as.character(cut(x, -4:4, labels = -4:3+0.5)))
dhat.x <- kernelSmooth(x = x, y = d, bw = 0.4, no.dedup = TRUE)
dhat.g <- kernelDiscreteDensitySmooth(x = g, y = d)
dhat.comp <- kernelDiscreteDensitySmooth(g, d, compact = TRUE)
plot(x, p, ylim = c(0, 1), bty = "n", type = "l", lty = 2)
points(x, dhat.x, col = "#00000044")
points(dhat.comp, col = 2, pch = 16, cex = 2)
lines(dhat.comp$x, dhat.comp$fhat, col = 4, pch = 16, lty = 3)
}
