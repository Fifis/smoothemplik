#' Adding leading or trailing zeros
#'
#' @param x A number or a numeric vector to format.
#' @param d The desired length of the string.
#'
#' @return A character vector.
#' @export
lead0 <- function(x, d) formatC(x, width = d, format = "d", flag = "0")

#' @rdname lead0
#' @export
trail0 <- function(x, d) sprintf(paste0("%.", d, "f"), x)

#' Converting matrices to LaTeX tables
#'
#' @param x A number vector or matrix
#' @param d Passed to trail0.
#'
#' @return No return, just prints the output.
#' @export
matrix2LaTeX <- function(x, d = 3) {
  if (is.null(dim(x))) x <- t(as.matrix(x))
  for (i in 1:nrow(x)) cat(gsub("-", "$-$", paste(trail0(x[i, ], d), collapse = " & ")), "\\\\\n")
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
