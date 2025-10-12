# From the manual: Inserting references in Rd and roxygen2 documentation (G. N. Boshnakov)
# Footnote 1: Any function for package Rdpack will do. This is to avoid getting a warning from `R CMD check`.
#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood version 0.0.16 (2025-11-11).\nPackage under active development; core functions subject to change.")
}
