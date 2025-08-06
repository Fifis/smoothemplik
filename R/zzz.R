#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood version 0.0.15 (2025-08-05).")
  packageStartupMessage("Package under active development; core functions subject to change.")
}
