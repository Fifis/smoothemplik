#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.13 (2025-02-05).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}
