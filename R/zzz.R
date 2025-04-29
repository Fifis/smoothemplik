#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.14 (2025-04-30).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}
