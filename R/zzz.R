#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.12 (2024-06-13).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}
