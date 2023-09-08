#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.9 (2023-09-08).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}

