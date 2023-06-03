#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.8 (2023-06-04).")
  packageStartupMessage("This is *NOT* a stable version. Core functions subject to change.")
}

