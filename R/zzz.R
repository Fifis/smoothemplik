#' @importFrom Rdpack reprompt

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Smoothed Empirical Likelihood estimation version 0.0.7 (2023-03-29).")
  packageStartupMessage("This is *NOT* a stable version. Core function subject to change.")
}

