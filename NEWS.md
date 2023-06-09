# smoothemplik 0.0.9

- Added a vignette on non-parametric methods
- Added 4th-order C++ versions of all kernels (for bias reduction) and their convolutions
- Auto-detect cross-validation type (density or least-squares) based on the input
- Changed the default behaviour of Silverman's rule of thumb: use a robust estimator of SD (IQR/1.34)

# smoothemplik 0.0.8

- Rewrote the C++ smoothers using `RcppArmadillo` for speed-up, refactored the kernel-related code
- Support for the 2nd-order uniform, triangular, Epanechnikov, and quartic kernel

# smoothemplik 0.0.7

- Added functions for parallelised numerical differentiation
- Rewrote multi-variate weighted empirical likelihood functions to allow for Taylor approximations of the empirical likelihood function of any order.

# smoothemplik 0.0.6

- Feature: getSELWeights now renormalises the weights to unity after trimming (by default, can be overridden via `renormalise = FALSE`).

# smoothemplik 0.0.5

- Initial release.

