# smoothemplik 0.0.11 (2023-MM-DD)

- Added EUPL licence

# smoothemplik 0.0.10 (2023-12-05)

- Added support for weighted kernel density and regression observation
- Added support for sparse weight matrices via `sparseKernelWeightsCPP` to save memory
- Added observation de-duplication to speed up the non-parametric functions
- Added low-level C++ parallelisation and chunking to limit maximum RAM use (via `RcppParallel`)
- Added leave-one-out kernel smoothing support for custom output grids
- Reworked mixed kernel density and regression workflow making use of the block data structure
- Fixed the bug that eliminated time savings in multicore mixed estimation
- Fixed the bug in DCV where duplicated observations were not merged (now the LOO is the TRUE LOO: duplicates of the same observation are also left out)
- Improved the initial value choice for cross-validation (using a log-scale grid around the rule-of-thumb value)
- Sped up C++ kernel functions with `RcppArmadillo` (20--50% speed gains through better iterations, code structure, and condition checks)

# smoothemplik 0.0.9 (2023-09-08)

- Added a vignette on non-parametric methods
- Added 4th-order C++ versions of all kernels (for bias reduction) and their convolutions
- Auto-detect cross-validation type (density or least-squares) based on the input
- Changed the default behaviour of Silverman's rule of thumb: use a robust estimator of SD (IQR/1.34)
- Prepared a stub for discontinuous densities
- Moving the gradient-related functions to a new package, `pnd`

# smoothemplik 0.0.8 (2023-06-03)

- Rewrote the C++ smoothers using `RcppArmadillo` for speed-up, refactored the kernel-related code
- Support for the 2nd-order uniform, triangular, Epanechnikov, and quartic kernel

# smoothemplik 0.0.7 (2023-03-29)

- Added functions for parallelised numerical differentiation
- Rewrote multi-variate weighted empirical likelihood functions to allow for Taylor approximations of the empirical likelihood function of any order.

# smoothemplik 0.0.6

- Feature: getSELWeights now renormalises the weights to unity after trimming (by default, can be overridden via `renormalise = FALSE`).

# smoothemplik 0.0.5 (2021-10-12)

- Initial release.

