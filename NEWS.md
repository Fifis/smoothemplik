# smoothemplik 0.next (2024-XX-XX)

These bug fixes and features are scheduled for the upcoming releases.

- Add a vignette for non-parametric methods to GitHub, finish the mixed-smoothing part
- If `x` contains duplicates, `DCV(x, bw = bw.grid, weights = w)` complains that no duplicates were found (see the example)
- Test the previous output of `weightedEL()` and `cemplik()` with the new version, add tests
- `maximiseSEL()` fails if no weights are provided -- check and auto-denerate
- De-duplicate at kernel weights already (via `.prepareKernel()`), return the attribute
- For `.prepareKernel()` AND mixed kernel: check if the max. column-wise gap between observations is >= than the bandwidth, otherwise write an informative message
- `bw.CV()` bug: `bw.CV(x, y = y, kernel = "quartic", order = 4)` and `bw.CV(x, y = y, kernel = "quartic")`
- LOO estimation: instead of dropping unique (X, Y) observations, leave each conditioning points (only X)
- Add weight support to `kernelDiscreteDensitySmooth()`
- Add `RcppParallel::setThreadOptions(numThreads = "auto")` as the 1st line of parallel-capable functions, use `setDTthreads` also
- Write test cases for C++ functions with and without speed-ups
- Eliminate matrices in smoothing completely, try only parallel loops
- Create a default value for `memsave` and when to invoke it (based on `nx*ng`)
- Remove parallelisation over workers via setThreadOptions when there is outer parallelisation in `.kernelMixed()`
- Make DCV either sparse or memsave, not both; reflect the changes in `bw.CV()`
- Fix the DCV code with convolutions (especially the quartic one)
- Move the de-duplication of the xout grid inside `kernelSmooth`
- Check all instances of `kernelSmooth()`, `kernelDensity()`, `kernelWeights()`, and everything that used the obsolete arguments
- The check `CV = "DCV"` and `is.null(y)` seems redundant.
- Speed up the computation of `kernelWeights` when `x = xgrid` by exploiting the symmetry of kernels.
- Check where `kernelWeights` and `standardise` were used.
- Fix the optimiser control argument in `bw.CV()`.
- Extend the CV to non-Gaussian cases.
- `smoothEmplikDiscrete()`: if the split variable does not take contiguous values from 1 to the number of categories, estimation fails.
- In `maximiseSEL()` and `constrSmoothEmplik()`, add names to `$restricted`
- CV: implement leave-K-out CV for speed
- In `kernelMixedSmooth()`: if LOO, do not de-duplicate `xout`, copy it from `arg$x` (currently mitigated via `deduplicate.xout = FALSE`)
- All LOO to the C++ density function
- Add custom kernels to Silverman's rule of thumb (with roughness != 1)
- Check: if the kernel is finite-support and bandwidth is smaller than the largest gap between two observations, then
- Like in the SEL application: de-duplicate the input matrix, replace with weights; allow the user to disable it
- Merging cells: allow arbitrary variables (including continuous ones) for proximity.
- `kernelSmooth()` and `kernelDensity()` should have an argument for increasing small bandwidths in case of zero weights to match the largest gap divided by 2 (times 1.1 to have at least some coverage)
- Check analytical expressions for all combinations of kernels, convolutions, and orders
- Create convolution for kernel orders 4 and 6;
- Check with R 3.0.0
- Reproduce the CKT (2019) results with the `shift` argument (i.e. test the shift)
- Create a summary class for SEL; print numerical gradients of lambdas; print the number of converged inner optimisation problems
- Address all remaining issues with with `todor`
- Check release with `lintr::lint_package()`, `R CMD check --as-cran`
- Add tests reproducing simple hard-coded examples


# smoothemplik 0.0.12 (2024-06-13)

- Fixed a bug in `prepareKernel()` where a valid `y` vector with attributes would not pass the check
- Implemented a more accurate check for lambda being close to the boundary based on the relative search interval length in `weightedEL()`
- `weightedEL()` preserves the names of the input vector in `wts`
- Sped up `ctracelr()` by using the previous lambda value in the search (~4 times faster)
- The output of `mllog()` now has column names because it was confusing without them
- The output of `svdlm()` is now a vector, not a 1-column matrix
- Replaced certain instances of `sapply()` with `vapply()` in smoothing functions
- Added unit tests for some functions


# smoothemplik 0.0.11 (2024-01-13)

- Added EUPL licence
- Initialised tests for unit testing
- Fixed the bug in DCV when weights were not taken into account
- Removed simulation functions for the paper (to be provided separately)
- Added examples for most functions


# smoothemplik 0.0.10 (2023-12-05)

- Feature: support for weighted kernel density and regression observation
- Feature: support for sparse weight matrices via `sparseKernelWeightsCPP()` to save memory
- Feature: observation de-duplication to speed up the non-parametric functions
- Feature: low-level C++ parallelisation and chunking to limit maximum RAM use (via `RcppParallel`)
- Feature: leave-one-out kernel smoothing support for custom output grids
- Reworked mixed kernel density and regression workflow making use of the block data structure
- Bug fix: now there are time savings in multi-core mixed estimation
- Bug fix: in DCV, duplicated observations were not merged (now the LOO is the true LOO: duplicates of the same observation are also left out).
- Improved the initial value choice for cross-validation (using a log-scale grid around the rule-of-thumb value).
- Sped up C++ kernel functions with `RcppArmadillo` (20--50% speed gains through better iterations, code structure, and condition checks).


# smoothemplik 0.0.9 (2023-09-08)

- Added a vignette on non-parametric methods.
- Added 4th-order C++ versions of all kernels (for bias reduction) and their convolutions.
- Auto-detect cross-validation type (density or least-squares) based on the input.
- Changed the default behaviour of Silverman's rule of thumb: use a robust estimator of SD (IQR/1.34).
- Prepared a stub for discontinuous densities.
- Moving the gradient-related functions to a new package, `pnd`.


# smoothemplik 0.0.8 (2023-06-03)

- Rewrote the C++ smoothers using `RcppArmadillo` for speed-up, refactored the kernel-related code.
- Feature: Support for the 2nd-order uniform, triangular, Epanechnikov, and quartic kernel.


# smoothemplik 0.0.7 (2023-03-29)

- Feature: added functions for parallelised numerical differentiation.
- Rewrote multi-variate weighted empirical likelihood functions to allow for Taylor approximations of the empirical likelihood function of any order.


# smoothemplik 0.0.6

- Feature: getSELWeights now renormalises the weights to unity after trimming (by default, can be overridden via `renormalise = FALSE`).


# smoothemplik 0.0.5 (2021-10-12)

- Initial release.

