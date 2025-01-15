# smoothemplik 0.next (2024-XX-XX)

These bug fixes and features are scheduled for the upcoming releases.

- BUG: `bw.CV(x, y = y, kernel = "quartic", order = 4)` and `bw.CV(x, y = y, kernel = "quartic")`
- BUG: `maximiseSEL()` fails if no weights are provided -- check and auto-generate
- BUG: Fix the DCV code with convolutions (especially the quartic one)
- BUG: `smoothEmplik()` with `attach.attributes = TRUE` should be converted to `"all"`
- BUG: The check `CV = "DCV"` and `is.null(y)` seems redundant.
- BUG: If `x` contains duplicates, `DCV(x, bw = bw.grid, weights = w)` complains that no duplicates were found (see the example)
- BUG: `smoothEmplikDiscrete()`: if the split variable does not take contiguous values from 1 to the number of categories, estimation fails.
- BUG: LOO estimation: instead of dropping unique (X, Y) observations, leave each conditioning points (only X)
- BUG: Fix the optimiser control argument in `bw.CV()`.
- SYNTAX: `kernelSmooth()`, being a local average, should have `na.rm` and check the inputs
- SYNTAX: In `kernelDiscreteDensitySmooth()`, remove the `table` attribute and change the test.
- SYNTAX: Create a summary class for SEL; print numerical gradients of lambdas; print the number of converged inner optimisation problems
- FEATURE: For sparseVectorToList, the default `trim(x)` should be such that the sum of sorted weights exceeds 0.99999999: `trim = \(w) min(which(cumsum(sort(w / sum(w), decreasing = TRUE)) > 1 - 1e-8))`
- FEATURE: Create convolution for kernel orders 4 and 6
- FEATURE: De-duplicate at kernel weights already (via `.prepareKernel()`), return the attribute
- FEATURE: For `.prepareKernel()` AND mixed kernel: check if the max. column-wise gap between observations is >= than the bandwidth, otherwise write an informative message
- FEATURE: Make DCV either sparse or memsave, not both; reflect the changes in `bw.CV()`
- FEATURE: Speed up the computation of `kernelWeights` when `x = xgrid` by exploiting the symmetry of kernels.
- FEATURE: Remove parallelisation over workers via setThreadOptions when there is outer parallelisation in `.kernelMixed()`
- FEATURE: Move the de-duplication of the xout grid inside `kernelSmooth`
- FEATURE: Create a default value for `memsave` and when to invoke it (based on `nx*ng`)
- FEATURE: Add weight support to `kernelDiscreteDensitySmooth()`
- FEATURE: Eliminate matrices in smoothing completely, try only parallel loops
- FEATURE: Extend the CV to non-Gaussian cases.
- SYNTAX: In `maximiseSEL()` and `constrSmoothEmplik()`, add names to `$restricted`
- FEATURE: CV: implement leave-K-out CV for speed
- FEATURE: In `kernelMixedSmooth()`: if LOO, do not de-duplicate `xout`, copy it from `arg$x` (currently mitigated via `deduplicate.xout = FALSE`)
- FEATURE: All LOO to the C++ density function
- FEATURE: Add custom kernels to Silverman's rule of thumb (with roughness != 1)
- FEATURE: Check: if the kernel is finite-support and bandwidth is smaller than the largest gap between two observations, then, set the bandwidth in that dimension to 1.1 times that gap. `kernelSmooth()` and `kernelDensity()` should have an argument for increasing small bandwidths in case of zero weights to match the largest gap divided by 2 (times 1.1 to have at least some coverage)
- FEATURE: Like in the SEL application: de-duplicate the input matrix, replace with weights; allow the user to disable it
- FEATURE: Merging cells: allow arbitrary variables (including continuous ones) for proximity.
- MISC: Test the previous output of `weightedEL()` and `cemplik()` with the new version, reproduce Owen's examples, add tests
- MISC: Check analytical expressions for all combinations of kernels, convolutions, and orders in Sage and Mathematica, publish the reproducing codes
- MISC: Reproduce the CKT (2019) results with the `shift` argument (i.e. test the shift)
- MISC: Add a vignette for non-parametric methods to GitHub, finish the mixed-smoothing part
- DEV: Check all instances of `kernelSmooth()`, `kernelDensity()`, `kernelWeights()`, and everything that used obsolete arguments in the examples
- DEV: Too much CPU time in the examples for `kernelDensity` and `kernelSmooth` (>5 s) (make dontrun?)
- DEV: Add `RcppParallel::setThreadOptions(numThreads = "auto")` as the 1st line of parallel-capable functions, use `setDTthreads` also (check how `data.table` does it)
- DEV: Write test cases for C++ functions with and without speed-ups
- DEV: Check compatibility with R 3.0.0
- DEV: Add tests reproducing simple hard-coded examples
- DEV: Check the release with `todor::todor_package()`, `lintr::lint_package()`, `R CMD check --as-cran`, and `goodpractice::gp()`


# smoothemplik 0.0.12 (2024-06-13)

- Fixed a bug in `prepareKernel()` where a valid `y` vector with attributes would not pass the check.
- Implemented a more accurate check for lambda being close to the boundary based on the relative search interval length in `weightedEL()`.
- `weightedEL()` preserves the names of the input vector in `wts`.
- Sped up `ctracelr()` by using the previous lambda value in the search (~4 times faster).
- The output of `mllog()` now has column names because it was confusing without them.
- The output of `svdlm()` is now a vector, not a 1-column matrix.
- Replaced certain instances of `sapply()` with `vapply()` in smoothing functions.
- Added unit tests for some functions.


# smoothemplik 0.0.11 (2024-01-13)

- Added EUPL licence.
- Initialised tests for unit testing.
- Fixed the bug in DCV when weights were not taken into account.
- Removed simulation functions for the paper (to be provided separately).
- Added examples for most functions.


# smoothemplik 0.0.10 (2023-12-05)

- Feature: support for weighted kernel density and regression observation.
- Feature: support for sparse weight matrices via `sparseKernelWeightsCPP()` to save memory.
- Feature: observation de-duplication to speed up the non-parametric functions.
- Feature: low-level C++ parallelisation and chunking to limit maximum RAM use (via `RcppParallel`).
- Feature: leave-one-out kernel smoothing support for custom output grids.
- Reworked mixed kernel density and regression workflow making use of the block data structure.
- Bug fix: now there are time savings in multi-core mixed estimation.
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

