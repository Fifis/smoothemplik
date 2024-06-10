<!-- badges: start -->
[![R-CMD-check](https://github.com/Fifis/smoothemplik/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Fifis/smoothemplik/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# smoothemplik

R package for estimation via Smoothed Empirical Likelihood.

This package provides direct functionality for the following articles:

* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2019). Inference in conditional moment restriction models when there is selection due to stratification. In *The Econometrics of Complex Survey Data* (Vol. 39, pp. 137–171). Emerald Publishing Limited. https://doi.org/10.1108/S0731-905320190000039010
* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2024). Missing endogenous variables in conditional moment restriction models. *In progress.*

The theory behind the method is provided in the following articles:

* Tripathi, G., & Kitamura, Y. (2003). Testing conditional moment restrictions. *The Annals of Statistics, 31(6)*, 2059–2095.  https://doi.org/10.1214/aos/1074290337
* Kitamura, Y., Tripathi, G., & Ahn, H. (2004). Empirical likelihood‐based inference in conditional moment restriction models. *Econometrica, 72(6)*, 1667–1714. https://doi.org/10.1111/j.1468-0262.2004.00550.x

Very preliminary and incomplete! Do not circulate!

## TODO

### High priority

* If `x` contains duplicates, `DCV(x, bw = bw.grid, weights = w)` complains that no duplicates were found (see the example)
* Test the previous output of `weightedEL()` and `cemplik()` with the new version, add tests
* De-duplicate at kernel weights already (via `.prepareKernel`), return the attribute!
* For `.prepareKernel()` AND mixed kernel: check if the max. column-wise gap between observations is >= than the bandwidth, otherwise write an informative message
* `bw.CV()` bug: `bw.CV(x, y = y, kernel = "quartic", order = 4)` and `bw.CV(x, y = y, kernel = "quartic")`
* LOO estimation: instead of dropping unique (X, Y) observations, leave each conditioning points (only X)
* Add weight support to `kernelDiscreteDensitySmooth`
* Add `RcppParallel::setThreadOptions(numThreads = "auto")` as the 1st line of parallel-capable functions, use `setDTthreads` also
* Write test cases for C++ functions with and without speed-ups
* Eliminate matrices in smoothing completely! Try only parallel loops!
* Create a default value for memsave, when to invoke it (based on nx*ng)
* Remove parallelisation over workers via setThreadOptions when there is outer parallelisation in `.kernelMixed()`
* Make DCV either sparse or memsave, not both; reflect the changes in `bw.CV()`
* Fix the DCV code with convolutions (especially the quartic one)
* Move the de-duplication of the xout grid inside `kernelSmooth`
* Check all instances of `kernelSmooth()`, `kernelDensity()`, `kernelWeights()`, and everything that used the obsolete arguments
* The check `CV = "DCV"` and `is.null(y)` seems redundat
* Check where `kernelWeights` and `standardise` were used.
* Fix the optimiser control argument in `bw.CV()`.
* Extend the CV to non-Gaussian cases.
* `smoothEmplikDiscrete()`: if the split variable does not take contiguous values from 1 to the number of categories, estimation fails.

### Medium priority

* CV: implement leave-K-out CV for speed
* In `kernelMixedSmooth()`: if LOO, do not de-duplicate `xout`, copy it from `arg$x` (currently mitigated via `deduplicate.xout = FALSE`)
* All LOO to the C++ density function
* Add custom kernels to Silverman's rule of thumb (with roughness != 1)
* Check: if the kernel is finite-support and bandwidth is smaller than the largest gap between two observations, then
* Like in the SEL application: de-duplicate the input matrix, replace with weights; allow the user to disable it
* Merging cells: allow arbitrary variables (including continuous ones) for proximity.
* `kernelSmooth()` and `kernelDensity()` should have an argument for increasing small bandwidths in case of zero weights to match the largest gap divided by 2 (times 1.1 to have at least some coverage)

### Low priority

* Create convolution for kernel orders 4 and 6;
* Reproduce the CKT (2019) results with the `shift` argument (i.e. test the shift)
* Create a summary class for SEL; print numerical gradients of lambdas; print the number of converged inner optimisation problems

## Note for package development

* Check release with `lintr::lint_package()`
* Add tests reproducing simple hard-coded examples
