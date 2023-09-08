# smoothemplik

R package for estimation via Smooth Empirical Likelihood

Functionality for the following articles:

* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2019). Inference in conditional moment restriction models when there is selection due to stratification. In *The Econometrics of Complex Survey Data* (Vol. 39, pp. 137–171). Emerald Publishing Limited. https://doi.org/10.1108/S0731-905320190000039010
* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2023). Missing endogenous variables in conditional moment restriction models. *In progress.*

Based on:

* Tripathi, G., & Kitamura, Y. (2003). Testing conditional moment restrictions. *The Annals of Statistics, 31(6)*, 2059–2095.  https://doi.org/10.1214/aos/1074290337
* Kitamura, Y., Tripathi, G., & Ahn, H. (2004). Empirical likelihood‐based inference in conditional moment restriction models. *Econometrica, 72(6)*, 1667–1714. https://doi.org/10.1111/j.1468-0262.2004.00550.x

Very preliminary and incomplete! Do not circulate!

## TODO

### High priority

* The parallelisation in .kernelMixed does not work
* Add custom kernels to Silverman's rule of thumb (with roughness != 1)
* Fix the DCV code with convolutions (especially the quartic one)
* Fix the bug: no speed-up in the `kernelMixedSmooth` example
* Check all instances of kernelSmooth, kernelDensity, kernelWeights, and everything that used the obsolete argument.
* Remove the `CV = "DCV"` or `"LSCV"` from the examples
* Check there `kernelWeights` and `standardise` were used.
* Extend the CV to non-Gaussian cases.
* `smoothEmplikDiscrete`: if the split variable does not take contiguous values from 1 to the number of categories, estimation fails.

### Medium priority

* Check: if the kernel is finite-support and bandwidth is smaller than the largest gap between two observations, then
* Merging cells: allow arbitrary variables (including continuous ones) for proximity.


### Low priority

* Remove AER from suggestions (too many dependencies)
* Create convolution for kernel orders 4 and 6;
* Reproduce the CKT (2019) results with the `shift` argument (i.e. test the shift)
* Create a summary class for SEL; print numerical gradients of lambdas; print the number of converged inner optimisation problems

