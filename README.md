# smoothemplik

R package for estimation via Smooth Empirical Likelihood

Very preliminary and incomplete! Do not circulate!

## TODO

### High priority

* `smoothEmplikDiscrete`: if the split variable does not take contiguous values from 1 to the number of categories, estimation fails.

### Medium priority

* Merging cells: allow arbitrary variables (including continuous ones) for proximity.
* Create a cluster in gradParallel for Windows machines inside the function.

### Low priority

* Remove AER from suggestions (too many dependencies)
* Reproduce the CKT (2019) results with the `shift` argument (i.e. test the shift)
* Create a summary class for SEL; print numerical gradients of lambdas; print the number of converged inner optimisation problems

