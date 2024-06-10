<!-- badges: start -->
[![R-CMD-check](https://github.com/Fifis/smoothemplik/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Fifis/smoothemplik/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/Fifis/smoothemplik/graph/badge.svg?token=LYXLUYWY5X)](https://codecov.io/gh/Fifis/smoothemplik)
<!-- badges: end -->

# smoothemplik

R package for estimation via Smoothed Empirical Likelihood (SEL).

<img src="https://kostyrka.lu/user/pages/05.programming/10.smoothemplik-package/smoothed-empirical-likelihood-r-package.png" alt="Smoothed Empirical Likelihood" width="400"/>

Smoothed Empirical Likelihood (SEL), also known as Conditional Empirical Likelihood (CEL) or Local Empirical Likelihood (LEL), is a powerful statistical method for estimation and hypothesis testing.
It shares the competitive advantages of Generalised Empirical Likelihood (GEL) and can estimate models defined by conditional moment restrictions (CMR).

For statistical inference, SEL does not require any explicit variance estimation because the smoothed empirical likelihood ratio (SELR) statistic is internally studentised.
Consequently, SEL-based confidence intervals and regions have accurate coverage probabilities.
To test a hypothesis using SEL, the user should compare the difference between the unrestricted and restricted SEL values with the critical values of the chi-squared distribution.
The SELR test achieves maximum average local power.

Some of the models that can be estimated via SEL are:
- Linear models with exogenous and endogenous variables: *Y = α + X'β + Z'γ + U*, *E(U | X, W) = 0*;
- Non-linear models: *E(Y | X) = f(X, θ)*, where *f* is known;
- Semi-parametric models: *E(Y | X) = X'β + h(X, θ)*, where *h* is unknown;
- Quantile regression models: *E[I(Y - X'β) - τ | X] = 0*.

Additionally, for a broad class of models, SEL can produce an estimator that is semi-parametrically efficient in the sense of [Chamberlain (1987)](https://doi.org/10.1016/0304-4076(87)90015-7).
Efficient estimation is the pipe dream of many empirical researchers because having an unbiased estimator with the smallest possible variance means *accuracy*.
This accuracy leads to four main benefits: (1) having point estimates closer to the true value, (2) having smaller variances and standard errors, (3) discovering more statistically significant effects when they exist, and (4) ensuring that a lack of significance is not due to the weaknesses of the employed estimation method.
For example, the popular ordinary-least-squares (OLS) estimator of a linear model is not efficient under heteroskedasticity – i.e. in virtually all real-world applications.
This happens because ‘noisier’ observations have a larger contribution to the objective function – the sum of squared residuals.
In SEL, the objective function is a sum of local ELR statistics, therefore, all observations contribute similarly to the objective regardless of the conditional variance.

We hope that this algorithmic implementation of SEL will make it more popular among researchers in various fields.

**This package is an evolving project. Comments and suggestions are welcome.**

## How it works?

This package is similar to the popular [gmm](https://CRAN.R-project.org/package=gmm) and [momentfit](https://CRAN.R-project.org/package=momentfit) packages by Pierre Chaussé.
The main input to the SEL optimiser is the function that computes the sample moment condition.
All that is needed is a formalised discrepancy between the model and the real observations.

For a linear model, the user can define the residuals like this:
```{r}
resFun <- function(theta, data) data$Y - data[, c("X1", "X2")] %*% theta
```
Then, `maximiseSEL()` takes it as the input, and the optimiser will find the value of `theta` that yields the ‘best’ fit of the model to the data.
Here, ‘best’ means that these residuals are closest to zero as measured by the sum of local ELR statistics for each observation.

Unlike functions for unconditional-moment-restriction models that typically require input matrices obtained by taking products of GMM instruments with the residuals, the input to SEL is just a function that computes a vector of residuals (or their close analogues).

## Literature

This package provides direct functionality for the following articles:

* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2019). Inference in conditional moment restriction models when there is selection due to stratification. In *The Econometrics of Complex Survey Data* (Vol. 39, pp. 137–171). Emerald Publishing Limited. [doi.org/10.1108/S0731-905320190000039010](https://doi.org/10.1108/S0731-905320190000039010)
* Cosma, A., Kostyrka, A. V., & Tripathi, G. (2024). *Missing endogenous variables in conditional moment restriction models.* (Working paper No. 2024-01). University of Luxembourg, Department of Economics and Management. [hdl.handle.net/10993/60100](hdl.handle.net/10993/60100)

The theory behind the method is provided in the following sources:

* Tripathi, G., & Kitamura, Y. (2003). Testing conditional moment restrictions. *The Annals of Statistics, 31(6)*, 2059–2095.  [10.1214/aos/1074290337](https://doi.org/10.1214/aos/1074290337).
* Kitamura, Y., Tripathi, G., & Ahn, H. (2004). Empirical likelihood‐based inference in conditional moment restriction models. *Econometrica, 72(6)*, 1667–1714. [10.1111/j.1468-0262.2004.00550.x](https://doi.org/10.1111/j.1468-0262.2004.00550.x).
* Owen, A. B. (2013). Self-concordance for empirical likelihood. *Canadian Journal of Statistics, 41*, 387–397.

## Installation

This package currently exists only on GitHub. To install it, run the following two commands:
```{r}
install.packages("devtools")
devtools::install_github("Fifis/smoothemplik")
```

To load this package, include this line in the code:
```{r}
library(smoothemplik)
```

## Licence

This software is released under the free/open-source [EUPL 1.2 licence](https://joinup.ec.europa.eu/collection/eupl/eupl-text-eupl-12).
