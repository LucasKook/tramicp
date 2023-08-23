<!-- badges: start -->
[![R-CMD-check](https://github.com/LucasKook/tramicp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/tramicp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Invariant causal prediction for transformation models

Package `tramicp` [1] implements invariant causal prediction (ICP) [2] for 
transformation models [3], including binary logistic regression, Weibull 
regression, the Cox model, linear regression and many others. The aim of ICP
is to discover the direct causes of a response given data from heterogeneous
experimental settings and a potentially large pool of candidate predictors.

# Installation

The `tramicp` package can be installed via
```r
remotes::install_github("LucasKook/tramicp")
```
or locally via `R CMD build`/`R CMD INSTALL`.

# Using package `tramicp`

Consider the following data simulated from a structural causal model:
```r
set.seed(-42)
n <- 5e2
E <- sample(0:1, n, TRUE)
X1 <- -E + rnorm(n)
Y <- as.numeric(0.5 * X1 > rlogis(n))
X2 <- Y + 0.8 * E + rnorm(n)
df <- data.frame(Y = Y, X1 = X1, X2 = X2, E = E)
```
The response `Y` is governed by a logistic regression model with parent
`X1`, `X2` is a child and both `X1` and `X2` are influenced by a binary
environment indicator `E`. `tramicp` can discover the parent of `Y` by
testing whether the score residuals for the models `Y ~ X1`, `Y ~ X2`,
and `Y ~ X1 + X2` are uncorrelated with (a residualized version of) `E`. Under 
the correctly specified model, `Y ~ X1`, this correlation will be zero.
To obtain an estimator for the parent set, ICP takes the intersection over all
sets for which the invariance hypothesis is failed to be rejected.

The code chunk below shows how to use `tramicp` on the data above.
```r
icp <- glmICP(Y ~ X1 + X2, data = df, env = ~ E, family = "binomial",
              test = "gcm.test", verbose = FALSE)
pvalues(icp, "set")
```
```
       Empty           X1           X2        X1+X2 
1.818449e-02 5.096300e-01 4.541276e-09 2.219956e-03 
```
Indeed, the only set which is not rejected is `X1`.

Here, `glmICP()` with `family = "binomial"` is used. The formula
`Y ~ X1 + X2` specifies the response (LHS) and all candidate predictors (RHS).
The environments are also specified as a formula (RHS only). Details on the
test and other options can be found in the manuscript and documentation of the
package.

# Reproducibility

This repository contains the code for reproducing the results in [1] in
the `inst` directory. Please follow the instructions in 
[the README](inst/README.md) to run the code.

# References

[1] Kook L., Saengkyongam S., Lundborg A., Hothorn T., Peter J. (2023) 
Model-based causal feature selection for general response types. Work in
progress.

[2] Peters, J., Bühlmann, P., & Meinshausen, N. (2016). Causal inference by 
using invariant prediction: identification and confidence intervals. Journal of 
the Royal Statistical Society Series B: Statistical Methodology, 78(5), 947-1012.

[3] Hothorn, T., Möst, L., & Bühlmann, P. (2018). Most likely transformations.
Scandinavian Journal of Statistics, 45(1), 110-134.
