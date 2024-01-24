<!-- badges: start -->
[![R-CMD-check](https://github.com/LucasKook/tramicp/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/LucasKook/tramicp/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Model-based causal feature selection for general response types

Package `tramicp` [1] implements invariant causal prediction (ICP) [2] for
transformation models [3], including binary logistic regression, Weibull
regression, the Cox model, linear regression and many others. Methods for other
generalized linear models are also provided. The aim of ICP is to discover the
direct causes of a response given data from heterogeneous experimental settings
and a potentially large pool of candidate predictors. Methodological details are
described in [the paper](https://doi.org/10.48550/arXiv.2309.12833).

# Installation

The development version of `tramicp` package can be installed via:
```r
# install.packages("remotes")
remotes::install_github("LucasKook/tramicp")
```

A stable version is available on [CRAN](https://CRAN.R-project.org/package=tramicp):
```r
install.packages("tramicp")
```

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

The code chunk below shows how to use TRAMICP on the data above.
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

# General usage

In full generality, TRAMICP is implemented in the `dicp()` function, which takes
the argument `modFUN` (model function). For instance, the `glmICP` call from
above is equivalent to `dicp(..., modFUN = glm, family = "binomial")`.

# Implemented model classes

Instead of using `dicp()`, `tramicp` directly implements several model classes
with an alias, as shown in the table below.

| **Function alias**  | **Corresponding `modFUN`** |
|---------------------|----------------------------|
| `BoxCoxICP()`       | `tram::BoxCox()`           | 
| `ColrICP()`         | `tram::Colr()`             |
| `cotramICP()`       | `cotram::cotram()`         |
| `CoxphICP()`        | `tram::Coxph()`            |
| `coxphICP()`        | `survival::coxph()`        |
| `glmICP()`          | `stats::glm()`             |
| `LehmannICP()`      | `tram::Lehmann()`          |
| `LmICP()`           | `tram::Lm()`               |
| `lmICP()`           | `stats::lm()`              |
| `PolrICP()`         | `tram::Polr()`             |
| `polrICP()`         | `MASS::polr()`             |
| `SurvregICP()`      | `tram::Survreg()`          |
| `survregICP()`      | `survival::survreg()`      |

Other implementations, such as additive TRAMs in `tramME`, can still be used via
the `dicp()` function, for instance, after loading `tramME`, `dicp(..., modFUN =
"BoxCoxME")` can be used.

Nonparametric ICP via the GCM test [4] and random forests for the two
regressions is implemented in the alias `rangerICP()`. Survival forests
are supported for right-censored observations and implemented in 
`survforestICP()`.

# Replication materials

This repository contains the code for reproducing the results in [1] in
the `inst` directory. Please follow the instructions in 
[the README](inst/README.md) to run the code.

# References

[1] Kook L., Saengkyongam S., Lundborg A., Hothorn T., Peter J. (2023) 
Model-based causal feature selection for general response types. arXiv preprint.
[doi:10.48550/arXiv.2309.12833](https://doi.org/10.48550/arXiv.2309.12833)

[2] Peters, J., Bühlmann, P., & Meinshausen, N. (2016). Causal inference by 
using invariant prediction: identification and confidence intervals. Journal of 
the Royal Statistical Society Series B: Statistical Methodology, 78(5), 947-1012.
[doi:10.1111/rssb.12167](http://dx.doi.org/10.1111/rssb.12167)

[3] Hothorn, T., Möst, L., & Bühlmann, P. (2018). Most likely transformations.
Scandinavian Journal of Statistics, 45(1), 110-134.
[doi:10.1111/sjos.12291](http://dx.doi.org/10.1111/sjos.12291)

[4] Shah, R. D., & Peters, J. (2020). The hardness of conditional independence
testing and the generalised covariance measure.
[doi:10.1214/19-aos1857](http://dx.doi.org/10.1214/19-aos1857)

