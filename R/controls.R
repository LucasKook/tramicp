
#' TRAMICP Controls
#'
#' @inheritParams dicp
#' @param baseline_fixed logical; whether or not the baseline transformation
#'     is allowed to vary with the environments. Only takes effect when
#'     \code{type} is one of \code{"wald"}, \code{"confint"}, or \code{"mcheck"}.
#' @param method See \code{\link[dHSIC]{dhsic.test}}.
#' @param B See \code{\link[dHSIC]{dhsic.test}}.
#' @param kernel See \code{\link[dHSIC]{dhsic.test}}.
#' @param vcov (Name of) function for computing the variance-covariance matrix of a model.
#' @param teststat See \code{\link[coin]{independence_test}}.
#' @param distribution See \code{\link[coin]{independence_test}}.
#' @param xtrafo See \code{\link[coin]{independence_test}}.
#' @param ytrafo See \code{\link[coin]{independence_test}}.
#' @param residuals (Name of) function for computing model residuals.
#' @param stop_if_empty_set_invariant Logical; \code{dicp} halts if the empty
#'     set is not rejected (the resulting intersection will always be empty).
#'     Default is \code{FALSE} and can be over-written by setting
#'     \code{options(stop_if_empty_set_invariant = TRUE)}.
#'
#' @return List of dicp controls
#'
#' @export
#'
#' @importFrom methods as
#' @importFrom stats as.formula binom.test binomial coef df filter glm logLik
#'     model.frame model.response optimize plogis predict qchisq qlogis rbinom
#'     residuals rlogis rmultinom rnorm sd simulate t.test terms update var.test
#'     vcov model.matrix dexp dlogis dnorm pexp pnorm qexp qnorm
#' @importFrom utils capture.output combn data setTxtProgressBar stack
#'     txtProgressBar write.csv
#' @importFrom coin trafo
#'
dicp_controls <- function(
    type = "residual", test = "gcm.test", baseline_fixed = TRUE,
    alpha = 0.05, method = "gamma", kernel = c("gaussian", "discrete"),
    B = 499, vcov = "vcov", teststat = "maximum", distribution = "asymptotic",
    xtrafo = trafo, ytrafo = trafo, residuals = "residuals",
    stop_if_empty_set_invariant = getOption("stop_if_empty_set_invariant", default = FALSE)
) {

  # Type of ICP
  type_fun <- .type_fun(type)
  vcov <- match.fun(vcov)
  residuals <- match.fun(residuals)

  # Type of test
  ctest <- if (!is.function(test)) test else "custom"
  test_info <- .test_fun(type, test, ctest)

  list(
    type = type, method = method, kernel = kernel, B = B,
    alpha = alpha, vcov = vcov, type_fun = .type_fun(type),
    ctest = ctest, test_name = test_info[[3]], test_fun = test_info[[2]],
    baseline_fixed = baseline_fixed, teststat = teststat,
    distribution = distribution, xtrafo = xtrafo, ytrafo = ytrafo,
    residuals = residuals,
    stop_if_empty_set_invariant = stop_if_empty_set_invariant
  )
}

.type_fun <- function(type) {
  switch(
    type,
    "residual" = .residual_invariance,
    "wald" = .wald_invariance,
    "partial" = .partial_invariance,
    "mcheck" = .gof_invariance,
    "confint" = .confint_invariance
  )
}

#' @importFrom multcomp Chisqtest
.test_fun <- function(type, test, ctest) {
  if (is.function(test))
    return(list(test = "custom", test_fun = test_fun, test_name = ctest))

  if (type %in% c("wald", "partial")) {
    ctest <- "wald"
    test_fun <- Chisqtest
  } else {
    test_fun <- switch(
      ctest,
      "HSIC" = .dhsic_test,
      "t.test" = .t_test,
      "var.test" = .var_test,
      "independence" = .indep_test,
      "combined" = .combined_test,
      "cor.test" = .cor_test,
      "spearman" = .spearman_test,
      "gcm.test" = .gcm_test,
      "custom" = identity
    )
  }

  if (type == "confint")
    ctest <- "confint"

  list(test = test, test_fun = test_fun, test_name = ctest)
}
