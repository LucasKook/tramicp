
#' TRAMICP Controls
#'
#' @inheritParams dicp
#' @param baseline_fixed Logical; whether or not the baseline transformation
#'     is allowed to vary with the environments. Only takes effect when
#'     \code{type} is \code{"wald"}.
#' @param method Only applies if \code{test = "HSIC"}. See
#'     \code{\link[dHSIC]{dhsic.test}}.
#' @param B For \code{test = "HSIC"}, see \code{\link[dHSIC]{dhsic.test}}.
#' @param kernel Only applies if \code{test = "HSIC"}. See \code{\link[dHSIC]{dhsic.test}}.
#' @param vcov (Name of) function for computing the variance-covariance matrix of a model.
#' @param teststat Only applies if \code{test = "independence"}.
#'     See \code{\link[coin]{independence_test}}.
#' @param distribution Only applies if \code{test = "independence"}.
#'     See \code{\link[coin]{independence_test}}.
#' @param xtrafo Only applies if \code{test = "independence"}.
#'     See \code{\link[coin]{independence_test}}.
#' @param ytrafo Only applies if \code{test = "independence"}.
#'     See \code{\link[coin]{independence_test}}.
#' @param residuals Character or function; (Name of) function for computing
#'     model residuals. The default is \code{stats::residuals} with methods
#'     dispatch.
#' @param crossfit Logical; toggle for cross fitting when \code{type = "residual"}.
#' @param stop_if_empty_set_invariant Logical; \code{dicp} halts if the empty
#'     set is not rejected (the resulting intersection will always be empty).
#'     Default is \code{FALSE} and can be over-written by setting
#'     \code{options(stop_if_empty_set_invariant = TRUE)}.
#' @param wald_test_interactions Logical; whether to test for interactions between
#'     residuals and environments when using \code{type = "wald"}
#'     (\code{wald_test_interactions = TRUE}, the default) or main effects only
#'     (\code{wald_test_interactions = FALSE}).
#'
#' @return List of dicp controls containing the evaluated arguments from above.
#'
#' @import tram
#'
#' @export
#'
dicp_controls <- function(
    type = "residual", test = "gcm.test", baseline_fixed = TRUE,
    alpha = 0.05, method = "gamma", kernel = c("gaussian", "discrete"),
    B = 499, vcov = "vcov", teststat = "maximum", distribution = "asymptotic",
    xtrafo = coin::trafo, ytrafo = coin::trafo, residuals = "residuals",
    crossfit = getOption("crossfit", default = FALSE),
    stop_if_empty_set_invariant = getOption("stop_if_empty_set_invariant", default = FALSE),
    wald_test_interactions = getOption("wald_test_interactions", default = TRUE)
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
    residuals = residuals, crossfit = crossfit,
    stop_if_empty_set_invariant = stop_if_empty_set_invariant,
    wald_test_interactions = wald_test_interactions
  )
}

.type_fun <- function(type) {
  switch(
    type,
    "residual" = .residual_invariance,
    "wald" = .wald_invariance,
    "partial" = .partial_invariance
  )
}

.test_fun <- function(type, test, ctest) {
  if (is.function(test))
    return(list(test = "custom", test_fun = test, test_name = ctest))

  if (type %in% c("wald", "partial")) {
    ctest <- "wald"
    test_fun <- multcomp::Chisqtest
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

  list(test = test, test_fun = test_fun, test_name = ctest)
}
