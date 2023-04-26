
#' TRAMICP Controls
#'
#' @inheritParams dicp
#' @param baseline_fixed logical; whether or not the baseline transformation
#'     is allowed to vary with the environments. Only takes effect when
#'     \code{type} is one of \code{"wald"}, \code{"confint"}, or \code{"mcheck"}.
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
    type = "residual", test = "independence", baseline_fixed = TRUE, alpha = 0.05
) {

  # Type of ICP
  type_fun <- .type_fun(type)

  # Type of test
  ctest <- if (!is.function(test)) test else "custom"
  test_info <- .test_fun(type, test, ctest)

  list(
    type = type, method = "gamma", kernel = c("gaussian", "discrete"), B = 200,
    alpha = 0.05, vcov = vcov, type_fun = .type_fun(type),
    ctest = ctest, test_name = test_info[[3]], test_fun = test_info[[2]],
    baseline_fixed = baseline_fixed, teststat = "maximum",
    distribution = "asymptotic", xtrafo = trafo, ytrafo = trafo
  )
}

.type_fun <- function(type) {
  switch(
    type,
    "residual" = .dicp_hsic,
    "wald" = .dicp_full,
    "partial" = .dicp_partial,
    "mcheck" = .modelcheck,
    "confint" = .crmethod
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
      "custom" = identity
    )
  }

  if (type == "confint")
    ctest <- "confint"

  list(test = test, test_fun = test_fun, test_name = ctest)
}
