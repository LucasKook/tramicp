# Model aliases

#' ICP for Box-Cox type regression models
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "boxcox")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E)
#'
BoxCoxICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "gcm.test", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = BoxCox,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for parametric survival models
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "weibull")
#' SurvregICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#'
SurvregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "gcm.test", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for parametric survival models
#' @rdname survival.survregICP
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @importFrom survival survreg
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "weibull")
#' if (require("survival")) {
#'   d$surv <- Surv(d$Y)
#'   survregICP(surv ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
survregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                        test = "gcm.test", controls = NULL, alpha = 0.05,
                        baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                        mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Continuous outcome logistic regression
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "colr")
#' ColrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, type = "wald", test = "wald")
#'
ColrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "gcm.test", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Colr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Cox proportional hazards regression
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "coxph")
#' CoxphICP(Y ~ X2, data = d, env = ~ E)
#'
CoxphICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "gcm.test", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Coxph,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Lehmann regression models
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "coxph")
#' LehmannICP(Y ~ X2, data = d, env = ~ E)
#'
LehmannICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "gcm.test", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Lehmann,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for normal linear regression
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "lm")
#' LmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#'
LmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                  test = "gcm.test", controls = NULL, alpha = 0.05,
                  baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                  mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for normal linear regression (using stats::lm)
#' @rdname statslmICP
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "lm")
#' lmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#'
lmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "gcm.test", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for cumulative ordinal regression
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "polr")
#' PolrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#'
PolrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "gcm.test", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Polr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for cumulative ordinal regression using \code{MASS::polr()}
#' @rdname masspolrICP
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @importFrom MASS polr
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "polr")
#' polrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#'
polrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "gcm.test", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     mandatory = NULL, ...) {
  call <- match.call()
  if (type == "confint")
    stop("`type = \"confint\"` not implemented for ICP with `MASS::polr()`.")
  modF <- polr
  if (type != "residual")
    modF <- function(formula, data, ...) polr(formula, data, Hess = TRUE, ...)
  if (is.character(test))
    test <- match.arg(test, .implemented_tests())
  if (is.null(controls))
    controls <- dicp_controls(match.arg(
      type, c("residual", "wald", "mcheck", "confint")), test, alpha = alpha)
  controls$vcov <- function(object) {
    cf <- coef(object)
    vcov <- vcov(object)
    vcov[names(cf), names(cf)]
  }
  ret <- dicp(formula = formula, data = data, env = env, modFUN = modF,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' @method residuals polr
#' @importFrom sandwich estfun
#'
residuals.polr <- function(object, ...) {
  K <- length(object$zeta)
  sc <- estfun(object)
  sc <- sc[, rev(seq_len(ncol(sc)))]
  -unname(rowSums(sc[, seq_len(K)]))
}

#' ICP for generalized linear models
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "binary")
#' glmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, family = "binomial")
#'
glmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "gcm.test", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   mandatory = NULL, ...) {
  call <- match.call()
  if (is.character(test))
    test <- match.arg(test, .implemented_tests())
  if (!is.null(fam <- list(...)$family) && (identical(fam, binomial) || fam == "binomial")) {
    resid <- "residuals.binglm"
  }
  if (is.null(controls))
    controls <- dicp_controls(type, test, alpha = alpha,
                              baseline_fixed = baseline_fixed,
                              residuals = resid)
  ret <- dicp(formula = formula, data = data, env = env, modFUN = glm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for count transformation models
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "cotram")
#' cotramICP(Y ~ X2, data = d, env = ~ E)
#'
cotramICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "gcm.test", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = cotram::cotram,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}
