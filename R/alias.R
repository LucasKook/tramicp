# Model aliases

#' @title Aliases for implemented model classes
#' @name implemented_model_classes
#'
#' @description ICP for Box-Cox-type transformed normal regression, parametric
#'     and semiparametric survival models, continuous outcome logistic
#'     regression, linear regression, cumulative ordered regression, generalized
#'     linear models; and nonparametric ICP via ranger. While TRAMICP based on
#'     quantile and survival random forests is also supported, for these methods
#'     it comes without theoretical guarantees as of yet.
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "boxcox", n = 300)
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "wald")
#'
BoxCoxICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "gcm.test", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::BoxCox,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for parametric survival models
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "weibull", n = 300)
#' SurvregICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' ### or
#' library("survival")
#' d$Y <- Surv(d$Y)
#' survregICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' CoxphICP(Y ~ X2, data = d, env = ~ E)
#' coxphICP(Y ~ X2, data = d, env = ~ E)
#' }
#'
SurvregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "gcm.test", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for parametric survival models
#' @rdname tramicp-alias
#' @inheritParams dicp
#'
#' @export
#'
survregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                        test = "gcm.test", controls = NULL, alpha = 0.05,
                        baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                        mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = survival::survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for the classical (semi-parametric) Cox model
#' @rdname tramicp-alias
#' @inheritParams dicp
#'
#' @export
#'
coxphICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "gcm.test", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = survival::coxph,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Continuous outcome logistic regression
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "colr", n = 300)
#' ColrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
ColrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "gcm.test", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Colr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Cox proportional hazards regression
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
CoxphICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "gcm.test", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Coxph,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for Lehmann regression models
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "coxph", n = 300)
#' LehmannICP(Y ~ X2, data = d, env = ~ E)
#' }
#'
LehmannICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "gcm.test", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Lehmann,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for normal linear regression
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "lm", n = 300)
#' LmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' ### or
#' lmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
LmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                  test = "gcm.test", controls = NULL, alpha = 0.05,
                  baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                  mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for normal linear regression (using stats::lm)
#' @rdname tramicp-alias
#' @inheritParams dicp
#'
#' @export
#'
lmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "gcm.test", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = stats::lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for cumulative ordinal regression
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "polr", n = 300)
#' PolrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' ### or
#' PolrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
PolrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "gcm.test", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = tram::Polr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for cumulative ordinal regression using \code{MASS::polr()}
#' @rdname tramicp-alias
#' @inheritParams dicp
#'
#' @export
#'
polrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "gcm.test", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     mandatory = NULL, ...) {
  call <- match.call()
  modF <- MASS::polr
  if (type != "residual")
    modF <- function(formula, data, ...) MASS::polr(formula, data, Hess = TRUE, ...)
  if (is.character(test))
    test <- match.arg(test, .implemented_tests())
  if (is.null(controls))
    controls <- dicp_controls(match.arg(
      type, c("residual", "wald")), test, alpha = alpha)
  controls$vcov <- function(object) {
    cf <- stats::coef(object)
    vcov <- stats::vcov(object)
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
#' @exportS3Method residuals polr
#'
residuals.polr <- function(object, ...) {
  K <- length(object$zeta)
  sc <- sandwich::estfun(object)
  sc <- sc[, rev(seq_len(ncol(sc)))]
  -unname(rowSums(sc[, seq_len(K)]))
}

#' ICP for generalized linear models
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "binary", n = 300)
#' glmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, family = "binomial")
#' }
#'
glmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "gcm.test", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   mandatory = NULL, ...) {
  call <- match.call()
  if (is.character(test))
    test <- match.arg(test, .implemented_tests())
  if (is.null(controls))
    controls <- dicp_controls(type, test, alpha = alpha,
                              baseline_fixed = baseline_fixed,
                              residuals = "residuals.binglm")
  ret <- dicp(formula = formula, data = data, env = env, modFUN = stats::glm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' ICP for count transformation models
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "cotram", n = 300)
#' cotramICP(Y ~ X2, data = d, env = ~ E)
#' }
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

#' nonparametric ICP with ranger GCM
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' d <- dgp_dicp(mod = "binary", n = 300)
#' rangerICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
rangerICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "gcm.test", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      mandatory = NULL, ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = RANGER,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' nonparametric ICP for right-censored observations with survival forest GCM
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(12)
#' d <- dgp_dicp(mod = "coxph", n = 3e2)
#' d$Y <- survival::Surv(d$Y, sample(0:1, 3e2, TRUE, prob = c(0.1, 0.9)))
#' survforestICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
survforestICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                          test = "gcm.test", controls = NULL, alpha = 0.05,
                          baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                          mandatory = NULL, ...) {
  call <- match.call()
  message("Note: `survforestICP()` does not come with theoretical guarantees.")
  ret <- dicp(formula = formula, data = data, env = env, modFUN = survforest,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}

#' nonparametric ICP with quantile forest GCM
#' @rdname tramicp-alias
#'
#' @inheritParams dicp
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(12)
#' d <- dgp_dicp(mod = "boxcox", n = 3e2)
#' qrfICP(Y ~ X1 + X2 + X3, data = d, env = ~ E)
#' }
#'
qrfICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "gcm.test", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   mandatory = NULL, ...) {
  call <- match.call()
  message("Note: `qrfICP()` does not come with theoretical guarantees.")
  ret <- dicp(formula = formula, data = data, env = env, modFUN = qrf,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, mandatory = mandatory, ... = ...)
  ret$call <- call
  ret
}
