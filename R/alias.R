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
#' BoxCoxICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, type = "wald", test = "wald")
#' BoxCoxICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, type = "wald", test = "wald",
#'     greedy = TRUE)
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
BoxCoxICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "independence", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = BoxCox,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
#' library("survival")
#' d$surv <- Surv(d$Y)
#' dicp(surv ~ X2, data = d, env = ~ E, type = "residual", modFUN = survreg,
#'     test = "indep")
#'
SurvregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "independence", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' @importFrom survival survreg
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(mod = "weibull")
#' library("survival")
#' d$surv <- Surv(d$Y)
#' ssurvregICP(surv ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#'
ssurvregICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                        test = "independence", controls = NULL, alpha = 0.05,
                        baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                        ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = survreg,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' d <- dgp_dicp(mod = "colr")
#' ColrICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' ColrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' ColrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' ColrICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
ColrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "independence", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Colr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' d <- dgp_dicp(mod = "coxph")
#' CoxphICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' CoxphICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' CoxphICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' CoxphICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
CoxphICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "independence", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Coxph,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' d <- dgp_dicp(mod = "coxph")
#' LehmannICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' LehmannICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' LehmannICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' LehmannICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
LehmannICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                       test = "independence", controls = NULL, alpha = 0.05,
                       baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                       ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Lehmann,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' d <- dgp_dicp(mod = "lm")
#' LmICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' LmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' LmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' LmICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#' dicp(Y ~ X2, data = d, env = ~ E, type = "confint", modFUN = lm)
#'
LmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                  test = "independence", controls = NULL, alpha = 0.05,
                  baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                  ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
  ret
}

#' ICP for normal linear regression (using stats::lm)
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @importFrom stats lm
#'
#' @examples
#' d <- dgp_dicp(mod = "lm")
#' slmICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' slmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' slmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' slmICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
slmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "independence", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = lm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' d <- dgp_dicp(mod = "polr")
#' PolrICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' PolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' PolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' PolrICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
PolrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                    test = "independence", controls = NULL, alpha = 0.05,
                    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                    ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = Polr,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
  ret
}

#' ICP for cumulative ordinal regression using \code{MASS::polr()}
#'
#' @inheritParams dicp
#' @inheritDotParams dicp
#'
#' @return Object of type \code{"dICP"}. See \code{\link[tramicp]{dicp}}
#' @export
#'
#' @importFrom MASS polr
#'
#' @examples
#' d <- dgp_dicp(mod = "polr")
#' mpolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' # Almost identical to `PolrICP()`
#' pvalues(mpolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep"), which = "set")
#' pvalues(PolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep"), which = "set")
#' pvalues(mpolrICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald"), which = "set")
#' # `"confint"` Not implemented yet (will throw exception):
#' # mpolrICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
mpolrICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                     test = "independence", controls = NULL, alpha = 0.05,
                     baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                     ...) {
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
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' set.seed(1334)
#' d <- dgp_dicp(mod = "binary")
#' glmICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald", family = "binomial")
#' glmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC", family = "binomial")
#' glmICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep", family = "binomial")
#' glmICP(Y ~ X2, data = d, env = ~ E, type = "confint", family = "binomial")
#' dicp(Y ~ X2, data = d, env = ~ E, modFUN = tramicp:::.mod_from_name("binary"))
#'
glmICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                   test = "independence", controls = NULL, alpha = 0.05,
                   baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                   ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = glm,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
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
#' set.seed(13312)
#' d <- dgp_dicp(mod = "cotram")
#' cotramICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' cotramICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' cotramICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' cotramICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#' glmICP(Y ~ X2, data = d, env = ~ E, family = "poisson", test = "indep")
#'
cotramICP <- function(formula, data, env, verbose = TRUE, type = "residual",
                      test = "independence", controls = NULL, alpha = 0.05,
                      baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
                      ...) {
  call <- match.call()
  ret <- dicp(formula = formula, data = data, env = env, modFUN = cotram::cotram,
              verbose = verbose, type = type, test = test, controls = controls,
              alpha = alpha, baseline_fixed = baseline_fixed, greedy = greedy,
              max_size = max_size, ... = ...)
  attr(ret, "call") <- call
  ret
}
