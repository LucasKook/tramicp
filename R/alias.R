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
#' d <- dgp_dicp(mod = "boxcox")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' BoxCoxICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
BoxCoxICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = BoxCox, ...)
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
#' d <- dgp_dicp(mod = "weibull")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "HSIC")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")
#' SurvregICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
#' library(survival)
#' d$surv <- Surv(d$Y)
#' dicp(surv ~ X2, data = d, env = ~ E, type = "residual", modFUN = survreg,
#'     test = "indep")
#'
SurvregICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Survreg, ...)
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
ColrICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Colr, ...)
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
CoxphICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Coxph, ...)
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
LehmannICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Lehmann, ...)
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
LmICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Lm, ...)
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
PolrICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = Polr, ...)
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
#' mpolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")$pvals
#' PolrICP(Y ~ X2, data = d, env = ~ E, type = "residual", test = "indep")$pvals
#' mpolrICP(Y ~ X2, data = d, env = ~ E, type = "wald", test = "wald")$pvals
#' # Not implemented yet (will throw exception):
#' # mpolrICP(Y ~ X2, data = d, env = ~ E, type = "confint")
#'
mpolrICP <- function(formula, data, env, type = "residual", test = "indep",
                    controls = NULL, ...) {
  dots <- list(...)
  if (type == "confint")
    stop("`type = \"confint\"` not implemented for ICP with `MASS::polr()`.")
  modF <- polr
  if (type != "residual")
    modF <- function(formula, data, ...) polr(formula, data, Hess = TRUE, ...)
  if (is.character(test))
    test <- match.arg(
      test, c("independence", "HSIC", "t.test", "var.test", "combined", "wald")
    )
  if (is.null(controls))
    controls <- dicp_controls(match.arg(
      type, c("residual", "wald", "mcheck", "confint")), test, alpha = alpha)
  controls$vcov <- function(object) {
    cf <- coef(object)
    vcov <- vcov(object)
    vcov[names(cf), names(cf)]
  }
  dicp(formula = formula, data = data, env = env, modFUN = modF, type = type,
       test = test, controls = controls, ...)
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
glmICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = glm, ...)
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
cotramICP <- function(formula, data, env, ...) {
  dicp(formula = formula, data = data, env = env, modFUN = cotram::cotram, ...)
}
