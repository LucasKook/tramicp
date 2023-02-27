#' Invariant causal prediction for transformation models
#'
#' @param formula Formula including response and shift terms
#' @param data Data.frame containing response and explanatory variables
#' @param env Character, name of environmental variable
#' @param modFUN Model function from package 'tram', i.e. \code{BoxCox},
#'     \code{Colr}, \code{Polr}, \code{Lm}, \code{Coxph}, \code{Survreg},
#'     \code{Lehmann}. Standard implementation \code{lm} is also supported.
#' @param verbose Logical, whether output should be verbose (default \code{TRUE})
#' @param type Character, type of invariance (\code{"residual"}, \code{"wald"},
#'     or \code{"partial"})
#' @param test Character, type of test to be used (\code{"HSIC"}, \code{"t.test"},
#'     \code{"var.test"}, \code{"wald"}) or custom function for testing
#'     invariance.
#' @param trt Character, supply only when \code{type = "partial"}, treatment
#'     variable. Ignored otherwise
#' @param baseline_fixed Logical, whether baseline transformation is fixed
#'     (\code{TRUE}) or allowed to vary across environments (\code{FALSE}).
#'     Defaults to \code{TRUE}, i.e. a fixed baseline transformation
#' @param alpha Level of invariance test, default
#' @param ... Further arguments passed to \code{modFUN}
#' @param controls Controls for the used tests, see \code{dicp_controls}
#'
#' @return Object of class \code{"dICP"}, containing the invariant set (if exists),
#'     pvalues from all invariance tests and the tests themselves
#'
#' @export
#' @importFrom dHSIC dhsic.test
#' @importFrom multcomp glht
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(n = 1e3, mod = "polr")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = "E", modFUN = Polr, type = "confint")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = "E", modFUN = Polr, type = "wald")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = "E",
#'    modFUN = tramicp:::.mod_from_name("polr"),
#'    type = "wald", weights = abs(rnorm(nrow(d))))
#' dicp(Y ~ X1 + X2 + X3, data = d, env = "E", modFUN = Polr, type = "residual")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = "E", modFUN = Polr, type = "residual",
#'      test = "indep")
#'
dicp <- function(
  formula, data, env, modFUN, trt = NULL, verbose = TRUE,
  type = c("residual", "wald", "partial", "mcheck", "confint"),
  test = "indep", controls = dicp_controls(match.arg(type), test),
  baseline_fixed = TRUE, alpha = 0.05, ...
) {

  # Process formula
  tms <- .get_terms(formula)
  resp <- tms$response
  me <- setdiff(tms$me, env)
  ps <- lapply(0:length(me), combn, x = length(me))

  # Options
  if (verbose)
    pb <- txtProgressBar(min = 0, max = length(ps), style = 3)

  # Run
  tests <- list()
  for (set in seq_along(ps)) {

    if (verbose)
      setTxtProgressBar(pb, set)

    ret <- apply(ps[[set]], 2, controls$type_fun, me = me, resp = resp,
                 set = set, baseline_fixed = baseline_fixed, env = env,
                 modFUN = modFUN, data = data, trt = trt, controls = controls,
                 ... = ...)

    tests <- c(tests, ret)

  }

  res <- .extract_results(tests)
  pvals <- structure(res[["pval"]], names = res[["set"]])
  inv <- try(.inv_set(res, alpha = alpha))
  if (inherits(inv, "try-error"))
    inv <- "Cannot be computed."
  ipv <- .indiv_pvals(me, pvals)

  structure(list(inv = if (identical(inv, character(0))) "Empty" else inv,
      pvals = pvals, ipvals = ipv, tests = tests
    ), class = "dICP", type = match.arg(type), test = controls$test_name,
    env = env, trt = trt)

}

# Pvalues for individual predictors being a causal parent
.indiv_pvals <- function(terms, pvals) {
  res <- lapply(terms, \(term) max(pvals[!grepl(term, names(pvals))]))
  structure(unlist(res), names = terms)
}

# GOF test
.modelcheck <- function(
    tx, me, resp, set, baseline_fixed, env, modFUN, data, trt, controls, ...
) {

  # Empty set skipped
  if (set == 1)
    return(structure(list(set = "1", test = list("p.value" = 0),
                          coef = NA, logLik = NA, tram = NA),
                     class = "dICPtest"))

  # Set up formula
  tset <- me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- as.formula(
    paste0(resp, ifelse(baseline_fixed, "", paste("|", env)), "~", meff)
  )

  # Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  # Test
  r <- matrix(residuals(m), ncol = 1)
  x <- model.matrix(update(mfm, NULL ~ .), data = data)
  x <- if (all(x[, 1] == 1)) x[, -1]
  tst <- controls$test_fun(r, x, controls)

  # Return
  structure(list(set = tset, test = tst, coef = coef(m),
                 logLik = logLik(m), tram = m$tram), class = "dICPtest")

}

# Overlapping confidence regions
#' @importFrom mlt as.mlt
.crmethod <- function(
    tx, me, resp, set, baseline_fixed, env, modFUN, data, trt, controls, ...
) {

  # Prepare formula
  tset <- if (set == 1) "1" else me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- as.formula(
    paste0(resp, ifelse(baseline_fixed, "", paste("|", env)), "~", meff)
  )

  # Fit model
  uevs <- unique(data[, env])
  dots <- list(...)
  dots$weights <- NULL
  ms <- lapply(
    uevs,
    \(e) {
      edat <- data[data[[env]] == e,]
      m <- do.call(modFUN, c(list(formula = mfm, data = edat), dots))
      if (inherits(m, "tram"))
        m <- as.mlt(m)
      m
    }
  )

  # Test
  tst <- list(p.value = optimize(.ci, interval = c(0, 1), maximum = TRUE,
                                 ms = ms, nenv = length(uevs))$maximum)

  structure(list(set = tset, test = tst, coef = lapply(ms, coef),
                 logLik = lapply(ms, logLik), tram = ms[[1]]$tram),
            class = "dICPtest")

}

# Residual test
.dicp_hsic <- function(
    tx, me, resp, set, baseline_fixed, env, modFUN, data,
    trt, controls, ...
) {

  # Skip empty set
  if (set == 1)
    return(structure(list(set = "1", test = list("p.value" = 0),
                          coef = NA, logLik = NA, tram = NA),
                     class = "dICPtest"))

  # Prepare formula
  tset <- me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- as.formula(
    paste0(resp, ifelse(baseline_fixed, "", paste("|", env)), "~", meff)
  )

  # Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  # Test
  r <- matrix(residuals(m), ncol = 1)
  # e <- matrix(as.numeric(data[[env]]) - 1, ncol = 1)
  e <- model.matrix(as.formula(paste0("~ 0 + ", env)), data = data)
  tst <- controls$test_fun(r, e, controls)

  structure(list(set = tset, test = tst, coef = coef(m),
                 logLik = logLik(m), tram = m$tram), class = "dICPtest")

}

# Full invariance
.dicp_full <- function(
    tx, me, resp, set, baseline_fixed, env, modFUN, data, trt, controls, ...
) {

  # Empty set treated separately
  if (set == 1) {
    tset <- "1"
    meff <- ifelse(baseline_fixed, env, "1")
    mint <- ""
  } else {
    tset <- me[tx]
    meff <- if (baseline_fixed) paste0(c(me[tx], env), collapse = "+") else
      paste0(me[tx], collapse = "+")
    mint <- paste0(c(paste0(me[tx], ":", env)), collapse = "+")
  }

  # Prepare formula
  mfm <- as.formula(
    paste0(resp, ifelse(baseline_fixed, "", paste0("|", env)),
           "~", meff, if (mint != "") "+", mint)
  )

  # Fit
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  if (inherits(m, "tram"))
    m <- as.mlt(m)
  cfs <- names(coef(m))
  tcfs <- union(grep(paste0(":", env), cfs, value = TRUE),
                grep(env, cfs, value = TRUE))

  # Test
  tst <- try(summary(glht(m, linfct = paste(tcfs, "== 0"),
                          vcov = controls$vcov),
                     test = controls$test_fun()), silent = FALSE)

  # Catch failure cases
  if (inherits(tst, "try-error")) {
    empty_res <- list(test = list(p.value = NA), set = me[tx])
    return(empty_res)
  }

  tst$set <- tset

  # Return
  tst

}

# Partial invariance
.dicp_partial <- function(tx, me, resp, set, baseline_fixed, env, modFUN,
                          data, trt, controls, ...) {

  if (set == 1) {
    tset <- "1"
    meff <- paste0(trt, "+", env)
    mint1 <- paste0(trt, ":", env)
    mint2 <- ""
  } else {
    tset <- me[tx]
    meff <- paste0(me, collapse = "+")
    mint1 <- paste0(c(paste0(trt, ":", tset), paste0(env, ":", tset),
                      paste0(trt, ":", env)), collapse = "+")
    mint2 <- paste0(paste0(trt, ":", env, ":", tset), collapse = "+")
  }
  mfm <- as.formula(paste0(resp, "~", meff, "+", mint1, if (mint2 != "") "+",
                           mint2))
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  cfs <- names(coef(m))
  tcfs <- grep(paste0(":", env), grep(paste0(if (set != 1) ":", trt),
                                      cfs, value = TRUE), value = TRUE)
  tst <- try(summary(glht(m, linfct = paste(tcfs, "== 0"),
                          vcov = controls$vcov),
                     test = controls$test_fun()), silent = FALSE)
  if (inherits(tst, "try-error")) {
    empty_res <- list(test = list(p.value = NA), set = me[tx])
    return(empty_res)
  }

  tst$set <- me[tx]
  tst

}

.extract_results <- function(res) {
  ret <- unlist(lapply(res, \(x) {
    ret <- .get_pvalue(x$test)
    names(ret) <- paste(x$set, collapse = "+")
    ret
  }))
  data.frame(set = names(ret), pval = unname(ret))
}

# Tests -------------------------------------------------------------------

.dhsic_test <- function(r, e, controls) {
  dhsic.test(r, e, alpha = controls$alpha, method = controls$method,
             B = controls$B)
}

.t.test <- function(r, e, controls) {
  t.test(r ~ e)
}

.var.test <- function(r, e, controls) {
  var.test(r ~ e)
}

.combined <- function(r, e, controls) {
  shift <- t.test(r ~ e)
  scale <- var.test(r ~ e)
  list("p.value" = 2 * min(shift$p.value, scale$p.value),
       shift = shift, scale = scale)
}

#' @importFrom coin independence_test pvalue
.rand.test <- function(r, e, controls) {
  tst <- independence_test(r ~ e)
  list(p.value = pvalue(tst), test = tst)
}

.get_pvalue <- function(x) {
  ret <- x[["p.value"]]
  if (inherits(x, "gtest"))
    ret <- c(x[["pvalue"]])
  if (is.null(ret))
    warning("supplied test has no entry called `p.value`")
  ret
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
      "t.test" = .t.test,
      "var.test" = .var.test,
      "indep" = .rand.test,
      "combined" = .combined,
      "custom" = identity
    )
  }

  if (type == "confint")
    ctest <- "confint"

  list(test = test, test_fun = test_fun, test_name = ctest)
}

#' Controls
#'
#' @inheritParams dicp
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
#'
dicp_controls <- function(type, test) {

  # Type of ICP
  type_fun <- .type_fun(type)

  # Type of test
  ctest <- if (!is.function(test)) test else "custom"
  test_info <- .test_fun(type, test, ctest)

  list(
    method = "gamma", kernel = c("gaussian", "discrete"), B = 200,
    alpha = 0.05, vcov = vcov, type_fun = .type_fun(type),
    ctest = ctest, test_name = test_info[[3]], test_fun = test_info[[2]]
  )
}
