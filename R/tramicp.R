#' Invariant causal prediction for transformation models
#'
#' @param formula Formula including response and shift terms
#' @param data Data.frame containing response and explanatory variables
#' @param env Formula specifying the exogenous environment variables
#' @param modFUN Model function from package 'tram', i.e. \code{BoxCox},
#'     \code{Colr}, \code{Polr}, \code{Lm}, \code{Coxph}, \code{Survreg},
#'     \code{Lehmann}. Standard implementations \code{lm}, \code{glm}, and
#'     \code{\link[MASS]{polr}} are also supported.
#' @param verbose Logical, whether output should be verbose (default \code{TRUE})
#' @param type Character, type of invariance (\code{"residual"}, \code{"wald"},
#'     or \code{"partial"})
#' @param test Character, type of test to be used (\code{"HSIC"}, \code{"t.test"},
#'     \code{"var.test"}, \code{"wald"}) or custom function for testing
#'     invariance.
#' @param controls Controls for the used tests and the overall procedure,
#'     see \code{dicp_controls}
#' @param alpha Level of invariance test, default
#' @param ... Further arguments passed to \code{modFUN}
#' @param baseline_fixed Fixed baseline transformation, see
#'     \code{\link[tramicp]{dicp_controls}}.
#' @param greedy Logical, whether to perform a greedy version of ICP (default is
#'     \code{FALSE})
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
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "confint")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
#'     weights = abs(rnorm(nrow(d))))
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual",
#'      test = "HSIC")
#'
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "confint", greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald", greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
#'     weights = abs(rnorm(nrow(d))), greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual", greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual",
#'      test = "HSIC", greedy = TRUE)
#'
dicp <- function(
  formula, data, env, modFUN, verbose = TRUE,
  type = c("residual", "wald", "mcheck", "confint"),
  test = "independence", controls = NULL, alpha = 0.05,
  baseline_fixed = TRUE, greedy = FALSE, ...
) {

  call <- match.call()

  ### Preliminary checks
  if (is.character(test))
    test <- match.arg(
      test, c("independence", "HSIC", "t.test", "var.test", "combined", "wald")
    )
  if (is.null(controls))
    controls <- dicp_controls(match.arg(type), test, alpha = alpha,
                              baseline_fixed = baseline_fixed)
  .check_args(formula, data, env, modFUN, type, test)

  ### Process formulae
  tms <- .get_terms(formula)
  resp <- tms$response
  etms <- .get_terms(env)
  me <- setdiff(tms$all, etms$all) # no env in main effects
  ps <- lapply(0:length(me), combn, x = length(me))

  # Options
  if (verbose && interactive())
    pb <- txtProgressBar(min = 0, max = length(ps), style = 3)

  if (!greedy) {
    # Run
    tests <- list()
    for (set in seq_along(ps)) {

      if (verbose && interactive())
        setTxtProgressBar(pb, set)

      ret <- apply(ps[[set]], 2, controls$type_fun, me = me, resp = resp,
                   set = set, env = etms, modFUN = modFUN, data = data,
                   controls = controls, ... = ...)

      tests <- c(tests, ret)

    }
  } else {
    if (length(me) <= 1)
      stop("Run greedy ICP only with more than one predictor.")
    # Run
    tests <- list()
    set <- seq_along(me)
    while (length(set) > 1) {

      if (verbose && interactive())
        setTxtProgressBar(pb, length(me) - length(set) + 1)

      sets <- combn(set, length(set) - 1)

      ret <- apply(sets, 2, controls$type_fun, me = me, resp = resp,
                   set = length(set), env = etms, modFUN = modFUN, data = data,
                   controls = controls, ... = ...)

      pvals <- lapply(ret, \(x) .get_pvalue(x$test))
      tests <- c(tests, ret)

      # if (verbose && interactive())
      #   cat("\nRemoving", setdiff(set, sets[, which.max(pvals)]), "\n")

      set <- sets[, which.max(pvals)[1]]
      # inv <- me[set]

      if (any(unlist(pvals) < 0.05)) { # && !all(unlist(pvals) < 0.05)) {
        if (verbose && interactive())
          cat("\nTerminated early.")
        set <- 0
      }
    }
  }

  res <- .extract_results(tests)
  pvals <- structure(res[["pval"]], names = res[["set"]])
  # if (!greedy) {
  inv <- try(.inv_set(res, alpha = controls$alpha))
  if (inherits(inv, "try-error"))
    inv <- "Cannot be computed."
  # }
  ipv <- .indiv_pvals(me, pvals)

  structure(list(candidate_causal_predictors = if (identical(inv, character(0)))
    "Empty" else inv, set_pvals = pvals, predictor_pvals = ipv,
    tests = tests), class = "dICP", type = match.arg(type),
    test = controls$test_name, env = env, greedy = greedy, call = call)

}

# Pvalues for individual predictors being a causal parent
.indiv_pvals <- function(terms, pvals) {
  res <- lapply(terms, \(term) suppressWarnings(
    max(pvals[!grepl(term, names(pvals))], na.rm = TRUE)))
  structure(unlist(res), names = terms)
}

# GOF test
.modelcheck <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  if (length(env$all) != 1)
    stop("`type = \"mcheck\"` not implemented for multivariable environments")

  env <- env$all

  # Empty set skipped
  if (set == 1)
    return(structure(list(set = "1", test = list("p.value" = 0),
                          coef = NA, logLik = NA, tram = NA),
                     class = "dICPtest"))

  # Set up formula
  tset <- me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- as.formula(
    paste0(resp, ifelse(controls$baseline_fixed, "", paste("|", env)), "~", meff)
  )

  # Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  # Test
  if ("glm" %in% class(m) && m$family$family == "binomial")
    r <- matrix(residuals.binglm(m), ncol = 1)
  else
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
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  if (length(env$all) != 1)
    stop("`type = \"confint\"` not implemented for multivariable environments")

  env <- env$all
  stopifnot(is.factor(data[[env]]))

  # Prepare formula
  tset <- if (set == 1) "1" else me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- as.formula(
    paste0(resp, ifelse(controls$baseline_fixed, "", paste("|", env)), "~", meff)
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
    tx, me, resp, set, env, modFUN, data,
    controls, ...
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
    paste0(resp, ifelse(
      controls$baseline_fixed, "", paste("|", paste0(env$all, collapse = "+"))), "~",
      meff)
  )

  # Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  # Test
  if ("glm" %in% class(m) && m$family$family == "binomial")
    r <- matrix(residuals.binglm(m), ncol = 1)
  else
    r <- matrix(residuals(m), ncol = 1)
  e <- .rm_int(model.matrix(as.formula(env$fml), data = data))
  tst <- controls$test_fun(r, e, controls)

  structure(list(set = tset, test = tst, coef = coef(m),
                 logLik = logLik(m), tram = m$tram), class = "dICPtest")

}

# Full invariance
.dicp_full <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  if (length(env$all) != 1)
    stop("`type = \"wald\"` not implemented for multivariable environments")

  env <- env$all

  # Empty set treated separately
  if (set == 1) {
    tset <- "1"
    meff <- ifelse(controls$baseline_fixed, env, "1")
    mint <- ""
  } else {
    tset <- me[tx]
    meff <- if (controls$baseline_fixed) paste0(c(me[tx], env), collapse = "+") else
      paste0(me[tx], collapse = "+")
    mint <- paste0(c(paste0(me[tx], ":", env)), collapse = "+")
  }

  # Prepare formula
  mfm <- as.formula(
    paste0(resp, ifelse(controls$baseline_fixed, "", paste0("|", env)),
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
.dicp_partial <- function(tx, me, resp, set, env, modFUN,
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

.t_test <- function(r, e, controls) {
  t.test(r ~ e)
}

.var_test <- function(r, e, controls) {
  var.test(r ~ e)
}

.combined_test <- function(r, e, controls) {
  shift <- t.test(r ~ e)
  scale <- var.test(r ~ e)
  list("p.value" = 2 * min(shift$p.value, scale$p.value),
       shift = shift, scale = scale)
}

#' @importFrom coin independence_test pvalue
.indep_test <- function(r, e, controls) {
  tst <- independence_test(
    r ~ e, teststat = controls$teststat, distribution = controls$distribution,
    xtrafo = controls$xtrafo, ytrafo = controls$ytrafo
  )
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
      "t.test" = .t_test,
      "var.test" = .var_test,
      "independence" = .indep_test,
      "combined" = .combined_test,
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
