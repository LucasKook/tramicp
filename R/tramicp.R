#' Invariant causal prediction for transformation models
#'
#' @param formula A \code{formula} including response and covariate terms.
#' @param data A \code{data.frame} containing response and explanatory variables.
#' @param env A \code{formula} specifying the environment variables (see details).
#' @param modFUN Model function from package 'tram', i.e.,
#'     \code{\link[tram]{BoxCox}}, \code{\link[tram]{Colr}},
#'     \code{\link[tram]{Polr}}, \code{\link[tram]{Lm}},
#'     \code{\link[tram]{Coxph}}, \code{\link[tram]{Survreg}},
#'     \code{\link[tram]{Lehmann}}. Standard implementations
#'     \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#'     \code{\link[survival]{survreg}}, \code{\link[survival]{coxph}},
#'     and \code{\link[MASS]{polr}} are also supported. See the corresponding
#'     alias \code{<model_name>ICP}, e.g., \code{\link{PolrICP}}.
#' @param verbose Logical, whether output should be verbose (default \code{TRUE}).
#' @param type Character, type of invariance (\code{"residual"}, \code{"wald"},
#'     or \code{"partial"}).
#' @param test Character, type of test to be used (\code{"HSIC"}, \code{"t.test"},
#'     \code{"var.test"}, \code{"wald"}) or custom function for testing
#'     invariance.
#' @param controls Controls for the used tests and the overall procedure,
#'     see \code{dicp_controls}.
#' @param alpha Level of invariance test, default \code{0.05}.
#' @param ... Further arguments passed to \code{modFUN}.
#' @param baseline_fixed Fixed baseline transformation, see
#'     \code{\link[tramicp]{dicp_controls}}.
#' @param greedy Logical, whether to perform a greedy version of ICP (default is
#'     \code{FALSE}).
#' @param max_size Numeric; maximum support size.
#'
#' @return Object of class \code{"dICP"}, containing
#'     \itemize{
#'     \item{\code{candidate_causal_predictors}: Character; intersection of all
#'     non-rejected sets,}
#'     \item{\code{set_pvals}: Numeric vector; set-specific p-values of the invariance
#'     test,}
#'     \item{\code{predictor_pvals}: Numeric vector; predictor-specific p-values,}
#'     \item{\code{tests}: List of invariance tests.}
#'     }
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
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "confint",
#'     greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
#'     greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
#'     weights = abs(rnorm(nrow(d))), greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual",
#'     greedy = TRUE)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual",
#'      test = "HSIC", greedy = TRUE)
#'
dicp <- function(
  formula, data, env, modFUN, verbose = TRUE,
  type = c("residual", "wald", "mcheck", "confint"),
  test = "independence", controls = NULL, alpha = 0.05,
  baseline_fixed = TRUE, greedy = FALSE, max_size = NULL, ...
) {

  call <- match.call()

  ### Preliminary checks
  if (is.character(test))
    test <- match.arg(test, .implemented_tests())
  if (is.null(controls))
    controls <- dicp_controls(match.arg(type), test, alpha = alpha,
                              baseline_fixed = baseline_fixed)
  .check_args(formula, data, env, modFUN, type, test)

  ### Process formulae
  tms <- .get_terms(formula)
  resp <- tms$response
  etms <- .get_terms(env)
  me <- setdiff(tms$all, etms$all) # no env in main effects
  if (is.null(max_size))
    max_size <- length(me)
  max_size <- min(max_size, length(me))
  ps <- lapply(0:max_size, combn, x = length(me))

  ### Options
  if (verbose && interactive())
    pb <- txtProgressBar(min = 0, max = length(ps), style = 3)

  ### Run invariant subset search
  out <- .invariant_subset_search(ps = ps, controls = controls, me = me,
                                  resp = resp, etms = etms, modFUN = modFUN,
                                  data = data, greedy = greedy,
                                  verbose = verbose, pb = pb, ... = ...)

  ### Process output
  tests <- out$tests
  res <- .extract_results(tests)

  ### Extract p-values for each set
  pvals <- structure(res[["pval"]], names = res[["set"]])
  if (!greedy) {
    inv <- try(.inv_set(res, alpha = controls$alpha))
    if (inherits(inv, "try-error"))
      inv <- "Cannot be computed."
  } else {
    inv <- me[sort(unique(unlist(out$MI)))]
  }

  ### Compute predictor-level p-values
  ipv <- .indiv_pvals(me, pvals)

  ### Return
  structure(list(candidate_causal_predictors = if (identical(inv, character(0)))
    "Empty" else inv, set_pvals = pvals, predictor_pvals = ipv,
    tests = tests), class = "dICP", type = match.arg(type),
    test = controls$test_name, env = env, greedy = greedy, call = call)

}

# Run invariant subset search
.invariant_subset_search <- function(ps, controls, me, resp, etms, modFUN,
                                     data, greedy, verbose, pb, ...) {

  if (!greedy) {
    ### Run
    tests <- list()
    MI <- NULL
    for (set in seq_along(ps)) {

      if (verbose && interactive())
        setTxtProgressBar(pb, set)

      ret <- apply(ps[[set]], 2, controls$type_fun, me = me, resp = resp,
                   set = set, env = etms, modFUN = modFUN, data = data,
                   controls = controls, ... = ...)

      tests <- c(tests, ret)

    }
  } else {
    ### Run
    tests <- list()
    MI <- list()
    lps <- .unlist_once(lapply(ps, \(x) apply(x, 2, \(y) y, simplify = FALSE)))
    for (set in seq_along(lps)) {

      if (verbose && interactive())
        setTxtProgressBar(pb, set)

      if (length(MI > 0) && any(unlist(MI) %in% lps[[set]])) {
        tests[[set]] <- .empty_output(me[lps[[set]]], 0)
        next
      }

      ret <- controls$type_fun(
        lps[[set]], me = me, resp = resp, set = set, env = etms,
        modFUN = modFUN, data = data, controls = controls, ... = ...
      )

      if (.get_pvalue(ret$test) > controls$alpha) {
        # cat("\nAdding", lps[[set]], "to MI\n")
        MI <- c(MI, lps[[set]])
      }

      UMI <- sort(unique(unlist(MI)))
      if (length(MI) > 0 && length(UMI) == length(me) && all(UMI == seq_along(me)))
        break

      tests[[set]] <- ret

    }
  }

  list(tests = tests, MI = MI)
}

# Type-functions ----------------------------------------------------------

# GOF test
.gof_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  env <- env$all

  ### Empty set skipped
  if (set == 1)
    return(.empty_output("Empty"))

  ### Set up formula
  tset <- me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- reformulate(meff, resp)

  ### Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  ### Test
  r <- matrix(controls$residuals(m), ncol = 1)
  x <- model.matrix(update(mfm, NULL ~ .), data = data)
  x <- if (all(x[, 1] == 1)) x[, -1, drop = FALSE]
  tst <- controls$test_fun(r, x, controls)

  ### Return
  structure(list(set = tset, test = tst, coef = coef(m), tram = m$tram),
            class = "dICPtest")

}

# Overlapping confidence regions
#' @importFrom mlt as.mlt
#' @importFrom stats reformulate
.confint_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  ### Checks
  env <- env$all
  stopifnot("`data[[env]]` needs to be a factor." = is.factor(data[[env]]))

  ### Prepare formula
  tset <- if (set == 1) "1" else me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- reformulate(meff, resp)

  ### Fit model
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

  ### Test
  tst <- list(p.value = optimize(.ci, interval = c(0, 1), maximum = TRUE,
                                 ms = ms, nenv = length(uevs))$maximum)
  ## Return
  if (set == 1) tset <- "Empty"
  structure(list(set = tset, test = tst, coef = lapply(ms, coef),
                 tram = ms[[1]]$tram), class = "dICPtest")

}

# Residual invariance
#' @importFrom ranger ranger
.residual_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  ### Prepare formula
  tset <- if (set == "1") 1 else me[tx]
  meff <- paste0(tset, collapse = "+")
  mfm <- reformulate(meff, resp)

  ### Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  ### Test
  r <- matrix(controls$residuals(m), ncol = 1)
  e <- .rm_int(model.matrix(as.formula(env$fml), data = data))
  if (controls$ctest == "gcm.test" & set != "1")
    e <- .ranger_gcm(e, meff, set, data, controls) # Fit random forest for GCM-type test
  tst <- controls$test_fun(r, e, controls)

  ### Return
  if (set == 1) tset <- "Empty"
  structure(list(set = tset, test = tst, coef = coef(m)), class = "dICPtest")

}

# Wald invariance
.wald_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, ...
) {

  env <- env$all

  ### Empty set treated separately
  if (set == 1) {
    tset <- "1"
    meff <- ifelse(controls$baseline_fixed, env, "1")
    mint <- ""
  } else {
    tset <- me[tx]
    meff <- if (controls$baseline_fixed)
      paste0(c(me[tx], env), collapse = "+")
    else
      paste0(me[tx], collapse = "+")
    mint <- paste0(c(paste0(me[tx], ":", env)), collapse = "+")
  }

  ### Prepare formula
  mfm <- as.formula(
    paste0(resp, ifelse(controls$baseline_fixed, "", paste0("|", env)),
           "~", meff, if (mint != "") "+", mint)
  )

  ### Fit
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  if (inherits(m, "tram"))
    m <- as.mlt(m)
  cfs <- names(coef(m))
  tcfs <- union(grep(paste0(":", env), cfs, value = TRUE),
                grep(env, cfs, value = TRUE))

  ### Test
  tst <- try(summary(glht(m, linfct = paste(tcfs, "== 0"),
                          vcov = controls$vcov),
                     test = controls$test_fun()), silent = FALSE)

  ### Catch failure cases
  if (inherits(tst, "try-error")) {
    return(.empty_output(me[tx], NA))
  }

  ### Return
  if (set == 1) tset <- "Empty"
  tst$set <- tset
  tst

}

# Partial invariance
.partial_invariance <- function(
    tx, me, resp, set, env, modFUN, data, trt, controls, ...
) {

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

  if (set == 1) tset <- "Empty"
  tst$set <- me[tx]
  tst

}
