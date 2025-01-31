#' Model-based causal feature selection for general response types
#'
#' @description
#' Function `dicp()` implements invariant causal prediction (ICP) for
#' transformation and generalized linear models, including binary logistic
#' regression, Weibull regression, the Cox model, linear regression and many
#' others. The aim of ICP is to discover the direct causes of a response given
#' data from heterogeneous experimental settings and a potentially large pool of
#' candidate predictors.
#'
#' @param formula A \code{formula} including response and covariate terms.
#' @param data A \code{data.frame} containing response and explanatory variables.
#' @param env A \code{formula} specifying the environment variables (see details).
#' @param modFUN Model function from 'tram' (or other packages), e.g.,
#'     \code{\link[tram]{BoxCox}}, \code{\link[tram]{Colr}},
#'     \code{\link[tram]{Polr}}, \code{\link[tram]{Lm}},
#'     \code{\link[tram]{Coxph}}, \code{\link[tram]{Survreg}},
#'     \code{\link[tram]{Lehmann}}. Standard implementations
#'     \code{\link[stats]{lm}}, \code{\link[stats]{glm}},
#'     \code{\link[survival]{survreg}}, \code{\link[survival]{coxph}},
#'     and \code{\link[MASS]{polr}} are also supported. See the corresponding
#'     alias \code{<model_name>ICP}, e.g., \code{\link{PolrICP}} or
#'     \code{?implemented_model_classes}. Models from 'lme4', 'tramME',
#'     'glmnet' and 'mgcv' are also supported.
#' @param verbose Logical, whether output should be verbose (default \code{TRUE}).
#' @param type Character, type of invariance (\code{"residual"} or \code{"wald"});
#'     see \code{Details}.
#' @param test Character, specifies the invariance test to be used when
#'     \code{type = "residual"}. The default is \code{"gcm.test"}. Other
#'     implemented tests are \code{"HSIC"}, \code{"t.test"}, \code{"var.test"},
#'     and \code{"combined"}. Alternatively, a custom function for testing
#'     invariance of the form \code{\(r, e, controls) {...}} can be supplied,
#'     which outputs a list with entry \code{"p.value"}.
#' @param controls Controls for the used tests and the overall procedure,
#'     see \code{\link{dicp_controls}}.
#' @param alpha Level of invariance test, default \code{0.05}.
#' @param ... Further arguments passed to \code{modFUN}.
#' @param baseline_fixed Fixed baseline transformation, see
#'     \code{\link[tramicp]{dicp_controls}}.
#' @param greedy Logical, whether to perform a greedy version of ICP (default is
#'     \code{FALSE}).
#' @param max_size Numeric; maximum support size.
#' @param mandatory A \code{formula} containing mandatory covariates, i.e.,
#'     covariates which by domain knowledge are believed to be parents
#'     of the response or are in another way required for the environment or
#'     model to be valid (for instance, conditionally valid environments or
#'     random effects in a mixed model).
#'
#' @details
#' TRAMICP iterates over all subsets of covariates provided in \code{formula}
#' and performs an invariance test based on the conditional covariance between
#' score residuals and environments in \code{env} (\code{type = "residual"}) or
#' the Wald statistic testing for the presence of main and interaction effects
#' of the environments (\code{type = "wald"}). The algorithm outputs the
#' intersection over all non-rejected sets as an estimate of the causal parents.
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
#' @references
#' Kook, L., Saengkyongam, S., Lundborg, A. R., Hothorn, T., & Peters, J.
#' (2024). Model-based causal feature selection for general response types.
#' Journal of the American Statistical Association, 1-12.
#' \doi{10.1080/01621459.2024.2395588}
#'
#' @export
#'
#' @examples
#' set.seed(12)
#' d <- dgp_dicp(n = 1e3, mod = "binary")
#' dicp(Y ~ X1 + X2 + X3,
#'   data = d, env = ~E, modFUN = "glm",
#'   family = "binomial", type = "wald"
#' )
#'
dicp <- function(
    formula, data, env, modFUN, verbose = TRUE,
    type = c("residual", "wald", "partial"),
    test = "gcm.test", controls = NULL, alpha = 0.05,
    baseline_fixed = TRUE, greedy = FALSE, max_size = NULL,
    mandatory = NULL, ...) {
  call <- match.call()
  type <- match.arg(type)

  ### Preliminary checks
  if (is.character(test)) {
    test <- match.arg(test, .implemented_tests())
  }
  if (is.null(controls)) {
    controls <- dicp_controls(type, test,
      alpha = alpha,
      baseline_fixed = baseline_fixed
    )
  }
  .check_args(formula, data, env, modFUN, type, test)

  ### Process formulae
  tms <- .get_terms(formula)
  resp <- tms$response
  etms <- .get_terms(env)
  me <- setdiff(tms$all, etms$all) # no env in main effects
  if (is.null(max_size)) {
    max_size <- length(me)
  }
  max_size <- min(max_size, length(me))
  ps <- lapply(0:max_size, utils::combn, x = length(me))

  ### Options
  if (verbose && interactive()) {
    pb <- utils::txtProgressBar(min = 0, max = length(ps), style = 3)
  }

  ### Run invariant subset search
  out <- .invariant_subset_search(
    ps = ps, controls = controls, me = me, resp = resp, etms = etms,
    modFUN = modFUN, data = data, greedy = greedy, verbose = verbose, pb = pb,
    mandatory = mandatory, ... = ...
  )

  ### Process output
  tests <- out$tests
  res <- .extract_results(tests)

  ### Extract p-values for each set
  pvals <- structure(res[["pval"]], names = res[["set"]])
  if (!greedy) {
    inv <- try(.inv_set(res, alpha = controls$alpha))
    if (inherits(inv, "try-error")) {
      inv <- "Cannot be computed."
    }
  } else {
    inv <- me[sort(unique(unlist(out$MI)))]
  }

  ### Compute predictor-level p-values
  ipv <- .indiv_pvals(me, pvals, alpha = alpha)
  if (identical(inv, character(0))) {
    inv <- "Empty"
  }
  ### Return
  structure(
    list(
      candidate_causal_predictors = inv, set_pvals = pvals,
      predictor_pvals = ipv, tests = tests, controls = controls,
      call = call
    ),
    class = "dICP", type = type,
    test = controls$test_name, env = env, greedy = greedy
  )
}

#' Return invariant sets
#'
#' @param object Object of class \code{"dICP"}.
#' @param with_pvalues Logical; whether to also return p-values of invariance
#'     tests for the non-rejected sets.
#'
#' @return Returns vector of all non-rejected sets. With
#'    \code{with_pvalues = TRUE}, a named vector of p-values is returned.
#'    Returns \code{named numeric(0)} if there are no invariant sets.
#'
#' @export
invariant_sets <- function(object, with_pvalues = FALSE) {
  pvals <- pvalues(object, "set")
  ret <- pvals[pvals > object$controls$alpha]
  if (with_pvalues) {
    return(ret)
  }
  names(ret)
}

# Invariant subset search -------------------------------------------------

.invariant_subset_search <- function(ps, controls, me, resp, etms, modFUN,
                                     data, greedy, verbose, pb, mandatory,
                                     ...) {
  if (!greedy) {
    ### Run
    tests <- list()
    MI <- NULL
    for (set in seq_along(ps)) {
      if (verbose && interactive()) {
        utils::setTxtProgressBar(pb, set)
      }

      ret <- apply(ps[[set]], 2, controls$type_fun,
        me = me, resp = resp,
        set = set, env = etms, modFUN = modFUN, data = data,
        controls = controls, mandatory = mandatory, ... = ...
      )

      if (set == 1 && controls$stop_if_empty_set_invariant &&
        .get_pvalue(ret[[1]]$test) > controls$alpha) {
        if (verbose && interactive()) {
          message("\nEmpty set is not rejected. Stopping.")
        }
        tests <- ret
        break
      }

      tests <- c(tests, ret)
    }
  } else {
    ### Run
    tests <- list()
    MI <- list()
    lps <- .unlist_once(lapply(ps, \(x) apply(x, 2, \(y) y, simplify = FALSE)))
    for (set in seq_along(lps)) {
      if (verbose && interactive()) {
        utils::setTxtProgressBar(pb, set)
      }

      if (length(MI > 0) && any(unlist(MI) %in% lps[[set]])) {
        tests[[set]] <- .empty_output(me[lps[[set]]], 0)
        next
      }

      ret <- controls$type_fun(
        lps[[set]],
        me = me, resp = resp, set = set, env = etms,
        modFUN = modFUN, data = data, controls = controls,
        mandatory = mandatory, ... = ...
      )
      tpv <- .get_pvalue(ret$test)
      if (!is.nan(tpv) && tpv > controls$alpha) {
        MI <- c(MI, lps[[set]])
      }

      UMI <- sort(unique(unlist(MI)))
      if (length(MI) > 0 && length(UMI) == length(me) && all(UMI == seq_along(me))) {
        break
      }

      tests[[set]] <- ret
    }
  }

  list(tests = tests, MI = MI)
}
