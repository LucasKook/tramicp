# S3 methods

#' ICP test
#' @noRd
#'
#' @param x Object of class \code{"dICPtest"}
#' @param ... Currently ignored
#'
#' @return Numeric vector of p-values from invariance tests
#' @method print dICPtest
#' @exportS3Method
#'
print.dICPtest <- function(x, ...) {

  cat("\n", "Set:", x$set, "\n\tp-value:",
      .get_pvalue(x$test))

}

#' Distributional ICP
#' @noRd
#'
#' @param x Object of class \code{"dICP"}
#' @param ... Currently ignored
#'
#' @return Invisible \code{NULL}
#' @method print dICP
#' @exportS3Method
#'
print.dICP <- function(x, ...) {

  summary(x)
  return(invisible(NULL))

}

#' Summary method dICP
#' @noRd
#'
#' @param object Object of class \code{"dICP"}
#' @param digits Numeric, number of digits to print, default is 3
#' @param print_all Print p-values for all tested sets
#' @param ... Currently ignored
#'
#' @return Invisible \code{NULL}
#' @exportS3Method
#'
#' @method summary dICP
#'
summary.dICP <- function(object, print_all = FALSE, digits = 3, ...) {

  ttitle <- paste(
    "\n", ifelse(attr(object, "greedy"), "Greedy model-based", "Model-based"),
    "Invariant Causal Prediction\n"
  )
  cat(ttitle)
  if (is.null(mod <- object[["tests"]][[1]]$model$tram))
    mod <- object[["tests"]][[1]]$tram
  if (!is.null(mod))
    cat(mod, "\n")
  if (length(deparse(tcall <- object$call)) < 10) {
    cat("\nCall: ")
    print(tcall)
    cat("\n")
  }
  # cat(" Invariance type:", attr(object, "type"), "\n")
  cat(" Invariance test:", attr(object, "test"), "\n")
  # cat(" Environment:", deparse(attr(object, "env")), "\n")
  # if (attr(object, "type") == "partial")
    # cat("Treatment:", attr(object, "trt"), "\n")
  if (print_all) {
    cat("\n Tested sets:\n")
    print(round(object[["set_pvals"]], digits = digits))
  }
  cat("\n Predictor p-values:\n")
  print(round(object[["predictor_pvals"]], digits = digits))
  cat("\n Set of plausible causal predictors:",
      object[["candidate_causal_predictors"]], "\n")
  return(invisible(NULL))

}

#' Extract set and predictor p-values from tramicp outputs
#' @param object Object of class \code{'dicp'}
#' @param which Which p-values to return, \code{"predictor"} returns p-values
#'     for individual predictors, \code{"set"} for each subset of the predictors,
#'     \code{"all"} returns a list of both
#'
#' @return Numeric vector (or list in case \code{which = "all"}) of p-values
#'
#' @details Predictor p-values are computed from the set p-values as follows:
#'     For each predictor j as the largest p-value of all sets not containing j.
#'
#' @examples
#' set.seed(123)
#' d <- dgp_dicp(n = 1e3, mod = "polr")
#' res <- polrICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, type = "wald")
#' pvalues(res, which = "predictor")
#' pvalues(res, which = "set")
#' pvalues(res, which = "all")
#'
#' @export
pvalues <- function(object, which = c("predictor", "set", "all")) {
  which <- match.arg(which)
  switch(which, "predictor" = object$predictor_pvals, "set" = object$set_pvals,
         "all" = object[c("predictor_pvals", "set_pvals")])
}
