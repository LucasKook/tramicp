# S3 methods

#' Distributional ICP test
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

  cat("\nDistributional Invariant Causal Prediction\n")
  mod <- object[["tests"]][[2]]$model$tram
  if (is.null(mod)) mod <- object[["tests"]][[2]]$tram
  cat(mod, "\n")
  cat("\nCall: ")
  print(attr(object, "call"))
  cat("\n")
  # cat(" Invariance type:", attr(object, "type"), "\n")
  cat(" Invariance test:", attr(object, "test"), "\n")
  # cat(" Environment:", deparse(attr(object, "env")), "\n")
  # if (attr(object, "type") == "partial")
    # cat("Treatment:", attr(object, "trt"), "\n")
  if (print_all) {
    cat("\n Tested sets:\n")
    print(round(object[["pvals"]], digits = digits))
  }
  cat("\n Predictor p-values:\n")
  print(round(object[["ipvals"]], digits = digits))
  cat("\n Set of plausible causal predictors:", object[["inv"]], "\n")
  return(invisible(NULL))

}
