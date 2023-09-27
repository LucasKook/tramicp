
#' Bootstrap stability for TRAMICP
#'
#' @param object Object of class \code{"dICP"}
#' @param B Numeric; number of bootstrap iterations
#' @param size Numeric; size of bootstrap samples
#' @param verbose Logical; print a progress bar (default: \code{FALSE})
#' @param return_all Logical; return all \code{"dICP"} objects (default:
#'     \code{FALSE})
#'
#' @return Table of output sets of candidate causal predictors
#' @export
#'
#' @examples
#' set.seed(12)
#' d <- dgp_dicp(n = 1e3, mod = "binary")
#' res <- glmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'     family = "binomial", test = "cor.test")
#' bootstrap_stability(res, B = 2)
#'
bootstrap_stability <- function(
    object, B = 1e2, size = NULL, verbose = FALSE, return_all = FALSE
) {

  ### Get the call and data
  call <- object$call
  d <- eval(call$data, envir = parent.frame())

  ### Default bootstrap size is nrow(d)
  if (is.null(size))
    size <- nrow(d)

  ### Progress bar
  if (verbose && interactive())
    pb <- utils::txtProgressBar(min = 0, max = B, style = 3)

  ### Repeat updated call B times with boostrap sample of the data
  res <- list()
  for (b in seq_len(B)) {
    if (verbose && interactive())
      utils::setTxtProgressBar(pb, b)
    idx <- sample.int(nrow(d), size, replace = TRUE)
    tmp <- stats::update(object, data = d[idx, ], verbose = FALSE)
    res[[b]] <- if (return_all) tmp else
      tmp[["candidate_causal_predictors", exact = TRUE]]
  }

  if (return_all)
    return(res)

  ### Tabulate the output
  table(unlist(lapply(res, paste0, collapse = "+")))
}

.reduce_size <- function(object) {
  object[["tests"]] <- NULL
  object
}
