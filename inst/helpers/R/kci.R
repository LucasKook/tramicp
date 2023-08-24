
#' Conditional independence testing for causal discovery in trams
#'
#' @param Y Response
#' @param E Environments
#' @param X Predictors
#' @param data Data
#' @param coin Conditional Independence Test
#' @param ... Further arguments supplied to \code{coin}
#' @param verbose Toggle berbose output; defaults to \code{TRUE}
#'
#' @return See KCI
#' @export
#'
#' @importFrom RCIT RCIT
#' @importFrom CondIndTests KCI
#' @importFrom GeneralisedCovarianceMeasure gcm.test
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' d <- dgp_random_dag(n = 5e1, mod = "colr", dag = random_dag(nanc = 2,
#'     ndec = 2, panc = 1, pdec = 1, penv = 1))
#' cdkci("Y", "E", paste0("X", 1:4), data = d, coin = CondIndTests::CondIndTest)
#' cdkci("Y", "E", paste0("X", 1:4), data = d)
#' }
#'
#'
cdkci <- function(Y, E, X, data, coin = RCIT, verbose = FALSE, ...) {

  ps <- lapply(0:length(X), combn, x = length(X))

  if (length(E) == 1 & length(unique(data[, E])) <= 5)
    data[, E] <- factor(data[, E])

  if (inherits(data[, Y], "response"))
    data[, Y] <- .R2vec(data[, Y])

  # Options
  if (verbose)
    pb <- txtProgressBar(min = 0, max = length(ps), style = 3)

  tests <- c()
  for (npred in seq_along(ps)) {

    if (verbose)
      setTxtProgressBar(pb, npred)

    ret <- apply(ps[[npred]], 2, function(set) {
      tset <- X[set]
      if (identical(tset, character(0))) {
        tset <- "Empty"
        pval <- 0
      } else {
        tst <- try(coin(as.numeric(data[, Y]), data[, E],
                        data[, tset], ... = ...), silent = TRUE)
        if (inherits(tst, "try-error"))
          tst <- coin(as.matrix(as.double(data[, Y])),
                      as.matrix(as.double(data[, E])),
                      as.matrix(data[, tset]), ... = ...)
        pval <- tst[["pvalue"]]
        if (is.null(pval)) pval <- tst[["p.value"]]
        if (is.null(pval)) pval <- tst[["p"]]

      }
      structure(pval, names = paste(tset, collapse = "+"))
    }, simplify = FALSE)
    tests <- c(tests, unlist(ret))
  }
  ints <- tramicp:::.inv_set(data.frame(set = names(tests), pval = unname(tests)))
  structure(tests, intersection = ifelse(identical(ints, character(0)), "Empty", ints))
}

.R2vec <- as.double
