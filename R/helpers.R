# Helpers

.check_args <- function(formula, data, env, modFUN, type, test) {
  stopifnot("`formula` should be a formula." = "formula" %in% class(formula))
  stopifnot("`env` should be a formula." = "formula" %in% class(env))
  stopifnot("`modFUN` should be a function or a character matching a
            function name." = is.function(modFUN) | is.character(modFUN))
  etms <- .get_terms(env)
  if (length(etms$all) != 1 && test %in% c("cor.test", "t.test"))
    stop("Invariance tests `\"cor.test\"` and `\"t.test\"`
         are not implemented for multivariable environments.")
}

.build_iform <- function(formula) {
  tms <- .get_terms(formula)
  stats::as.formula(paste0("~", paste0(tms$ie, collapse = "+")))
}

.get_terms <- function(formula) {
  if (is.null(formula))
    return(NULL)
  atms <- stats::terms(formula)
  tms <- attr(atms, "term.labels")
  resp <- all.vars(formula)[1]
  ridx <- grep("|", tms, fixed = TRUE)
  tms[ridx] <- paste0("(", tms[ridx], ")")
  ie <- grep(":", tms, value = TRUE)
  me <- grep(":", tms, value = TRUE, invert = TRUE)
  list(all = tms, me = me, ie = ie, response = resp, terms = atms, fml = formula)
}

.rm_int <- function(x) {
  if (all(x[, 1] == 1))
    return(x[, -1L, drop = FALSE])
  return(x)
}

.inv_set <- function(x, alpha = 0.05) {
  sets <- strsplit(x[["set"]], split = "\\+")
  idx <- which(x[["pval"]] > alpha)
  if (identical(idx, integer(0)))
    return(character(0))
  if (length(idx) == 1)
    return(sets[[idx]])
  Reduce(intersect, sets[idx])
}

.mod_from_name <- function(mod, prob = c(0.001, 0.999), spec = "correct") {
  if (spec != "link") {
    switch(
      mod,
      "polr" = \(formula, data, ...)
      tram::Polr(formula, data, ...),
      "weibull" = \(formula, data, ...)
      tram::Survreg(formula, data, ..., prob = prob),
      "lm" = \(formula, data, ...)
      tram::Lm(formula, data, ..., prob = prob),
      "coxph" = \(formula, data, ...)
      tram::Coxph(formula, data, prob = prob, ...),
      "colr" = \(formula, data, ...)
      tram::Colr(formula, data, ..., prob = prob),
      "boxcox" = \(formula, data, ...)
      tram::BoxCox(formula, data, ..., prob = prob),
      "cotram" = \(formula, data, ...)
      cotram::cotram(formula, data, ..., log_first = FALSE, prob = prob[2]),
      "binary" = \(formula, data, ...) {
        m <- stats::glm(formula, data, ..., family = stats::binomial(link = "logit"))
        structure(m, class = c("binglm", class(m)))
      }
    )
  } else {
    switch(
      mod,
      "polr" = \(formula, data, ...)
      tram::Polr(formula, data, method = "probit", ...),
      "weibull" = \(formula, data, ...)
      tram::Survreg(formula, data, dist = "loglogistic", ..., prob = prob),
      "lm" = \(formula, data, ...)
      tram::Survreg(formula, data, dist = "logistic", ..., prob = prob),
      "coxph" = \(formula, data, ...)
      tram::Colr(formula, data, prob = prob, ...),
      "colr" = \(formula, data, ...)
      tram::Coxph(formula, data, ..., prob = prob),
      "boxcox" = \(formula, data, ...)
      tram::Colr(formula, data, ..., prob = prob),
      "cotram" = \(formula, data, ...)
      cotram::cotram(formula, data, ..., method = "probit", log_first = FALSE,
                     prob = prob[2]),
      "binary" = \(formula, data, ...) {
        m <- stats::glm(formula, data, ..., family = stats::binomial(link = "probit"))
        structure(m, class = c("binglm", class(m)))
      }
    )
  }
}

.get_test <- function(type) {
  switch(type, "residual" = "HSIC", "wald" = "wald", "kci" = "KCI")
}

.setdiff <- function(x, y) {
  x <- x[x != "empty"]
  y <- y[y != "empty"]
  ret <- try(setdiff(x, y))
  if (inherits(ret, "try-error"))
    return(character(0))
  ret
}

.intersect_intervals <- function(v1, v2) {
  if (identical(v1, numeric(0)) | identical(v2, numeric(0)))
    return(numeric(0))
  if (v2[1] > v1[2] | v1[1] > v2[2])
    return(numeric(0))
  else
    c(max(v1[1], v2[1]), min(v1[2], v2[2]))
}

intersect_intervals <- function(...) {
  dots <- list(...)
  if (length(dots) == 1)
    return(dots[[1]])
  else .intersect_intervals(dots[[1]], do.call(intersect_intervals, dots[-1]))
}

.df2vecs <- function(df) {
  apply(df, 1, c, simplify = FALSE)
}

#' @method residuals binglm
residuals.binglm <- function(object, ...) {
  resp <- stats::model.response(stats::model.frame(object))
  success <- if (is.factor(resp)) levels(resp)[2] else sort(unique(resp))[2]
    y <- c(0, 1)[1 + as.numeric(resp == success)]
  y - stats::predict(object, type = "response")
}

.check_depth <- function(x) {
  if (is.null(x)) {
    0L
  }
  else if (is.atomic(x)) {
    1L
  }
  else if (is.list(x)) {
    depths <- as.integer(unlist(lapply(x, .check_depth)))
    1L + max(depths, 0L)
  }
  else {
    stop("`x` must be a vector")
  }
}

.unlist_once <- function(x) {
  if (is.data.frame(x)) {
    return(lapply(unname(unlist(x)), function(x) c(x)))
  } else if (.check_depth(x) <= 2L) {
    return(x)
  } else {
    unlist(x, recursive = FALSE)
  }
}

.extract_results <- function(res) {
  ret <- unlist(lapply(res, \(x) {
    ret <- .get_pvalue(x$test)
    names(ret) <- paste(x$set, collapse = "+")
    ret
  }))
  data.frame(set = names(ret), pval = unname(ret))
}

.get_pvalue <- function(x) {
  ret <- x[["p.value"]]
  if (inherits(x, "gtest"))
    ret <- c(x[["pvalue"]])
  if (is.null(ret))
    warning("supplied test has no entry called `p.value`")
  ret
}

# Pvalues for individual predictors being a causal parent
.indiv_pvals <- function(terms, pvals, alpha) {
  res <- lapply(terms, \(term) suppressWarnings(
    max(pvals[!grepl(term, names(pvals), fixed = TRUE)], na.rm = TRUE)))
  ret <- structure(unlist(res), names = terms)
  if (all(ret < alpha))
    ret[] <- 1
  ret
}

.sub_smooth_terms <- function(tm) {
  gsub("\\w+\\(([^,)]+).*\\)", "\\1", tm)
}

.empty_output <- function(set, pv = 0) {
  structure(list(set = set, test = list("p.value" = pv, test = NA), coef = NA,
                 tram = NA), class = "dICPtest")
}

.compute_residuals <- function(formula, data, meff, mand, set, controls, modFUN,
                               env, ...) {
  if (controls$crossfit) {
    idx <- sample.int(nrow(data), floor(nrow(data) / 2))
    rs <- numeric(nrow(data))
    es <- numeric(nrow(data))
    for (sgn in c(-1, 1)) {
      train <- data[sgn * idx, ]
      test <- data[- sgn * idx, ]
      m <- do.call(modFUN, c(list(formula = formula, data = train), list(...)))

      ### Test
      r <- matrix(controls$residuals(m, newdata = test), ncol = 1)
      e <- .rm_int(stats::model.matrix(stats::as.formula(env$fml), data = data))
      if (controls$ctest == "gcm.test" & set != "1") # Fit RF for GCM-type test
        e <- .ranger_gcm_cf(e[sgn * idx, ], c(meff, mand), set, train, controls,
                            e[-sgn * idx, ], test)
      else e <- e[- sgn * idx]
      rs[- sgn * idx] <- r
      es[- sgn * idx] <- e
    }
    r <- rs
    e <- es
  } else {
    m <- do.call(modFUN, c(list(formula = formula, data = data), list(...)))

    ### Test
    r <- matrix(controls$residuals(m), ncol = 1)
    e <- .rm_int(stats::model.matrix(stats::as.formula(env$fml), data = data))
    if (controls$ctest == "gcm.test") # Fit RF for GCM-type test
      e <- .ranger_gcm(e, c(meff, mand), set, data, controls)
  }

  list(r = r, e = e, m = m)
}

.ranger_gcm <- function(e, meff, set, data, controls) {
  if (NCOL(e) != 1L) {
    resids <- apply(e, 2, .ranger_gcm, meff = meff, set = set, data = data,
                    controls = controls)
    return(resids)
  }
  data[["ranger_e"]] <- e
  mfe <- stats::reformulate(sapply(meff, .sub_smooth_terms), "ranger_e")
  rf <- RANGER(mfe, data = data)
  stats::residuals(rf)
}

.ranger_gcm_cf <- function(e, meff, set, data, controls, etest, test) {
  if (NCOL(e) != 1L) {
    resids <- apply(e, 2, .ranger_gcm_cf, meff = meff, set = set, data = data,
                    controls = controls, etest = etest, test = test)
    return(resids)
  }
  data[["ranger_e"]] <- e
  mfe <- stats::reformulate(sapply(meff, .sub_smooth_terms), "ranger_e")
  rf <- RANGER(mfe, data = data)
  stats::residuals(rf, newdata = test, newy = etest)
}

.pplus <- function(terms) {
  paste0(terms, collapse = "+")
}

RANGER <- function(formula, data, ...) {
  response <- stats::model.response(stats::model.frame(formula, data))
  is_ordered <- is.ordered(response)
  is_binary <- is.factor(response) && !is_ordered
  tms <- .get_terms(formula)
  resp <- if (is_binary) as.numeric(response) - 1 else if (is_ordered)
    as.numeric(response) else response
  tmp <- list(data = data, response = resp, is_binary = is_binary,
              is_ordered = is_ordered)
  if (identical(tms$me, character(0))) {
    if (is_binary)
      return(structure(c(list(mean = mean(as.numeric(response) - 1)), tmp),
                       class = "ranger"))
    else
      return(structure(c(list(mean = mean(as.numeric(response))), tmp),
                       class = "ranger"))
  }
  ret <- ranger::ranger(formula, data, probability = is_binary | is_ordered,
                        ...)
  structure(c(ret, tmp), class = "ranger")
}

#' @exportS3Method residuals ranger
residuals.ranger <- function(object, newdata = NULL, newy = NULL, ...) {
  if (is.null(newdata))
    newdata <- object$data
  if (!is.null(newy))
    newy <- if (object$is_binary) as.numeric(newy) - 1 else if (object$is_ordered)
      as.numeric(newy) else newy
  if (is.null(newy))
    newy <- object$response
  if (!is.null(object$mean))
    return(newy - object$mean)
  preds <- stats::predict(object, data = newdata)$predictions
  if (object$is_ordered)
    preds <- preds %*% seq_len(ncol(preds))
  if (object$is_binary)
    preds <- preds[, 2]
  unname(newy - preds)
}

.intersect <- function(x, y) {
  ret <- try(intersect(x, y))
  if (inherits(ret, "try-error"))
    return(character(0))
  ret
}
