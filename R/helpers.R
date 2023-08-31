# Helpers

.check_args <- function(formula, data, env, modFUN, type, test) {
  stopifnot("`formula` should be a formula." = "formula" %in% class(formula))
  stopifnot("`env` should be a formula." = "formula" %in% class(env))
  stopifnot("`modFUN` should be a function or a character matching a
            function name." = is.function(modFUN) | is.character(modFUN))
  etms <- .get_terms(env)
  if (length(etms$all) != 1 && type %in% c("confint", "mcheck"))
    stop("Invariance types `\"confint\"` and `\"mcheck\"` are not
         implemented for multivariable environments.")
  if (length(etms$all) != 1 && test %in% c("cor.test", "t.test"))
    stop("Invariance tests `\"cor.test\"` and `\"t.test\"`
         are not implemented for multivariable environments.")
}

.build_iform <- function(formula) {
  tms <- .get_terms(formula)
  as.formula(paste0("~", paste0(tms$ie, collapse = "+")))
}

.get_terms <- function(formula) {
  if (is.null(formula))
    return(NULL)
  atms <- terms(formula)
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

#' @importFrom MASS polr
#' @import tram
.mod_from_name <- function(mod, prob = c(0.001, 0.999), spec = "correct") {
  if (spec != "link") {
    switch(
      mod,
      "polr" = \(formula, data, ...) {
        res <- try(polr(formula, data, Hess = TRUE, ...), silent = TRUE)
        if (inherits(res, "try-error"))
          res <- Polr(formula, data, ...)
        res
      },
      "weibull" = \(formula, data, ...)
      Survreg(formula, data, ..., prob = prob),
      "lm" = \(formula, data, ...)
      Lm(formula, data, ..., prob = prob),
      "coxph" = \(formula, data, ...)
      Coxph(formula, data, prob = prob, ...),
      "colr" = \(formula, data, ...)
      Colr(formula, data, ..., prob = prob),
      "boxcox" = \(formula, data, ...)
      BoxCox(formula, data, ..., prob = prob),
      "cotram" = \(formula, data, ...)
      cotram::cotram(formula, data, ..., log_first = FALSE, prob = prob[2]),
      "binary" = \(formula, data, ...) {
        m <- stats::glm(formula, data, ..., family = binomial(link = "logit"))
        structure(m, class = c("binglm", class(m)))
      }
    )
  } else {
    switch(
      mod,
      "polr" = \(formula, data, ...) {
        res <- try(polr(formula, data, method = "probit", Hess = TRUE, ...))
        if (inherits(res, "try-error"))
          res <- Polr(formula, data, method = "probit", ...)
        res
      },
      "weibull" = \(formula, data, ...)
      Survreg(formula, data, dist = "loglogistic", ..., prob = prob),
      "lm" = \(formula, data, ...)
      Survreg(formula, data, dist = "logistic", ..., prob = prob),
      "coxph" = \(formula, data, ...)
      Colr(formula, data, prob = prob, ...),
      "colr" = \(formula, data, ...)
      Coxph(formula, data, ..., prob = prob),
      "boxcox" = \(formula, data, ...)
      Colr(formula, data, ..., prob = prob),
      "cotram" = \(formula, data, ...)
      cotram::cotram(formula, data, ..., method = "probit", log_first = FALSE,
                     prob = prob[2]),
      "binary" = \(formula, data, ...) {
        m <- stats::glm(formula, data, ..., family = binomial(link = "probit"))
        structure(m, class = c("binglm", class(m)))
      }
    )
  }
}

.get_test <- function(type) {
  switch(type, "residual" = "HSIC", "mcheck" = "HSIC", "wald" = "wald",
         "kci" = "KCI")
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

.ci <- function(alpha, ms, nenv) {
  # TODO: Extend to more than two environments
  talpha <- alpha / nenv
  ells <- lapply(ms, \(m) list(c = coef(m), C = vcov(m)))
  nempty <- intersect.ellipses(ells[[1]]$c, ells[[2]]$c, ells[[1]]$C,
                               ells[[2]]$C, talpha)
  if (!nempty) return(-talpha) else talpha
}

# Taken from https://github.com/runesen/icph.git
intersect.ellipses <- function(c1, c2, C1, C2, alpha){

  if (alpha == 1)
    return(FALSE)

  if (alpha == 0)
    return(TRUE)

  m <- length(c1)
  q <- sqrt(qchisq(1 - alpha, df = m))

  if (m == 1) {
    intersect <- c(abs(c1 - c2) < (sqrt(C1) + sqrt(C2)) * q)
  } else{

    e1 <- eigen(C1)
    d1 <- e1$values
    d1 <- makePositive(d1)
    d1 <- sqrt(d1)
    U1 <- e1$vectors

    c <- 1/q * diag(1/d1) %*% t(U1) %*%(c2-c1)
    C <- diag(1/d1) %*% t(U1) %*% C2 %*% U1 %*% diag(1/d1)

    e <- eigen(C)
    U <- e$vectors
    d <- e$values
    d <- makePositive(d,silent=TRUE)
    d <- sqrt(d)

    y <- -t(U) %*% c
    y <- abs(y) # 0 expressed in coordinate system of l and rotated to the first quadrant

    if (sum((y / d) ^ 2) <= 1) {
      # y inside the ellipse
      intersect <- TRUE
    } else{
      # Newton-Rhapson iterations
      # f goes monotone, quadratically to -1, so sure and fast convergence
      f <- function(t) sum((d * y / (t + d ^ 2)) ^ 2) - 1
      df <- function(t) - 2 * sum((y * d) ^ 2 / (t + d ^ 2) ^ 3)

      t0 <- 0
      ft0 <- f(t0)

      while (ft0 > 1e-4) {
        t0 <- t0 - f(t0) / df(t0)
        ft0 <- f(t0)
      }

      x0 <- y * d ^ 2 / (d ^ 2 + t0) # projection of y onto (c,C)
      dist <- sqrt(sum((y - x0) ^ 2))
      intersect <- (dist < 1)
    }
  }
  intersect
}

# Taken from https://github.com/runesen/icph.git
makePositive <- function(v, silent = TRUE) {
  w <- which(v < 10 ^ (-14))
  if (length(w) > 0 &
      !silent)
    warning("Some eigenvalues are below 10^(-14) and will therefore be regularized")
  for (ww in w) {
    if (v[ww] < 0) {
      v[ww] <- v[ww] - 2 * v[ww] + 10 ^ (-14)
    } else {
      v[ww] <- v[ww] + 10 ^ (-14)
    }
  }
  return(v)
}

#' @method residuals binglm
residuals.binglm <- function(object, ...) {
  resp <- model.response(model.frame(object))
  success <- if (is.factor(resp)) levels(resp)[2] else sort(unique(resp))[2]
    y <- c(0, 1)[1 + as.numeric(resp == success)]
  y - predict(object, type = "response")
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
  structure(list(set = set, test = list("p.value" = pv, test = NA), coef = NA, tram = NA),
            class = "dICPtest")
}

.compute_residuals <- function(formula, data, meff, mand, set, controls, modFUN, env, ...) {
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
      e <- .rm_int(model.matrix(as.formula(env$fml), data = data))
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
    e <- .rm_int(model.matrix(as.formula(env$fml), data = data))
    if (controls$ctest == "gcm.test" & set != "1") # Fit RF for GCM-type test
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
  mfe <- reformulate(sapply(meff, .sub_smooth_terms), "e")
  if (dprob <- (length(unique(e)) == 2)) {
    fe <- as.factor(e)
    mfe <- update(mfe, fe ~ .)
  }
  rf <- ranger(mfe, data = data, probability = dprob)
  if (dprob) e - predict(rf, data = data)$predictions[, 2] else
    e - predict(rf, data = data)
}

#' @importFrom ranger ranger
.ranger_gcm_cf <- function(e, meff, set, data, controls, etest, test) {
  if (NCOL(e) != 1L) {
    resids <- apply(e, 2, .ranger_gcm_cf, meff = meff, set = set, data = data,
                    controls = controls, etest = etest, test = test)
    return(resids)
  }
  mfe <- reformulate(sapply(meff, .sub_smooth_terms), "e")
  if (dprob <- (length(unique(e)) == 2)) {
    fe <- as.factor(e)
    mfe <- update(mfe, fe ~ .)
  }
  rf <- ranger(mfe, data = data, probability = dprob)
  if (dprob) etest - predict(rf, data = test)$predictions[, 2] else
    etest - predict(rf, data = test)
}

.pplus <- function(terms) {
  paste0(terms, collapse = "+")
}

RANGER <- function(formula, data, ...) {
  response <- model.response(model.frame(formula, data))
  is_ordered <- is.ordered(response)
  is_binary <- is.factor(response) && !is_ordered
  tms <- .get_terms(formula)
  if (identical(tms$me, character(0))) {
    if (is_ordered)
      return(polr(formula, data))
    else if (is_binary)
      return(glm(formula, data, family = "binomial"))
    else
      return(lm(formula, data))
  }
  ret <- ranger(formula, data, probability = is_binary | is_ordered, ...)
  ret$data <- data
  ret$response <- if (is_binary) as.numeric(response) - 1 else if (is_ordered)
    as.numeric(response) else response
  ret$is_binary <- is_binary
  ret$is_ordered <- is_ordered
  ret
}

#' @exportS3Method residuals ranger
residuals.ranger <- function(object, ...) {
  if ("polr" %in% class(object))
    return(residuals.polr(object))
  else if ("glm" %in% class(object))
    return(residuals.binglm(object))
  else if ("lm" %in% class(object))
    return(residuals(object))
  preds <- predict(object, data = object$data)$predictions
  if (object$is_ordered)
    preds <- preds %*% seq_len(ncol(preds))
  if (object$is_binary)
    preds <- preds[, 2]
  object$response - preds
}

.intersect <- function(x, y) {
  ret <- try(intersect(x, y))
  if (inherits(ret, "try-error"))
    return(character(0))
  ret
}
