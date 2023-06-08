# Helpers

.check_args <- function(formula, data, env, modFUN, type, test) {
  stopifnot("`formula` should be of class 'formula'." = "formula" %in% class(formula))
  stopifnot("`env` should be of class 'formula'." = "formula" %in% class(env))
  stopifnot("`modFUN` should be a function or a character referring to the
            name of a function." = is.function(modFUN) | is.character(modFUN))
  etms <- .get_terms(env)
  if (length(etms$all) != 1 && type %in% c("confint", "mcheck", "wald"))
    stop("Invariance types `\"confint\"`, `\"mcheck\"` and `\"wald\"` are not
         implemented for multivariable environments.")
  if (length(etms$all) != 1 && test %in% c("gcm.test", "cor.test", "t.test"))
    stop("Invariance tests `\"gcm.test\"`, `\"cor.test\"` and `\"t.test\"`
         are not implemented for multivariable environments.")
}

.build_iform <- function(formula) {
  tms <- .get_terms(formula)
  as.formula(paste0("~", paste0(tms$ie, collapse = "+")))
}

.get_terms <- function(formula) {
  atms <- terms(formula)
  tms <- attr(atms, "term.labels")
  resp <- all.vars(formula)[1]
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

#' @import tram
.mod_from_name <- function(mod, prob = c(0.001, 0.999), spec = "correct") {
  if (spec != "link") {
    switch(
      mod,
      "polr" = \(formula, data, ...)
      polr(formula, data, Hess = TRUE, ...),
      # Polr(formula, data, ...),
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
      "polr" = \(formula, data, ...)
      polr(formula, data, method = "probit", Hess = TRUE, ...),
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
  success <- levels(resp)[2]
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
.indiv_pvals <- function(terms, pvals) {
  res <- lapply(terms, \(term) suppressWarnings(
    max(pvals[!grepl(term, names(pvals), fixed = TRUE)], na.rm = TRUE)))
  structure(unlist(res), names = terms)
}

.sub_smooth_terms <- function(tm) {
  gsub("s\\((.*?)\\)", "\\1", tm)
}
