
# Type-functions ----------------------------------------------------------

# Residual invariance
.residual_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, mandatory, ...
) {

  mand <- .get_terms(mandatory)$all

  ### Prepare formula
  tset <- if (set == "1") 1 else me[tx]
  meff <- .pplus(tset)
  mfm <- stats::reformulate(c(meff, mand), resp)

  ### (Cross-) fit models
  resids <- .compute_residuals(mfm, data, tset, mand, set, controls, modFUN,
                               env, ...)
  m <- resids$m
  r <- resids$r
  e <- resids$e
  tst <- controls$test_fun(r, e, controls)

  # plot(r ~ e, main = paste0(tset, collapse = "+"), col = rgb(.1, .1, .1, .1))
  # legend("top", ifelse(tst$p.value > 0.05, "accept", "reject"), bty = "n")
  # abline(lm(r ~ e), lwd = 1.5)

  ### Return
  if (set == 1) tset <- "Empty"
  structure(list(set = tset, test = tst, coef = stats::coef(m), tram = m$tram,
                 rYX = r, rEX = e), class = "dICPtest")

}

# Wald invariance
.wald_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, mandatory, ...
) {

  env <- env$all
  mand <- .get_terms(mandatory)$all

  ### Empty set treated separately
  if (set == 1) {
    tset <- "1"
    meff <- .pplus(c(if (controls$baseline_fixed) env else "1", mand))
    mint <- ""
  } else {
    tset <- me[tx]
    meff <- if (controls$baseline_fixed) .pplus(c(me[tx], env, mand)) else
      .pplus(c(me[tx], mand))
    mint <- .pplus(c(paste0(c(me[tx], mand), ":", env)))
  }

  if (!controls$wald_test_interactions)
    mint <- ""

  ### Prepare formula
  mfm <- stats::as.formula(
    paste0(resp, ifelse(
      controls$baseline_fixed, "", paste0("|", .pplus(env))),
      "~", meff, if (mint != "") "+", mint)
  )

  ### Fit
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  if (inherits(m, "tram"))
    m <- mlt::as.mlt(m)
  cfs <- names(stats::coef(m))
  tcfs <- .env_coefs(cfs, env, data)

  ### Test
  tst <- try(summary(
    multcomp::glht(m, linfct = paste(tcfs, "== 0"), vcov = controls$vcov),
    test = controls$test_fun()), silent = FALSE
  )

  ### Catch failure cases
  if (inherits(tst, "try-error")) {
    return(.empty_output(if (set == 1) "Empty" else me[tx], NA))
  }

  ### Return
  if (set == 1) tset <- "Empty"
  tst$set <- tset
  tst$tram <- m$tram
  tst

}

# Coefficients corresponding to environment main effects and interactions
.env_coefs <- function(cfs, env, data) {
  ecfs <- colnames(.rm_int(stats::model.matrix(stats::reformulate(env), data)))
  ecfs <- c(ecfs, unlist(lapply(env, \(e) if (is.factor(data[[e]]))
    paste0(e, levels(data[[e]])))))
  cfs[vapply(strsplit(cfs, ":", fixed = TRUE), \(cf) any(cf %in% ecfs),
             logical(1))]
}

# Partial invariance
.partial_invariance <- function(
    tx, me, resp, set, env, modFUN, data, trt, controls, mandatory, ...
) {

  mand <- .get_terms(mandatory)$all
  env <- env$all

  if (set == 1) {
    tset <- "1"
    meff <- .pplus(c(mand, env))
    mint1 <- .pplus(paste0(mand, ":", env))
    mint2 <- NULL
  } else {
    tset <- me[tx]
    meff <- .pplus(c(me[tx], mand, env))
    mint1 <- .pplus(c(paste0(mand, ":", tset), paste0(env, ":", tset)))
    mint2 <- .pplus(paste0(mand, ":", env, ":", tset))
  }
  mfm <- stats::reformulate(c(meff, mint1, mint2), resp)
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  cfs <- names(stats::coef(m))
  tcfs <- grep(paste0(":", env), grep(paste0(if (set != 1) ":", mand),
                                      cfs, value = TRUE), value = TRUE)
  tst <- try(summary(multcomp::glht(m, linfct = paste(tcfs, "== 0"),
                                    vcov = controls$vcov),
                     test = controls$test_fun()), silent = FALSE)
  if (inherits(tst, "try-error")) {
    empty_res <- list(test = list(p.value = NA), set = me[tx])
    return(empty_res)
  }

  if (set == 1) tset <- "Empty"
  tst$set <- tset
  tst

}
