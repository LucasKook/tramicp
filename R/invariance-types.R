
# Type-functions ----------------------------------------------------------

# GOF test
.gof_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, mandatory, ...
) {

  env <- env$all
  mand <- .get_terms(mandatory)$all

  ### Empty set skipped
  if (set == 1)
    return(.empty_output("Empty"))

  ### Set up formula
  tset <- me[tx]
  meff <- .pplus(tset)
  mfm <- reformulate(c(meff, mand), resp)

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
    tx, me, resp, set, env, modFUN, data, controls, mandatory, ...
) {

  ### Checks
  env <- env$all
  stopifnot("`data[[env]]` needs to be a factor." = is.factor(data[[env]]))
  mand <- .get_terms(mandatory)$all

  ### Prepare formula
  tset <- if (set == 1) "1" else me[tx]
  meff <- .pplus(tset)
  mfm <- reformulate(c(meff, mand), resp)

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
.residual_invariance <- function(
    tx, me, resp, set, env, modFUN, data, controls, mandatory, ...
) {

  mand <- .get_terms(mandatory)$all

  ### Prepare formula
  tset <- if (set == "1") 1 else me[tx]
  meff <- .pplus(tset)
  mfm <- reformulate(c(meff, mand), resp)

  ### Fit model
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  ### Test
  r <- matrix(controls$residuals(m), ncol = 1)
  e <- .rm_int(model.matrix(as.formula(env$fml), data = data))
  if (controls$ctest == "gcm.test" & set != "1") # Fit RF for GCM-type test
    e <- .ranger_gcm(e, c(meff, mand), set, data, controls)
  tst <- controls$test_fun(r, e, controls)

  ### Return
  if (set == 1) tset <- "Empty"
  structure(list(set = tset, test = tst, coef = coef(m), tram = m$tram),
            class = "dICPtest")

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
    meff <- .pplus(c(ifelse(controls$baseline_fixed, env, "1"), mand))
    mint <- ""
  } else {
    tset <- me[tx]
    meff <- if (controls$baseline_fixed) .pplus(c(me[tx], env, mand)) else
      .pplus(c(me[tx], mand))
    mint <- .pplus(c(paste0(c(me[tx], mand), ":", env)))
  }

  ### Prepare formula
  mfm <- as.formula(
    paste0(resp, ifelse(
      controls$baseline_fixed, "", paste0("|", .pplus(env))),
      "~", meff, if (mint != "") "+", mint)
  )

  ### Fit
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  if (inherits(m, "tram"))
    m <- as.mlt(m)
  cfs <- names(coef(m))
  tcfs <- union(unlist(sapply(paste0(":", env), \(pat) grep(pat, cfs, value = TRUE))),
                unlist(sapply(env, \(pat) grep(pat, cfs, value = TRUE))))

  ### Test
  tst <- try(summary(
    glht(m, linfct = paste(tcfs, "== 0"), vcov = controls$vcov),
    test = controls$test_fun()), silent = FALSE
  )

  ### Catch failure cases
  if (inherits(tst, "try-error")) {
    return(.empty_output(me[tx], NA))
  }

  ### Return
  if (set == 1) tset <- "Empty"
  tst$set <- tset
  tst$tram <- m$tram
  tst

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
  mfm <- reformulate(c(meff, mint1, mint2), resp)
  m <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))
  cfs <- names(coef(m))
  tcfs <- grep(paste0(":", env), grep(paste0(if (set != 1) ":", mand),
                                      cfs, value = TRUE), value = TRUE)
  tst <- try(summary(glht(m, linfct = paste(tcfs, "== 0"),
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
