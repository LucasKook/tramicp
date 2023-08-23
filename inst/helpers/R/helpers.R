
#' @importFrom tibble tibble
.binom_test <- function(x, n) {
  tst <- try(binom.test(x = x, n = n))
  if (inherits(tst, "try-error"))
    tst <- list(estimate = NA, conf.int = c(NA, NA))
  tibble(frac = tst$estimate, lwr = tst$conf.int[1], upr = tst$conf.int[2])
}

fmt <- function(x, break.eps = 1e-04, break.middle = 0.01, na.form = "NA",  ...) {
  format1Pval <- function(pv) {
    if (is.na(pv)) {
      na.form
    }
    else if (pv < break.eps) {
      paste("<", format(break.eps, scientific = FALSE))
    }
    else {
      largep <- pv >= break.middle
      format(pv, digits = 1 + largep, nsmall = 1 + largep,
             scientific = FALSE, ...)
    }
  }
  vapply(X = x, FUN = format1Pval, FUN.VALUE = "", USE.NAMES = TRUE)
}

.auc_invariance <- function(tx, me, resp, set, env, modFUN, data, controls,
                            mandatory, ...) {

  env <- env$all
  mand <- tramicp:::.get_terms(mandatory)$all

  ### Empty set treated separately
  if (set == 1) {
    tset <- "1"
    meffx <- tramicp:::.pplus(c("1", mand))
    meff <- tramicp:::.pplus(c(ifelse(controls$baseline_fixed, env, "1"), mand))
    mint <- ""
  } else {
    tset <- me[tx]
    meffx <- tramicp:::.pplus(c(me[tx], mand))
    meff <- if (controls$baseline_fixed) tramicp:::.pplus(c(me[tx], env, mand)) else
      tramicp:::.pplus(c(me[tx], mand))
    mint <- tramicp:::.pplus(c(paste0(c(me[tx], mand), ":", env)))
  }

  ### Prepare formula
  mfx <- reformulate(meffx, resp)
  mfm <- as.formula(
    paste0(resp, ifelse(
      controls$baseline_fixed, "", paste0("|", tramicp:::.pplus(env))),
      "~", meff, if (mint != "") "+", mint)
  )

  ### Fit
  mX <- do.call(modFUN, c(list(formula = mfx, data = data), list(...)))
  mE <- do.call(modFUN, c(list(formula = mfm, data = data), list(...)))

  ### Test
  suppressMessages({
    rX <- roc(data[[resp]], predict(mX, newdata = data), direction = "<")
    rE <- roc(data[[resp]], predict(mE, newdata = data), direction = "<")
    tst <- list(test = roc.test(rX, rE, method = "delong"))
  })

  ### Return
  if (set == 1) tset <- "Empty"
  tst$set <- tset
  tst$tram <- mE$tram
  tst

}

AUCANA <- function(condition, dat, fixed_objects = NULL) {
  mFUN <- simtramicp:::.mod_from_name(condition$mod, spec = fobs$spec)
  paY <- attr(dat, "paY")
  chE <- attr(dat, "chE")
  res <- lapply(fixed_objects$tests, \(ttest) {
    ttype <- fixed_objects$types[ttest] # .get_test(ttype)
    kbw <- if (condition$mod == "polr") 0.01 else 0
    oicp <- attr(dat, "oracle_icp")
    pvals <- if (ttype == "kci") {
      tmp <- cdkci(fixed_objects$resp, fixed_objects$env, fixed_objects$preds,
                   data = dat)
      inv <- attr(tmp, "intersection")
      tmp
    } else {
      if (ttest == "cor.test") {
        ctrls <- dicp_controls("residual", "cor.test", alpha = 0.05,
                               residuals = "residuals")
        ctrls$type_fun <- .auc_invariance
      } else
        ctrls <- NULL
      tmp <- dicp(as.formula(fixed_objects$fml), data = dat,
                  env = reformulate(fixed_objects$env), modFUN = mFUN,
                  type = ttype, test = ttest, weights = rep(1, nrow(dat)),
                  baseline_fixed = fixed_objects$blfix, controls = ctrls,
                  verbose = FALSE)
      inv <- tmp$candidate
      pvalues(tmp, which = "set")
    }
    tibble(paY = paste0(paY, collapse = "+"), chE = paste0(chE, collapse = "+"),
           oicp = paste0(oicp, collapse = "+"), inv = paste0(inv, collapse = "+"),
           npaY = length(attr(dat, "paY")), type = ttype, test = ttest,
           pvals = list(stack(pvals)))
  })
  bind_rows(res)
}

vis <- function(mods = tmods, tests = ttests) {
  tmp <- prdat %>%
    mutate(fwer = as.numeric(fwer > 0)) %>%
    dplyr::filter(mod %in% mods, test %in% tests) %>%
    pivot_longer(c(jaccard, fwer), names_to = "output", values_to = "value") %>%
    group_by(n, mod, test, type, method, output, dag) %>%
    summarize(value = mean(value)) %>%
    ungroup()

  pd <- full_join(tmp %>% dplyr::filter(method == "inv"),
                  tmp %>% dplyr::filter(method == "oicp", test == "wald") %>% mutate(test = "oicp")) %>%
    mutate(output = factor(output, levels = c("jaccard", "fwer"),
                           labels = c("Jaccard", "FWER")),
           mod = factor(mod, levels = mods, labels = names(mods)),
           test = factor(test, levels = c(tests, "oicp"), labels = c(names(tests), "Oracle")))

  pd %>%
    ggplot(aes(x = n, y = value, color = test, fill = test)) +
    stat_summary(geom = "errorbar", width = 0.075) +
    stat_summary(geom = "ribbon", linewidth = NA, alpha = 0.1, show.legend = FALSE) +
    stat_summary(geom = "line", show.legend = FALSE) +
    scale_x_log10() +
    facet_grid(output ~ mod, scale = "free_y") +
    geom_hline(yintercept = 0.05, color = "transparent") +
    labs(color = "Invariance test", linetype = element_blank(), y = "Fraction") +
    scale_color_manual(values = c(cols, "Oracle" = "gray60")) +
    scale_fill_manual(values = c(cols, "Oracle" = "gray60"))
}
