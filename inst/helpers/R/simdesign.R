GEN <- function(condition, fixed_objects = NULL) {
  dag <- random_dag(
    nenv = condition$nenv, nanc = condition$nanc, ndec = condition$ndec,
    penv = condition$penv, panc = condition$panc, pdec = condition$pdec,
    cfb = fixed_objects$cfb,
    cfe = .rcfx(nanc + ndec, penv, FALSE, fixed_objects$sde)
  )
  d <- structure(NA, class = "try-error")
  while (inherits(d, "try-error")) {
    d <- dgp_random_dag(n = condition$n, mod = condition$mod,
      K = ifelse(condition$mod == "polr", fixed_objects$polrK, fixed_objects$K),
      rm_censoring = fixed_objects$rmc, standardize = fixed_objects$stdz,
      errDistAnY = fixed_objects$errDistAnY, mixAnY = fixed_objects$mixAnY,
      errDistDeY = fixed_objects$errDistDeY, mixDeY = fixed_objects$mixDeY,
      dag = dag
    )
  }
  d
}

GENFIX <- function(condition, fixed_objects = NULL) {
  d <- structure(NA, class = "try-error")
  while (inherits(d, "try-error")) {
    d <- dgp_random_dag(n = condition$n, mod = condition$mod,
      K = ifelse(condition$mod == "polr", fixed_objects$polrK, fixed_objects$K),
      rm_censoring = fixed_objects$rmc, standardize = fixed_objects$stdz,
      errDistAnY = fixed_objects$errDistAnY, mixAnY = fixed_objects$mixAnY,
      errDistDeY = fixed_objects$errDistDeY, mixDeY = fixed_objects$mixDeY,
      dag = fixed_objects$dags[[condition$dag]]
    )
  }
  d
}

ANA <- function(condition, dat, fixed_objects = NULL) {
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
      if (condition$mod == "polr") {
        ctrls <- dicp_controls(type = ttype, test = ttest,
                               baseline_fixed = fixed_objects$blfix)
        ctrls$vcov <- function(object) {
          cf <- coef(object)
          vcov <- vcov(object)
          vcov[names(cf), names(cf)]
        }
      }
      else
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

SUM <- function(condition, results, fixed_objects = NULL) {
  bind_rows(results) |>
    mutate(
      paY_all = str_split(paY, "\\+"),
      oicp_all = str_split(oicp, "\\+"),
      inv_all = str_split(inv, "\\+"),
      oicp_numpa = map2(paY_all, oicp_all, ~ length(tramicp:::.intersect(.x, .y))),
      inv_numpa = map2(paY_all, inv_all, ~ length(tramicp:::.intersect(.x, .y))),
      oicp_nonpa = map2(oicp_all, paY_all, ~ length(tramicp:::.setdiff(.x, .y))),
      inv_nonpa = map2(inv_all, paY_all, ~ length(tramicp:::.setdiff(.x, .y)))
    ) |> unnest(c(oicp_numpa, inv_numpa, oicp_nonpa, inv_nonpa)) |>
    group_by(npaY, type, test) |>
    reframe(
      mean_oicp_empty = mean(oicp == "empty"),
      mean_inv_empty = mean(inv == "empty"),
      mean_oicp_numpa = mean(oicp_numpa) / npaY,
      mean_inv_numpa = mean(inv_numpa) / npaY,
      mean_oicp_nonpa = (mean(oicp_nonpa) / (length(fixed_objects$preds) - npaY)),
      mean_inv_nonpa = (mean(inv_nonpa) / (length(fixed_objects$preds) - npaY))
    ) |>
    group_by(npaY, type, test) |>
    summarize_if(is.numeric, mean) |>
    group_by(npaY, type, test) |>
    gather(key = "key", value = "val", mean_oicp_empty, mean_inv_empty,
           mean_oicp_numpa, mean_inv_numpa, mean_oicp_nonpa, mean_inv_nonpa) |>
    unite(desc, type, test, npaY, key) |>
    spread(desc, val) |> ungroup()
}

read_sim_results <- function(path) {
  res <- list.files(path, pattern = "results-row-[0-9]+.rds", full.names = TRUE) |>
    map(readRDS)
  bind_rows(map(res, ~ bind_cols(.x$condition, bind_rows(.x$results))))
}

sum_sim_results <- function(res, fobs = NULL) {
  res |>
    group_by(n, mod, dag, nanc, panc, ndec, pdec, penv) |>
    nest() |>
    mutate(summary = map(data, ~ SUM(., data, fobs))) |>
    unnest(c(summary)) |>
    gather(key = "key", value = "val", matches("mean")) |>
    separate(key, into = c("type", "test", "npaY", "mean", "oracle", "metric")) |>
    dplyr::select(-mean) |>
    group_by(n, mod, dag, nanc, panc, ndec, pdec, penv, type, test, npaY, oracle, metric) |>
    reframe(val = mean(val)) |>
    ungroup()
}

vis_sim_results <- function(out, fobs = NULL, nrow = 1, ncol = NULL,
                            lbls = .get_labels(),
                            types = c("wald", "residual", "kci"),
                            mods = c("cotram", "polr", "weibull")) {
  ps <- lapply(unique(out$npaY), \(x) {
    ggplot(out |> dplyr::filter(npaY == x, mod %in% mods, type %in% types) |>
             mutate(test = factor(test, levels = fobs$tests)),
           aes(x = n, y = frac, color = metric, linetype = method)) +
      geom_line() +
      geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.05) +
      facet_grid(test ~ mod, labeller = as_labeller(lbls)) +
      labs(subtitle = paste0("Ground truth: ", x, " parent(s)")) +
      scale_color_discrete(labels = c(
        "empty" = "Empty",
        "numpa" = "Parents",
        "nonpa" = "Non-parents"
      )) +
      scale_linetype(labels = c(
        "inv" = "TRAMICP",
        "oicp" = "Oracle ICP"
      )) +
      labs(x = "sample size", y = "Fraction of output",
           color = "output", linetype = "method") +
      scale_x_log10()
  })
  ggpubr::ggarrange(plotlist = ps, common.legend = TRUE, nrow = nrow, ncol = ncol)
}

.get_labels <- function() {
  # Models
  mods <- c("coxph", "polr", "weibull", "lm", "colr", "cotram", "boxcox", "binary")
  lmods <- c("Coxph", "Polr", "Survreg", "Lm", "Colr", "cotram", "BoxCox", "Binary")
  names(lmods) <- mods

  # Types
  types <- c("wald", "residual", "mcheck", "confint", "kci", "full")
  ltypes <- c("Wald", "Residual", "GOF", "Confint", "KCI", "Wald")
  names(ltypes) <- types

  # Test
  test <- c("wald", "HSIC", "HSIC", "confint", "indep", "KCI", "full")
  ltest <- c("Wald", "HSIC", "HSIC", "Confint", "Permutation", "KCI", "Wald")
  names(ltest) <- test

  c(lmods, ltypes, ltest)
}

#' @importFrom MASS polr
.mod_from_name <- function(mod, prob = c(0.001, 0.999), spec = "correct") {
  if (spec != "link") {
    switch(
      mod,
      "polr" = \(formula, data, ...) {
        res <- try(polr(formula, data, Hess = TRUE, ...))
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