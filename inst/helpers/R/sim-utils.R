simall <- function(d = dgp_random_dag(n), fml = Y ~ X1 + X2 + X3,
                   modFUN = Polr, pfml = A ~ X1 + X2, type = "residual",
                   test = "HSIC", wtype = c("none", "oracle", "estimate"),
                   blfix = TRUE, pFUN = .policy_clipped, resp = "Y", env = "E",
                   preds = paste0("X", 1:3), kbw = 0) {

  wtype <- match.arg(wtype)

  ps <- switch(
    wtype,
    "none" = rep(1, nrow(d)),
    "oracle" = sapply(1:n, function(x) pFUN(d$X1[x], d$X2[x])[as.numeric(d$A[x])]),
    "estimate" = {
      ps <- predict(glm(pfml, data = d, family = binomial), type = "response")
      apply(cbind(1 - ps, ps) * model.matrix(~ 0 + A, data = d), 1, max)
    }
  )

  w <- 1 / ps

  res <- if (type == "kci") {
    cdkci(resp, env, preds, data = d, gammaApprox = TRUE, GP = FALSE, width = kbw)
  } else if (type %in% c("residual", "wald", "mcheck")) {
    dicp(fml, data = d, env = "E", modFUN = modFUN, trt = "A", type = type,
         test = test, weights = w, baseline_fixed = blfix, verbose = FALSE)[["pvals"]]
  } else stop("Not implemented.")
  res
}

sim <- function(n = 1e3, dgp = dgp_dicp, fml = Y ~ X1 + X2 + X3,
                modFUN = Polr, pfml = A ~ X1 + X2, type = "residual",
                test = "HSIC", wtype = c("none", "oracle", "estimate"),
                blfix = TRUE, pFUN = .policy_clipped, resp = "Y", env = "E",
                preds = paste0("X", 1:3), kbw = 0) {

  wtype <- match.arg(wtype)
  d <- dgp(n = n)

  ps <- switch(
    wtype,
    "none" = rep(1, nrow(d)),
    "oracle" = sapply(1:n, function(x) pFUN(d$X1[x], d$X2[x])[as.numeric(d$A[x])]),
    "estimate" = {
      ps <- predict(glm(pfml, data = d, family = binomial), type = "response")
      apply(cbind(1 - ps, ps) * model.matrix(~ 0 + A, data = d), 1, max)
    }
  )

  w <- 1 / ps

  res <- if (type == "kci") {
    cdkci(resp, env, preds, data = d, gammaApprox = TRUE, GP = FALSE,
          width = kbw)
  } else if (type %in% c("residual", "wald", "mcheck")) {
    dicp(fml, data = d, env = "E", modFUN = modFUN, trt = "A", type = type,
         test = test, weights = w, baseline_fixed = blfix, verbose = FALSE)[["pvals"]]
  } else stop("Not implemented.")
  res
}

#' @importFrom dplyr bind_rows mutate starts_with
#' @importFrom tidyr gather
.summarize_results <- function(res, n) {
  dplyr::bind_rows(res) |>
  dplyr::mutate(n = n) |>
  tidyr::gather("set", "frac_accept", dplyr::starts_with("X"))
}

.visualize_results <- function(out) {
  ggplot2::ggplot(out, ggplot2::aes(x = ordered(n), y = frac_accept,
                                    color = set, group = set)) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0.95, lty = 3) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(x = "sample size", y = "proportion accepted",
         color = "set")
}

sim_run <- function(
  nsim = 2, n = c(50, 100), dgp = dgp_dicp, modFUN = Polr,
  fml = Y ~ X1 + X2 + X3, type = "residual", test = "HSIC",
  pfml = A ~ X1 + X2, wtype = "none", blfix = TRUE, fname = "dgp1",
  pFUN = .policy_clipped, save = TRUE, plot = TRUE, outdir = "inst"
) {

  res1 <- lapply(n, \(tn) {
    res <- lapply(1:nsim, function(x) {
      ret <- try(sim(n = tn, dgp = dgp, fml = fml, modFUN = modFUN,
                     type = type, pfml = pfml, wtype = wtype, blfix = blfix,
                     pFUN = pFUN),
                 silent = FALSE)
      if (inherits(ret, "try-error"))
        return(NULL)
      ret
    })
    cat("Finished n =", tn, "\n")
    tres <- do.call("rbind", res)

    if (save)
      write.csv(
        tres,
        file.path(outdir, paste0(type, "-", fname, "-sim-pvals_n", tn, ".csv")),
        row.names = FALSE, quote = FALSE
      )

    colMeans(tres >= 0.05, na.rm = TRUE)
  })

  out1 <- .summarize_results(res1, n)

  if (plot)
    print(.visualize_results(out1))
  if (save) {
    ggplot2::ggsave(file.path(outdir, paste0(type, "-", fname, "-sim-tram_icp.pdf")),
           height = 5, width = 5)
    write.csv(
      out1,
      file.path(outdir, paste0(type, "-", fname, "-sim-tram_icp.csv")),
      row.names = FALSE, quote = FALSE
    )
  }

  out1

}

.policy_clipped <- function(x1, x2, ...) {
  # p0 <- plogis(rnorm(n))
  p0 <- plogis((x1 + x2) / 2)
  p0 <- ifelse(p0 < 0.1, 0.1, p0)
  p0 <- ifelse(p0 > 0.9, 0.9, p0)
  cbind(p0, 1 - p0)
}

.onesim <- function(.par) {
}

simfun <- function(mod, n, nanc, panc, ndec, pdec, penv, types, ws, blfix,
                   resp, env, preds, fml, niter) {
  params <- list(mod = mod, n = n, nanc = nanc, panc = panc, ndec = ndec,
                 pdec = pdec, penv = penv)
  cat("Run:", paste(names(params), unlist(params), sep = "="), "\n")
  plan(multisession(workers = 28))
  res <- future_lapply(1:nsim, \(iter) {
    # library(magrittr)
    mFUN <- tramicp:::.mod_from_name(mod)
    d <- structure(NA, class = "try-error")
    while (inherits(d, "try-error")) {
      d <- try(dgp_random_dag(n = n, mod = mod, K = ifelse(mod == "polr", 8, 6),
                              dag = random_dag(nanc = nanc, ndec = ndec,
                                               panc = panc, pdec = pdec,
                                               penv = penv)))
    }
    paY <- attr(d, "paY")
    chE <- attr(d, "chE")
    res <- lapply(types, \(ttype) {
      ttest <- .get_test(ttype)
      kbw <- if (mod == "polr") 0.01 else 0
      oicp <- attr(d, "oracle_icp")
      pvals <- simall(d = d, fml = as.formula(fml), modFUN = mFUN, type = ttype,
                      wtype = ws, blfix = blfix, resp = resp,
                      env = env, preds = preds)
      inv <- tramicp:::.get_invariant_set(pvals)
      tibble(paY = paste0(paY, collapse = "+"), chE = paste0(chE, collapse = "+"),
             oicp = paste0(oicp, collapse = "+"), inv = inv,
             npaY = length(attr(d, "paY")), type = ttype, test = ttest,
             pvals = list(stack(pvals)))
    })
    bind_rows(res)
  }, future.seed = TRUE)
  out <- bind_rows(res)
  if (!dir.exists("sim_results"))
    dir.create("sim_results")
  write_csv(out |> select(-pvals), file.path("sim_results", paste0(
    "results_", paste(names(params), unlist(params), sep = "", collapse = "_"),
    ".csv")))
  write_csv(out |> unnest(c(pvals)), file.path("sim_results", paste0(
    "raw_results_", paste(names(params), unlist(params), sep = "", collapse = "_"),
    ".csv")))
  out
}

read_from_single <- function(files) {
  tibble(file = files) |>
    mutate(dat = map(file, ~ read_csv(.x, show_col_types = FALSE))) |>
    separate(file, into = c("jnk1", "jnk2", "mod", "n", "nanc", "panc", "ndec",
                            "pdec", "penv"),
             sep = "_") |>
    unnest(c(dat)) |>
    mutate_at(c("n", "nanc", "panc", "ndec", "pdec", "penv"), parse_number) |>
    mutate(mod = str_remove(mod, "mod"))
}
