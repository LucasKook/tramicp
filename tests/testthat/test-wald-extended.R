combs <- c("wald" = "wald")

dtypes <- list(
  "boxcox" = BoxCoxICP,
  "weibull" = SurvregICP,
  "colr" = ColrICP,
  "coxph" = CoxphICP,
  "lm" = LmICP,
  "polr" = PolrICP,
  "cotram" = cotramICP
)

test_that("main function works", {
  set.seed(123)
  d <- dgp_dicp(n = 1e3, mod = "polr")

  ### Main function with ordinal outcome and tram::Polr
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr,
                type = ttype, test = ttest, verbose = FALSE,
                baseline_fixed = FALSE)
    expect_length(pvalues(res, "set"), 2^3)
    expect_length(pvalues(res, "predictor"), 3)
    expect_type(res$candidate_causal_predictors, "character")
  })

  ### Weights
  ww <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
             weights = abs(rnorm(nrow(d))), verbose = FALSE, baseline_fixed = FALSE)
  expect_length(pvalues(ww, "set"), 2^3)
  expect_length(pvalues(ww, "predictor"), 3)

  ### Mandatory
  md <- dicp(Y ~ X1 + X3, data = d, env = ~ E, mandatory = ~ X2, modFUN = Polr,
             type = "wald", verbose = FALSE, baseline_fixed = FALSE)
  expect_length(pvalues(md, "set"), 2^2)
  expect_length(pvalues(md, "predictor"), 2)

})

test_that("All aliases work", {

  set.seed(123)
  library("survival")

  ### All aliases
  lapply(seq_along(dtypes), \(didx) {
    dtype <- names(dtypes[didx])
    FUN <- dtypes[[didx]]
    d <- dgp_dicp(mod = dtype, cfb = c(-22, 8))
    dotest <- seq_along(combs)
    if (dtype == "weibull")
      dotest <- seq_along(combs)[-1:-2]
    lapply(dotest, \(tcomb) {
      ttype <- names(combs)[tcomb]
      ttest <- unname(combs[tcomb])
      res <- FUN(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                 test = ttest, verbose = FALSE, baseline_fixed = FALSE)
      print(res)
      expect_length(pvalues(res, "set"), 2^2)
    })
  })

})
