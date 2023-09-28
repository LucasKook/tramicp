combs <- c("wald" = "wald", "residual" = "independence", "residual" = "cor.test")

dtypes <- list(
  "boxcox" = BoxCoxICP,
  "weibull" = SurvregICP,
  "colr" = ColrICP,
  "coxph" = CoxphICP,
  "lm" = LmICP,
  "polr" = PolrICP,
  "cotram" = cotramICP,
  "binary" = \(...) glmICP(..., family = "binomial")
)

test_that("main function works with mandatory terms", {
  set.seed(123)
  d <- dgp_dicp(n = 1e3, mod = "polr")

  ### Main function with ordinal outcome and tram::Polr
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- dicp(Y ~ X1 + X3, data = d, env = ~ E, modFUN = Polr,
                type = ttype, test = ttest, mandatory = ~ X2,
                verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
    expect_length(pvalues(res, "predictor"), 2)
    expect_type(res$candidate_causal_predictors, "character")
  })

})

test_that("All aliases work with mandatory terms", {

  set.seed(123)
  library("survival")

  ### All aliases
  lapply(seq_along(dtypes), \(didx) {
    dtype <- names(dtypes[didx])
    FUN <- dtypes[[didx]]
    d <- dgp_dicp(mod = dtype)
    dotest <- seq_along(combs)
    if (dtype == "weibull")
      dotest <- seq_along(combs)[-1:-2]
    lapply(dotest, \(tcomb) {
      ttype <- names(combs)[tcomb]
      ttest <- unname(combs[tcomb])
      res <- FUN(Y ~ X1, data = d, env = ~ E, type = ttype,
                 test = ttest, verbose = FALSE, mandatory = ~ X2)
      expect_length(pvalues(res, "set"), 2)
      if (dtype == "lm") {
        res <- lmICP(Y ~ X1, data = d, env = ~ E, type = ttype,
                      test = ttest, verbose = FALSE, mandatory = ~ X2)
      } else if (dtype == "cotram") {
        res <- glmICP(Y ~ X1, data = d, env = ~ E, type = ttype,
                      test = ttest, verbose = FALSE, family = "poisson",
                      mandatory = ~ X2)
      } else if (dtype == "polr" && tcomb > 1) {
        res <- polrICP(Y ~ X1, data = d, env = ~ E, type = ttype,
                        test = ttest, verbose = FALSE,
                        mandatory = ~ X2)
      }
      expect_length(pvalues(res, "set"), 2^1)
    })
  })

})
