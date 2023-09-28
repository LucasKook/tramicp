combs <- c("wald" = "wald", "residual" = "independence", "residual" = "HSIC")

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

test_that("main function works", {
  set.seed(123)
  d <- dgp_dicp(n = 1e3, mod = "polr")

  ### Main function with ordinal outcome and tram::Polr
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr,
                type = ttype, test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^3)
    expect_length(pvalues(res, "predictor"), 3)
    expect_type(res$candidate_causal_predictors, "character")
  })

  ### Weights
  ww <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
             weights = abs(rnorm(nrow(d))), verbose = FALSE)
  expect_length(pvalues(ww, "set"), 2^3)
  expect_length(pvalues(ww, "predictor"), 3)

})

test_that("All aliases work", {

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
      res <- FUN(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                 test = ttest, verbose = FALSE)
      expect_length(pvalues(res, "set"), 2^2)
      if (dtype == "lm") {
        res <- lmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                    test = ttest, verbose = FALSE)
      } else if (dtype == "cotram") {
        res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                      test = ttest, verbose = FALSE, family = "poisson")
      } else if (dtype == "polr" && tcomb > 1) {
        res <- polrICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                        test = ttest, verbose = FALSE)
      }
      expect_length(pvalues(res, "set"), 2^2)
    })
  })
})

test_that("Output of cotramICP and glmICP", {

  ### cotram and poisson glm
  set.seed(13312)
  d <- dgp_dicp(mod = "cotram")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- cotramICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                     test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
    res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                  test = ttest, verbose = FALSE, family = "poisson")
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### binary glm
  set.seed(1334)
  d <- dgp_dicp(mod = "binary")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                  test = ttest, verbose = FALSE, family = "binomial")
    expect_length(pvalues(res, "set"), 2^2)
    expect_true(res$candidate_causal_predictors %in% c("X2", "Empty"))
  })
})

test_that("Multi-environment GCM works", {
  set.seed(1234)
  d <- dgp_dicp(mod = "binary")
  expect_no_error(glmICP(Y ~ X1 + X2, data = d, env = ~ E + X3, family = "binomial"))
})

test_that("argument checks work", {
  d <- dgp_dicp(mod = "boxcox")
  expect_error(BoxCoxICP("Y ~ X1", d, ~ E))
  expect_error(BoxCoxICP(Y ~ X1, d, "E"))
  expect_error(dicp(Y ~ X1, d, ~ E, modFUN = NA))
  expect_error(dicp(Y ~ X1, d, ~ E + X2, modFUN = "BoxCox", test = "cor.test"))
  expect_error(dicp(Y ~ X1, d, ~ E + X2, modFUN = "BoxCox", test = "t.test"))
})
