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
                type = ttype, test = ttest, verbose = FALSE, greedy = TRUE)
    expect_length(pvalues(res, "predictor"), 3)
    expect_type(res$candidate_causal_predictors, "character")
  })

  ### Weights
  ww <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
             weights = abs(rnorm(nrow(d))), verbose = FALSE, greedy = TRUE)
  expect_length(pvalues(ww, "predictor"), 3)

})

test_that("All aliases work", {
  set.seed(123)
  library("survival")
  ### All aliases
  expect_no_error({
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
                   test = ttest, verbose = FALSE, greedy = TRUE)
        if (dtype == "lm") {
          res <- lmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                       test = ttest, verbose = FALSE, greedy = TRUE)
        } else if (dtype == "cotram") {
          res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                        test = ttest, verbose = FALSE, family = "poisson",
                        greedy = TRUE)
        } else if (dtype == "polr" && tcomb > 1) {
          res <- polrICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                         test = ttest, verbose = FALSE, greedy = TRUE)
        }
      })
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
                     test = ttest, verbose = FALSE, greedy = TRUE)
    res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                  test = ttest, verbose = FALSE, family = "poisson", greedy = TRUE)
  })

  ### binary glm
  set.seed(1334)
  d <- dgp_dicp(mod = "binary")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                  test = ttest, verbose = FALSE, family = "binomial", greedy = TRUE)
    expect_true(res$candidate_causal_predictors %in% c("X2", "Empty"))
  })
})
