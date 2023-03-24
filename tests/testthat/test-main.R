combs <- c("confint" = "independence", "wald" = "wald",
           "residual" = "independence", "residual" = "HSIC")

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

  ### BoxCox
  d <- dgp_dicp(mod = "boxcox")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- BoxCoxICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                     test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### Survreg
  d <- dgp_dicp(mod = "weibull")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- SurvregICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                      test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### survival::survreg
  library("survival")
  d$surv <- Surv(d$Y)
  lapply(seq_along(combs)[-1:-2], \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- dicp(surv ~ X1 + X2, data = d, env = ~ E, type = ttype,
                test = ttest, verbose = FALSE, modFUN = survreg)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### Colr
  d <- dgp_dicp(mod = "colr")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- ColrICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                   test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### Coxph
  d <- dgp_dicp(mod = "coxph")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- CoxphICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                    test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### Lm and lm
  d <- dgp_dicp(mod = "lm")
  lapply(seq_along(combs), \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- LmICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                 test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
    res <- dicp(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                test = ttest, verbose = FALSE, modFUN = lm)
    expect_length(pvalues(res, "set"), 2^2)
  })

  ### MASS::polr
  d <- dgp_dicp(mod = "polr")
  lapply(seq_along(combs)[-1], \(tcomb) {
    ttype <- names(combs)[tcomb]
    ttest <- unname(combs[tcomb])
    res <- mpolrICP(Y ~ X1 + X2, data = d, env = ~ E, type = ttype,
                    test = ttest, verbose = FALSE)
    expect_length(pvalues(res, "set"), 2^2)
  })

  expect_error(mpolrICP(Y ~ X1 + X2, data = d, env = ~ E, type = "confint",
                        verbose = FALSE))

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
