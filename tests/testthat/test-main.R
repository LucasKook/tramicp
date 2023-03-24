
test_that("main function works", {
  set.seed(123)
  d <- dgp_dicp(n = 1e3, mod = "polr")

  lapply(c("confint", "wald", "residual"), \(ttype) {
    res <- dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = ttype)
    expect_length(res$pvals, 2^3)
    expect_length(res$ipvals, 3)
    expect_type(res$candidate_causal_predictors, "character")
  })

  dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald")
  dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "wald",
       weights = abs(rnorm(nrow(d))))
  dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual")
  dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = Polr, type = "residual",
       test = "HSIC")
})
