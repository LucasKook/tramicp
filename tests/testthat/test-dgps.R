
test_that("Coefs from DGP can be recovered", {
  mods <- c("weibull", "lm", "coxph", "colr", "polr")
  tcfx <- c(1, 0)

  lapply(mods, \(tmod) {
    set.seed(1)
    tmodFUN <- tramicp:::.mod_from_name(tmod)
    dgp <- function(n = 1e4, mod = tmod, ...)
      dgp_dicp(n = n, mod = mod, cfx = tcfx, ... = ...)
    d <- dgp()
    expect_true(max(abs(coef(tmodFUN(Y ~ X1 + X2, data = d))) - abs(tcfx)) < 0.05)
  })
})
