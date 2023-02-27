# Test DGPs

mods <- c("weibull", "lm", "coxph", "colr", "polr", "binary")
tcfx <- c(1, 0)

lapply(mods, \(tmod) {
  tmodFUN <- tramicp:::.mod_from_name(tmod)
  dgp <- function(n = 1e4, mod = tmod, ...)
    dgp_dicp(n = n, mod = mod, cfx = tcfx, ... = ...)
  d <- dgp()
  confint(tmodFUN(Y ~ X1 + X2, data = d))
})
