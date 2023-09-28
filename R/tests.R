
.implemented_tests <- function() {
  c("independence", "HSIC", "t.test", "var.test", "combined", "wald",
    "cor.test", "spearman", "gcm.test")
}

# Tests -------------------------------------------------------------------

.gcm_test <- function(r, e, controls) {
  dR <- NCOL(r)
  dE <- NCOL(e)
  nn <- NROW(r)
  if (dR > 1 || dE > 1) {
    R_mat <- matrix(r, nrow = nn, ncol = dE) * e
    sigma <- crossprod(R_mat) / nn - tcrossprod(colMeans(R_mat))
    eig <- eigen(sigma)
    siginvhalf <- eig$vectors %*% diag(eig$values^(-1/2)) %*% t(eig$vectors)
    tstat <- siginvhalf %*% colSums(R_mat) / sqrt(nn)
    p.value <- stats::pchisq(sum(tstat^2), df = dE, lower.tail = FALSE)
  }
  else {
    R <- r * e
    R.sq <- R^2
    meanR <- mean(R)
    tstat <- sqrt(nn) * meanR / sqrt(mean(R.sq) - meanR^2)
    p.value <- 2 * stats::pnorm(abs(tstat), lower.tail = FALSE)
  }
  return(list(p.value = p.value, test.statistic = tstat,
              reject = (p.value < controls$alpha)))
}

.dhsic_test <- function(r, e, controls) {
  dHSIC::dhsic.test(r, e, alpha = controls$alpha, method = controls$method,
                    B = controls$B)
}

.t_test <- function(r, e, controls) {
  stats::t.test(r ~ e)
}

.cor_test <- function(r, e, controls) {
  stats::cor.test(r, e)
}

.var_test <- function(r, e, controls) {
  stats::var.test(r ~ e)
}

.combined_test <- function(r, e, controls) {
  shift <- stats::t.test(r ~ e)
  scale <- stats::var.test(r ~ e)
  list("p.value" = 2 * min(shift$p.value, scale$p.value),
       shift = shift, scale = scale)
}

.spearman_test <- function(r, e, controls) {
  tst <- coin::spearman_test(r ~ e)
  list(p.value = coin::pvalue(tst), test = tst)
}

.indep_test <- function(r, e, controls) {
  tst <- coin::independence_test(
    r ~ e, teststat = controls$teststat, distribution = controls$distribution,
    xtrafo = controls$xtrafo, ytrafo = controls$ytrafo
  )
  list(p.value = coin::pvalue(tst), test = tst)
}
