
.implemented_tests <- function() {
  c("independence", "HSIC", "t.test", "var.test", "combined", "wald",
    "cor.test", "spearman")
}

# Tests -------------------------------------------------------------------

.dhsic_test <- function(r, e, controls) {
  dhsic.test(r, e, alpha = controls$alpha, method = controls$method,
             B = controls$B)
}

.t_test <- function(r, e, controls) {
  t.test(r ~ e)
}

#' @importFrom stats cor.test
.cor_test <- function(r, e, controls) {
  cor.test(r, e)
}

.var_test <- function(r, e, controls) {
  var.test(r ~ e)
}

.combined_test <- function(r, e, controls) {
  shift <- t.test(r ~ e)
  scale <- var.test(r ~ e)
  list("p.value" = 2 * min(shift$p.value, scale$p.value),
       shift = shift, scale = scale)
}

#' @importFrom coin spearman_test
.spearman_test <- function(r, e, controls) {
  tst <- spearman_test(r ~ e)
  list(p.value = pvalue(tst), test = tst)
}

#' @importFrom coin independence_test pvalue
.indep_test <- function(r, e, controls) {
  tst <- independence_test(
    r ~ e, teststat = controls$teststat, distribution = controls$distribution,
    xtrafo = controls$xtrafo, ytrafo = controls$ytrafo
  )
  list(p.value = pvalue(tst), test = tst)
}
