
.implemented_tests <- function() {
  c("independence", "HSIC", "t.test", "var.test", "combined", "wald",
    "cor.test", "spearman", "gcm.test")
}

# Tests -------------------------------------------------------------------

.gcm_test <- function(r, e, controls) {
  nn <- NA
  nsim <- controls$B
  if (NCOL(r) > 1 || NCOL(e) > 1) {
    d_X <- NCOL(r)
    d_Y <- NCOL(e)
    nn <- NROW(r)
    R_mat <- rep(r, times = d_Y) * as.numeric(as.matrix(
      e)[, rep(seq_len(d_Y), each = d_X)])
    dim(R_mat) <- c(nn, d_X * d_Y)
    R_mat <- t(R_mat)
    R_mat <- R_mat/sqrt((rowMeans(R_mat^2) - rowMeans(R_mat)^2))
    test.statistic <- max(abs(rowMeans(R_mat))) * sqrt(nn)
    test.statistic.sim <- apply(abs(R_mat %*% matrix(
      rnorm(nn * nsim), nn, nsim)), 2, max) / sqrt(nn)
    p.value <- (sum(test.statistic.sim >= test.statistic) + 1)/(nsim + 1)
  }
  else {
    nn <- ifelse(is.null(dim(r)), length(r), dim(r)[1])
    R <- r * e
    R.sq <- R^2
    meanR <- mean(R)
    test.statistic <- sqrt(nn) * meanR/sqrt(mean(R.sq) - meanR^2)
    p.value <- 2 * pnorm(abs(test.statistic), lower.tail = FALSE)
  }
  return(list(p.value = p.value, test.statistic = test.statistic,
              reject = (p.value < controls$alpha)))
}

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
