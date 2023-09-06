# Show non (pZ, h)-identifiability
# LK, SS, Nov22

set.seed(0)

# Dependencies ------------------------------------------------------------

library("MASS")
devtools::load_all()
library("tidyverse")
library("coin")
library("generalhoslem")

oup <- file.path("inst", "results", Sys.Date())
if (!dir.exists(oup))
  dir.create(oup)

# FUNs --------------------------------------------------------------------

# CDF
# theta1 / (1 + theta1)         | theta1 / (1 + theta2)                 | theta1 * beta3 / (1 + theta1 * beta3) |
# theta2 / (1 + theta1 * beta2) | theta2 * beta2 / (1 + theta2 * beta2) | theta2 * beta3 / (1 + theta2 * beta3) |
# 1                             | 1                                     | 1                                     |

# PDF
# theta1 / (1 + theta1)         | theta1 / (1 + theta2)                 | theta1 * beta3 / (1 + theta1 * beta3) |
# R2 - R1
# R3 - R2 - R1

mk_distr <- function(gamma, xi, type = c("pdf", "cdf")) {
  type <- match.arg(type)
  ret <- rbind(outer(seq_along(gamma), seq_along(xi), \(x, y) {
    gamma[x] * xi[y] / (1 + gamma[x] * xi[y])
  }), 1)
  if (type == "pdf")
    ret <- apply(ret, 2, \(p) diff(c(0, p)))
  structure(ret, gamma = gamma, xi = xi)
}

copula <- function(x, y, alpha = 1.1) {
  exp(-((-log(x))^alpha + (-log(y))^alpha)^(1/alpha))
}

jcdf2pdf <- function(cdf) {
  pdf <- cdf
  pdf[] <- 0
  pdf[1,] <- diff(c(0, cdf[1,]))
  pdf[,1] <- diff(c(0, cdf[,1]))
  for (x in 2:nrow(cdf)) {
    for (y in 2:ncol(cdf)) {
      pdf[x, y] <- cdf[x, y] - cdf[x - 1, y] - cdf[x, y - 1] + cdf[x - 1, y - 1]
    }
  }
  pdf
}

ct2df <- function(x) {
  nm <- names(x)
  x_tab <- as.table(as.matrix(x))
  rownames(x_tab) <- nm
  x_df <- as.data.frame(x_tab)
  ret <- data.frame(lapply(x_df[seq_len(ncol(x_df) - 1)], function(col) {
    rep(col, x_df$Freq)
  }))
  colnames(ret) <- c("X", "Y")
  ret$Y <- ordered(ret$Y)
  ret
}

### Find solution for fixed marginal of Y
find_sol <- function(dM, tP = 3, try = 1e3) {
  K <- dim(dM)[1]
  sols <- list()
  for (iter in seq_len(1e3)) {
    gamma <- round(seq(runif(1, 0.1, 0.8), runif(
      1, 1.5, 2.5), length.out = K - 1), 1)
    xi <- 1 / c(1, round(runif(tP - 1, min = 1.1, max = 5), 1))
    dC <- mk_distr(gamma, xi)
    sol <- try(solve(dC) %*% dM)
    if (inherits(sol, "try-error"))
      next
    isSOL <- all(sol > 0) & all(sol < 1)
    if (isSOL)
      sols <- c(sols, list(list(gamma = attr(dC, "gamma"), xi = attr(dC, "xi"))))
  }
  sols
}

gen_dat <- function(
    n = 1e3, tK = 3, tP = 3, alphas = c(1.1, 1.2, 1.3), dX1 = rep(1/tP, tP),
    g1 = c(0.5, 1.5), b1 = c(1, 1.4, 1.8), g2 = c(0.5, 1.6), b2 = c(1, 4.6, 2.4)
) {
  # 1. Generate density Y | X1
  dYX1 <- mk_distr(g1, b1)

  # 2. Compute marginal density of Y
  dY <- dYX1 %*% dX1

  # 3. Generate density Y | X2
  dYX2 <- mk_distr(g2, b2)
  # print(dYX2)
  # print(Matrix::rankMatrix(dYX2))

  # 4. Compute marginal density of X2
  dX2 <- solve(dYX2) %*% dY

  # 5. Joint density of Y, X1
  jdYX1 <- dYX1 * matrix(dX1, nrow = tK, ncol = tP, byrow = TRUE)

  # 6. Conditional density of X1 | Y
  dX1Y <- t(jdYX1) / matrix(dY, nrow = tP, ncol = tK, byrow = TRUE)

  # 7. Joint density of Y, X2
  jdYX2 <- dYX2 * matrix(dX2, nrow = tK, ncol = tK, byrow = TRUE)

  # 8. Conditional density of X2 | Y
  dX2Y <- t(jdYX2) / matrix(dY, nrow = tK, ncol = tK, byrow = TRUE)

  # 9. Conditional copula
  pX1Y <- apply(dX1Y, 2, cumsum) # CDF X1 | Y
  pX2Y <- apply(dX2Y, 2, cumsum) # CDF X2 | Y

  pCop <- lapply(1:tK, \(y) outer(
    pX1Y[, y], pX2Y[, y], copula, alpha = alphas[y]))

  dCop <- lapply(pCop, jcdf2pdf)

  # Joint -------------------------------------------------------------------

  dJ <- array(unlist(dCop), c(tP, tK, tK)) *
    array(rep(dY, each = tP * tK), dims <- c(tP, tK, tK))
  dimnames(dJ) <- dnms <- list(paste0("X1", 1:tP), paste0("X2", 1:tK), paste0(
    "Y", 1:tK))

  # Sample ------------------------------------------------------------------

  smpl <- array(rmultinom(n = 1, size = n, prob = c(dJ)), dim = dims,
                dimnames = dnms)
  dfj <- do.call("rbind", lapply(seq_len(tK), \(x) cbind(ct2df(smpl[, , x]), x)))
  colnames(dfj) <- c("X1", "X2", "Y")
  dfj$X1 <- factor(dfj$X1, ordered = FALSE, labels = c("A", "B", "C"))
  dfj$X2 <- factor(dfj$X2, ordered = FALSE, labels = c("A", "B", "C"))
  dfj$Y <- factor(dfj$Y, ordered = TRUE)

  # Return
  dfj

}

# Simulation --------------------------------------------------------------

nsim <- 1e3
tn <- 1e4
tg1 <- c(0.5, 1.5)
tb1 <- c(1, 1.4, 1.8)

# find_sol(mk_distr(tg1, tb1) %*% rep(1/3, 3), try = 3e3)
tg2 <- tg1 # c(0.8, 2.4)
tb2 <- tb1 # c(1, 0.76, 0.22)

talp <- rep(1, 3) # c(1.1, 1.2, 1.3)

t1 <- c(log(tg1), log(1 / tb1[-1]))
t2 <- c(log(tg2), log(1 / tb2[-1]))

pb <- txtProgressBar(min = 0, max = nsim, style = 3)
res <- lapply(seq_len(nsim), \(iter) {
  setTxtProgressBar(pb, iter)
  d <- gen_dat(n = tn, g1 = tg1, b1 = tb1, g2 = tg2, b2 = tb2, alphas = talp)
  m1 <- polr(Y ~ X1, data = d)
  m2 <- polr(Y ~ X2, data = d)
  d$r1 <- tramicp:::residuals.polr(m1)
  d$r2 <- tramicp:::residuals.polr(m2)
  list(
    bias_x1 = (c(m1$zeta, coef(m1)) - t1) / t1,
    bias_x2 = (c(m2$zeta, coef(m2)) - t2) / t2,
    pval_x1 = pvalue(independence_test(r1 ~ X1, data = d, ytrafo = rank)),
    pval_x2 = pvalue(independence_test(r2 ~ X2, data = d, ytrafo = rank))
  )
}) %>% transpose()

bias_x1 <- bind_rows(res$bias_x1)
# write_csv(bias_x1, file.path(oup, "bias_x1.csv"))
boxplot(bias_x1, ylab = "bias", xlab = "coef")

bias_x2 <- bind_rows(res$bias_x2)
# write_csv(bias_x2, file.path(oup, "bias_x2.csv"))
boxplot(bias_x2, ylab = "bias", xlab = "coef")

(rf1 <- mean((pvals1 <- unlist(res$pval_x1)) < 0.05))
(rf2 <- mean((pvals2 <- unlist(res$pval_x2)) < 0.05))

# write_csv(data.frame(p = pvals1), file.path(oup, "pvals_x1.csv"))
# write_csv(data.frame(p = pvals2), file.path(oup, "pvals_x2.csv"))

# Visualization -----------------------------------------------------------

lx1 <- bias_x1 %>%
  dplyr::rename(
    "vartheta[11]" = `1|2`,
    "vartheta[12]" = `2|3`,
    "beta[12]" = X1B,
    "beta[13]" = X1C,
  ) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "estimate")

lx2 <- bias_x2 %>%
  dplyr::rename(
    "vartheta[21]" = `1|2`,
    "vartheta[22]" = `2|3`,
    "beta[22]" = X2B,
    "beta[23]" = X2C,
  ) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "estimate")

pd <- full_join(lx1, lx2) %>%
  mutate(parameter = factor(parameter, levels = c("vartheta[11]", "vartheta[12]", "beta[12]", "beta[13]",
                                                 "vartheta[21]", "vartheta[22]", "beta[22]", "beta[23]")))

ggplot(pd, aes(x = parameter, y = estimate)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "relative bias", x = element_blank()) +
  scale_x_discrete(labels = parse(text = levels(pd$parameter))) +
  geom_hline(yintercept = 0, linetype = 3, color = 2)

ggsave(file.path("inst", "figures", "bias.pdf"), height = 3.5, width = 4.5)

# GOF ---------------------------------------------------------------------

nsim <- 1e2
tn <- 1e4
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
pvals <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)
  d <- gen_dat(n = tn, g1 = tg1, b1 = tb1, g2 = tg2, b2 = tb2, alphas = talp)
  m1 <- polr(Y ~ X1, data = d)
  m2 <- polr(Y ~ X2, data = d)
  c(m1 = pulkrob.chisq(m1, "X1")$p.value, m2 = pulkrob.chisq(m2, "X2")$p.value)
}) %>% bind_rows()

pvals %>% summarize_all(~ mean(.x < 0.05))
