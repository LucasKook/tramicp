### GCM example with score residuals
### LK 2023

set.seed(123)

# Deps --------------------------------------------------------------------

library("tramME")
devtools::load_all()
devtools::load_all("inst/helpers")
library("tidyverse")
library("patchwork")
library("multcomp")

# Smooth ------------------------------------------------------------------

dgp <- function(n = 1e3, g0 = \(z) qchisq(pnorm(z), df = 3)) {
  E <- sample(0:1, n, TRUE)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  f1 <- \(x) 0.5 * x + sin(pi * x) * exp(-x^2/2)
  Y <- g0(rnorm(n) + f1(X1) - X2 / 2)
  X3 <- Y + rnorm(n) - E
  data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E)
}

### Behaviour of GCM test for LmME
res <- lapply(seq_len(3e2), \(iter) {
  if (iter %% 20 == 0)
    cat("Iteration", iter, "done\n")
  df <- dgp(3e2, identity)
  m <- mgcv::gam(Y ~ s(X1, k = 30) + s(X2), data = df)
  m2 <- lm(Y ~ E * (X1 + X2), data = df)
  wald <- summary(glht(m2, c("E = 0", "E:X1 = 0", "E:X2 = 0"), vcov = vcov),
                  test = Chisqtest())
  c(COR = cor.test(residuals(m), df$E)$p.value, GCM = .gcm_test(
    residuals(m), df$E - predict(ranger(
      factor(E) ~ X1 + X2, data = df, probability = TRUE),
      data = df)$predictions[, 2], controls = list(alpha = 0.05))$p.value,
    WALD = c(wald$test$pvalue))
}) %>% bind_rows()

### Vis
p0 <- res %>%
  pivot_longer(names_to = "test", values_to = "pval", everything()) %>%
  ggplot(aes(x = pval, linetype = test)) +
  stat_ecdf() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "p-value", linetype = "Invariance test", y = "ECDF")
p0

### Vis DGP
df <- dgp(3e2, identity)
p1 <- ggplot(df, aes(x = X1, y = Y)) +
  geom_point(alpha = 0.3, size = rel(1)) +
  stat_function(fun = ~ 0.5 * .x + sin(pi * .x) * exp(-.x^2/2), color = 4, linewidth = 1.2) +
  theme_bw()

(p1 + labs(tag = "A")) + (p0 + labs(tag = "B"))
ggsave("inst/figures/gcm-smooth-app.pdf", height = 3.5, width = 8)

### P-value reported in the paper
binom.test(sum(res$WALD < 0.05), nrow(res), p = 0.05)
