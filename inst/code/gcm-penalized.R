### Test rates by perturbation
### LK 2023

set.seed(1234)

# DEPs --------------------------------------------------------------------

devtools::load_all()
library("tidyverse")
library("tramnet")
library("glmnet")

### Make out-directory
odir <- file.path("inst", "results", Sys.Date())
if (!dir.exists(odir))
  dir.create(odir, recursive = TRUE)

# FUNs --------------------------------------------------------------------

mME <- BoxCox
fam <- mlt:::.Normal()
pZ <- fam$p
qZ <- fam$q
rZ <- \(n) qZ(runif(n))
dd2d <- fam$dd2d

p0 <- \(y) pchisq(y, df = 3)
q0 <- \(p) qchisq(p, df = 3)

f1 <- \(x) 0.5 * x
g0 <- \(z) q0(pZ(z))

tp0 <- 10

g0 <- identity

dgp <- function(n = 1e3, p0 = tp0) {
  E <- sample(0:1, n, TRUE)
  X1 <- E + rnorm(n)
  X2 <- rnorm(n)
  Y <- g0(rZ(n) + f1(X1) + X2)
  X3 <- Y + rnorm(n) - E
  XZ <- matrix(rnorm(n * p0), nrow = n, ncol = p0)
  data.frame(Y = Y, X1 = X1, X2 = X2, X3 = X3, E = E, XZ = XZ)
}

mglmnet <- function(x, y, alpha, lambda) {
  m <- glmnet(x, y, alpha = alpha, lambda = lambda)
  m$x <- x
  m$y <- y
  class(m) <- c("mglmnet", class(m))
  m
}

residuals.mglmnet <- function(m) {
  m$y - predict(m, newx = m$x)
}

pvals <- lapply(1:3e2, \(iter) {
  if (iter %% 20 == 0)
    cat("--- done", iter, "iterations ---\n")
  df <- dgp(5e1)
  rf <- ranger(reformulate(c("X1", "X2", paste0("XZ.", 1:tp0)), "factor(E)"),
               data = df, probability = TRUE)
  df$rE <- df$E - predict(rf, data = df)$predictions[, 2]
  m0 <- mME(Y ~ 1, data = df, prob = c(0.001, 0.999))
  X <- scale(as.matrix(df[, c(2:3, 6:ncol(df))]))
  m <- mglmnet(X, df$Y, alpha = 0, lambda = 10)
  m2 <- lm(Y ~ E * (X1 + X2 + XZ.1 + XZ.2 + XZ.3 + XZ.4 + XZ.5 + XZ.6 + XZ.7 +
                      XZ.8 + XZ.9 + XZ.10), data = df)
  wald <- summary(glht(m2, c("E = 0", "E:X1 = 0", "E:X2 = 0", paste0(
    "E:XZ.", 1:tp0, " = 0")), vcov = vcov), test = Chisqtest())
  tribble(~ "test",                                                ~ "pvalue",
          "COR",   cor.test(residuals(m), df$E)$p.value,
          "GCM",   .gcm_test(residuals(m), df$rE, list(alpha = 0.05))$p.value,
          "WALD", c(wald$test$pvalue)
  )
}) %>% bind_rows()

ggplot(pvals, aes(x = pvalue, linetype = test, group = test)) +
  stat_ecdf() +
  geom_abline(intercept = 0, slope = 1) +
  theme_bw() +
  labs(x = "p-value", y = "ECDF", linetype = "Invariance test")

ggsave(file.path("inst", "figures", "gcm-penalized-breaks-correlation-rest.pdf"),
       height = 3, width = 5)

### P-value reported in the paper
pvals %>%
  nest_by(test) %>%
  summarise(pval = binom.test(sum(data < 0.05), nrow(data), p = 0.05,
                              alternative = "greater")$p.value)
