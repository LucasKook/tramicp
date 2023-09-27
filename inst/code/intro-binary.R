## Intro invariance example with binary outcome
## LK 2023

set.seed(-42)

# DEPs --------------------------------------------------------------------

devtools::load_all()
devtools::load_all("inst/helpers")
library("ranger")
library("pROC")
library("tidyverse")
theme_set(theme_bw())

# DGP ---------------------------------------------------------------------

dgp <- function(n = 5e2) {
  E <- sample(0:1, n, TRUE)
  X1 <- -E + rnorm(n)
  Y <- as.numeric(0.5 * X1 > rlogis(n))
  X2 <- Y + 0.8 * E + rnorm(n)
  data.frame(Y = Y, X1 = X1, X2 = X2, E = E)
}

d <- dgp()

icp <- glmICP(Y ~ X1 + X2, data = d, env = ~ E, family = "binomial",
              test = "gcm.test", verbose = FALSE)
pvalues(icp, "set")

m0 <- glm(Y ~ 1, data = d, family = "binomial")
r0 <- residuals.binglm(m0)
re0 <- d$E - mean(d$E)

m1 <- glm(Y ~ X1, data = d, family = "binomial")
r1 <- residuals.binglm(m1)
re1 <- d$E - predict(ranger(factor(E) ~ X1, data = d, probability = TRUE),
                     data = d)$predictions[, 2]

m2 <- glm(Y ~ X2, data = d, family = "binomial")
r2 <- residuals.binglm(m2)
re2 <- d$E - predict(ranger(factor(E) ~ X2, data = d, probability = TRUE),
                     data = d)$predictions[, 2]

m12 <- glm(Y ~ X1 + X2, data = d, family = "binomial")
r12 <- residuals.binglm(m12)
re12 <- d$E - predict(ranger(factor(E) ~ X1 + X2, data = d, probability = TRUE),
                      data = d)$predictions[, 2]

pd <- data.frame(
  set = factor(rep(c("1", "X1", "X2", "X1+X2"), each = nrow(d)),
               levels = c("1", "X1", "X2", "X1+X2")),
  R = c(r0, r1, r2, r12),
  rE = c(re0, re1, re2, re12)
)

lbl <- pvalues(icp, "set")
ldf <- structure(paste0(names(lbl), " (p-value ", fmt(lbl), ")"), names = names(lbl))

ggplot(pd %>% dplyr::filter(set != "1"), aes(x = rE, y = R)) +
  geom_point(alpha = 0.5, size = rel(1.2)) +
  facet_grid(~ set, labeller = as_labeller(ldf)) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(y = "Score residuals", x = "Residualized environment")

ggsave("inst/figures/binary-invariance.pdf", height = 2.7, width = 7)
