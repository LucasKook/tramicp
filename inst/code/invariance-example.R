# Invariance with score residuals
# LK 2023

set.seed(24101968)

# DEPs --------------------------------------------------------------------

devtools::load_all()
devtools::load_all("inst/helpers")
library("ranger")
library("cotram")
library("tidyverse")

# GEN ---------------------------------------------------------------------

dag <- random_dag(nanc = 2, ndec = 2, panc = 0.8, penv = 0.8, pdec = 0.5)
pdf("inst/figures/dag.pdf", width = 3.5, height = 3.5)
Rgraphviz::plot(as(abs(dag$dag), "graphNEL"))
dev.off()

df <- dgp_random_dag(n = 4e2, mod = "cotram", dag = dag)
mf <- cotram
# plot(df)

# RUN ---------------------------------------------------------------------

### Invariant model
fm_inv <- reformulate(dag$paY, "Y")
m_inv <- mf(fm_inv, data = df)
rf_inv <- ranger(reformulate(dag$paY, "E"), data = df, probability = TRUE)
df$e_inv <- as.numeric(df$E) - 1 - predict(rf_inv, data = df)$predictions[, 2]

### Non-invariant model
fm_non <- reformulate(ss <- "X1", "Y")
m_non <- mf(fm_non, data = df)
rf_non <- ranger(reformulate(ss, "E"), data = df, probability = TRUE)
df$e_non <- as.numeric(df$E) - 1 - predict(rf_non, data = df)$predictions[, 2]

### Run full TRAMICP
icp <- cotramICP(Y ~ X1 + X2 + X3 + X4, data = df, env = ~ E, type = "residual",
                 test = "gcm.test", verbose = FALSE)
icp
pvals <- pvalues(icp, "set")

pi <- pvals[paste0(dag$paY, collapse = "+")]
pni <- pvals[paste0(ss, collapse = "+")]

# VIS ---------------------------------------------------------------------

pdat <- data.frame(
  E = c(df$e_inv, df$e_non),
  residuals = c(residuals(m_inv), residuals(m_non)),
  mod = rep(c("inv", "non"), each = nrow(df))
)

ggplot(pdat, aes(x = E, y = residuals)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ mod, labeller = as_labeller(c(
    "inv" = paste0("Invariant set (p-value ", fmt(pi), ")"),
    "non" = paste0("Non-invariant set (p-value ", fmt(pni), ")")))) +
  theme_bw() +
  labs(x = "Residualized environment", y = "Score residuals")

ggsave("inst/figures/invariance-example.pdf", height = 2.5, width = 5, scale = 1.2)
