### Nonparametric TRAM-GCM with survival forests
### LK 2024

set.seed(2410)
save <- TRUE

# Deps --------------------------------------------------------------------

devtools::load_all()
devtools::load_all("inst/helpers")
library("ranger")
library("tidyverse")
library("survival")
library("tram")

# FUNs --------------------------------------------------------------------

dgp <- function(n = 1e2) {
  X1 <- 0.1 * rlogis(n)
  E <- 2 * X1 + 0.1 * rlogis(n)
  X2 <- 2 * E + 0.1 * rlogis(n)
  rY <- 1.5 + 3 * X1 + 3 * X2 * (X2 < 0.2) +
    2 * (X2 > 0.2) * (X1 > 0) + 0.1 * rlogis(n)
  Y <- (plogis(rY) + 1)^3
  C <- 4 # (X2 < 0.2) * 4 + (X2 > 0.2) * 3.5
  Yobs <- (cens <- as.numeric(Y < C)) * Y + as.numeric(C < Y) * C
  X3 <- -6 * E + 2 * Y + 0.1 * rlogis(n)
  data.frame(Y = Yobs, X1 = X1, X2 = X2, X3 = X3, E = E, cens = cens,
             surv = Surv(Yobs, cens), Ync = Y)
}

plot(dgp()[, -7])

jaccard <- function(x, y) length(intersect(x, y)) / length(union(x, y))

.z_test <- function(x) {
  tst <- t.test(x)
  tibble(mean = tst$estimate, lwr = tst$conf.int[1], upr = tst$conf.int[2])
}

# RUN ---------------------------------------------------------------------

nsim <- 300
tn <- 8e1
preds <- paste0("X", 1:3)

pb <- txtProgressBar(0, nsim, style = 3)
out <- lapply(1:nsim, \(iter) {
  setTxtProgressBar(pb, iter)
  d <- dgp(n = tn)
  suppressMessages({
    icp <- list(
      survforest = survforestICP(reformulate(preds, "surv"), data = d, env = ~ E, verbose = FALSE),
      coxph = coxphICP(reformulate(preds, "surv"), data = d, env = ~ E, verbose = FALSE),
      ranger = rangerICP(reformulate(preds, "Y"), data = d, env = ~ E, verbose = FALSE)
    )
  })
  dplyr::bind_rows(lapply(icp, pvalues, which = "set")) |>
    dplyr::mutate(mod = names(icp), output = unlist(
      lapply(icp, \(x) tramicp:::.pplus(x$cand))))
}) |> dplyr::bind_rows()

# VIZ ---------------------------------------------------------------------

sdat <- out |> group_by(mod) |>
  mutate(lop = str_split(output, "\\+"),
         level = unlist(map(lop, ~ as.numeric("X3" %in% .x))),
         jac = unlist(map(lop, ~ jaccard(.x, c("X1", "X2"))))) |>
  summarize(lev = .binom_test(sum(level), length(level)), pow = .z_test(jac))

power <- ggplot(sdat, aes(x = "", y = pow$mean, color = mod, ymin = pow$lwr,
                          ymax = pow$upr, fill = mod)) +
  geom_pointrange(position = position_dodge(0.2), size = rel(0.3)) +
  theme_bw() +
  theme(text = element_text(size = 13.5), axis.ticks.length.x = ggplot2::unit(0, "cm")) +
  labs(x = element_blank(), y = "Mean Jaccard similarity", color = "model", fill = "model") +
  scale_color_brewer(palette = "Dark2")

level <- ggplot(sdat, aes(x = "", y = lev$frac, color = mod, ymin = lev$lwr,
                          ymax = lev$upr, fill = mod)) +
  geom_hline(yintercept = 0.05, color = "darkred", linetype = 3) +
  geom_pointrange(position = position_dodge(0.2), size = rel(0.3)) +
  # lims(y = c(0, 0.1)) +
  theme_bw() +
  theme(text = element_text(size = 13.5), axis.ticks.length.x = ggplot2::unit(0, "cm")) +
  labs(x = element_blank(), y = "Fraction of incorrect output", color = "Model", fill = "Model") +
  scale_color_brewer(palette = "Dark2")

ggpubr::ggarrange(level, power, common.legend = TRUE)

# Save --------------------------------------------------------------------

if (save) {
  write_csv(out, "inst/results/nonparametric-raw.csv")
  ggsave("inst/figures/nonparametric.pdf", height = 3, width = 7)
}
