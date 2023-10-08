# Case study: SUPPORT2 data
# LK 2023

# Dependencies ------------------------------------------------------------

library("tramicp")
library("tidyverse")
library("forcats")
library("survival")
library("tram")

bpath <- file.path("inst", "results", "case-study")

# Options -----------------------------------------------------------------

### Don't fit/test interactions for TRAMICP-Wald
options("wald_test_interactions" = FALSE)

# Load data ---------------------------------------------------------------

load(url("https://hbiostat.org/data/repo/support2.sav"))

### Preprocessing
nrow(support2)
dat <- support2 |>
  mutate(surv = Surv(d.time, death)) |>
  select(all_of(c("surv", "num.co", "sex", "race", "scoma", "dzgroup", "ca",
                  "diabetes", "dementia", "age"))) %>%
  mutate(num.co = fct_lump_min(factor(num.co), 100),
         scoma = factor(scoma), age = sqrt(age)) %>%
  mutate_if(is.factor, ~ factor(.x, levels = levels(.x),
                                labels = make.names(levels(.x)),
                                ordered = FALSE)) %>%
  na.omit()
nrow(dat) # 43 patients removed due to missingness

# The set of all predictors is not invariant ------------------------------

summary(coxph(surv ~ scoma + dzgroup + ca + sex + race + num.co +
                age + diabetes + dementia, data = dat))

# Evidence of age and cancer being direct causes of time-to-death ---------

fm <- surv ~ scoma + dzgroup + ca + age + diabetes + dementia + sex + race
fe <- ~ num.co
if (file.exists(file.path(bpath, "gcm.rds"))) {
  gcm <- readRDS(file.path(bpath, "gcm.rds"))
} else {
  gcm <- coxphICP(formula = fm, data = dat, env = fe, test = "gcm.test")
  saveRDS(gcm, file.path(bpath, "gcm.rds"))
}
gcm
(wald <- coxphICP(formula = fm, data = dat, env = fe, type = "wald"))

# Multiple environments ---------------------------------------------------

tfe <- ~ num.co + race
tfm <- surv ~ scoma + dzgroup + ca + age + diabetes + dementia + sex
if (file.exists(file.path(bpath, "tgcm.rds"))) {
  tgcm <- readRDS(file.path(bpath, "tgcm.rds"))
} else {
  tgcm <- coxphICP(formula = tfm, data = dat, env = tfe, test = "gcm.test")
  saveRDS(tgcm, file.path(bpath, "tgcm.rds"))
}
tgcm
(twald <- coxphICP(formula = tfm, data = dat, env = tfe, type = "wald"))


# Incorporating prior knowledge about direct causes -----------------------

ffm <- surv ~ scoma + dzgroup + ca + sex + race
fmand <- ~ age + diabetes + dementia
if (file.exists(file.path(bpath, "fgcm.rds"))) {
  fgcm <- readRDS(file.path(bpath, "fgcm.rds"))
} else {
  fgcm <- coxphICP(formula = ffm, data = dat, env = fe, test = "gcm.test",
                         mandatory = fmand)
  saveRDS(fgcm, file.path(bpath, "fgcm.rds"))
}
fgcm
(fwald <- coxphICP(formula = ffm, data = dat, env = fe, type = "wald",
                   mandatory = fmand))

# Create Table ------------------------------------------------------------

tab <- lapply(list(gcm, wald, tgcm, twald, fgcm, fwald),
              \(mod) {
                tmp <- pvalues(gcm)
                tmp[] <- NA
                pvals <- pvalues(mod)
                tmp[names(pvals)] <- pvals
                tmp
              }) |> bind_rows()
tab
knitr::kable(tab, booktabs = TRUE, format = "latex", digits = 3L)

# Create Figure -----------------------------------------------------------

pdat <- data.frame(set = names(pvalues(gcm, "set")),
                   gcm = pvalues(gcm, "set"),
                   wald = pvalues(wald, "set")) |>
  mutate(caage = grepl("age", set) & grepl("ca", set),
         invariant = gcm > 0.05 | wald > 0.05)

p1 <- ggplot(pdat, aes(x = gcm, y = wald, shape = caage)) +
  geom_point(size = 2, alpha = 0.8) +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = 0, ymax = 1, fill = "darkblue", alpha = 0.1) +
  annotate("rect", ymin = 0, ymax = 0.05, xmin = 0, xmax = 1, fill = "darkred", alpha = 0.1) +
  theme_bw() +
  labs(x = "TRAM-GCM p-value", y = "TRAM-Wald p-value", shape = "Set contains both ca and age") +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 4),
                     labels = c("TRUE" = "yes", "FALSE" = "no")) +
  theme(legend.position = "top", text = element_text(size = 13.5)) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  coord_cartesian(xlim = c(0, 0.25), ylim = c(0, 0.25))

p2 <- ggplot(pdat, aes(x = gcm, y = wald, shape = caage)) +
  geom_point(size = 2, alpha = 0.4) +
  theme_bw() +
  annotate("rect", xmin = 0, xmax = 0.05, ymin = 0, ymax = 1, fill = "darkblue", alpha = 0.1) +
  annotate("rect", ymin = 0, ymax = 0.05, xmin = 0, xmax = 1, fill = "darkred", alpha = 0.1) +
  labs(x = "TRAM-GCM p-value", y = "TRAM-Wald p-value", shape = "Set contains both ca and age") +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 4),
                     labels = c("TRUE" = "yes", "FALSE" = "no")) +
  theme(legend.position = "top", text = element_text(size = 12.5)) +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  scale_x_log10() + scale_y_log10()

library("ggpubr")
ggarrange(p1, p2, common.legend = TRUE)
ggsave("inst/figures/casestudy.pdf", height = 4, width = 7.5)

# Appendix: Informative censoring -----------------------------------------

### Informative right-censoring in random fractions
if (!file.exists(file.path(bpath, "censoring.csv"))) {
  probs <- seq(1, 9, 1)/10
  eventobs <- which(dat$surv[, 2] == 1)
  nsim <- 10
  set.seed(241068)
  out <- lapply(seq_along(probs), \(idx) {
    lapply(1:nsim, \(iter) {
      dat$nd <- dat$surv[, 2]
      to_censor <- sample(eventobs, floor(probs[idx] * length(eventobs)))
      dat$nd[to_censor] <- 0
      dat$nsurv <- Surv(dat$surv[, 1], dat$nd)
      cfm <- nsurv ~ scoma + dzgroup + ca + age + diabetes + dementia + sex + race
      cat(paste0("\nTotal fraction of censored observations is ", round(mean(dat$nd == 0) * 100, 2), "%\n"))
      cicp <- coxphICP(formula = cfm, data = dat, env = fe, type = "wald")
      tibble(frac = probs[idx], output = paste0(cicp$cand, collapse = "+"),
             bind_rows(pvalues(cicp, "set")))
    }) |> bind_rows()
  }) |> bind_rows()
  write_csv(out, file.path(bpath, "censoring.csv"))
} else {
  out <- read_csv(file.path(bpath, "censoring.csv"))
}

### Generate tables
o2 <- out |> group_by(frac) |> count(output) |>
  pivot_wider(names_from = output, values_from = n, values_fill = 0)
o2

o2 |> select(frac, Empty, ca, `ca+age`) |>
  knitr::kable(booktabs = TRUE, format = "latex", digits = 3L)
