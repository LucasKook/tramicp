#!/usr/bin/env Rscript
# Visualize simulation results
# LK 2023

settings <- c("main", "app", "hidden", "link", "wald-extended")

setting <- settings[as.numeric(commandArgs(TRUE))[1]]
if (is.na(setting))
  setting <- settings[1]

# deps --------------------------------------------------------------------

library("tidyverse")
devtools::load_all()
devtools::load_all("inst/helpers")

theme_set(theme_bw() + theme(
  legend.position = "top", strip.background = element_rect(
    fill = "transparent")))

# FUNs --------------------------------------------------------------------

.setdiff <- function(x, y) {
  x <- x[x != "empty"]
  y <- y[y != "empty"]
  ret <- try(setdiff(x, y))
  if (inherits(ret, "try-error"))
    return(character(0))
  ret
}

.intersect <- function(x, y) {
  ret <- try(intersect(x, y))
  if (inherits(ret, "try-error"))
    return(character(0))
  ret
}

okabe_ito <- function(n) {
  unname(c(orange = "#E69F00",
           green = "#009E73",
           # amber = "#F5C710",
           # blue = "#0072B2",
           red = "#D55E00",
           purple = "#CC79A7",
           `light blue` = "#56B4E9",
           `light blue` = "#56B4E9",
           grey = "#999999",
           black = "#000000")[1:n])
}

# params ------------------------------------------------------------------

mods <- c("Binary" = "binary", "Cotram" = "cotram", "Weibull" = "weibull")
amods <- c("BoxCox" = "boxcox", "Colr" = "colr", "Coxph" = "coxph", "Lm" = "lm",
           "Polr" = "polr")
mtests <- c("TRAM-Wald" = "wald", "TRAM-GCM" = "gcm.test", "RCIT" = "KCI")
tests <- c(mtests, "Correlation" = "cor.test", "ROC" = "roc.test", "TRAM-Wald extended" = "wald.test")
cols <- okabe_ito(length(tests))
names(cols) <- names(tests)

if (setting == "main") {
  ### Paths
  bpath <- "inst/results/results_main"
  tmods <- mods
  ttests <- tests[-4]
  out <- "sim-main.pdf"
  th <- 4
  tw <- 7

  fobs <- readRDS(file.path(bpath, "fobs.rds"))
  res <- read_sim_results(file.path(bpath, "results"))
  res$oicp[res$oicp == "empty"] <- "Empty"

  bpath <- "inst/results/results_binary-roc"
  bfobs <- readRDS(file.path(bpath, "fobs.rds"))
  bres <- read_sim_results(file.path(bpath, "results"))
  bres$oicp[bres$oicp == "empty"] <- "Empty"
  bres <- dplyr::filter(bres, test == "cor.test")
  bres$test <- "roc.test"

  res <- full_join(res, bres)
} else if (setting == "app") {
  bpath <- "inst/results/results_appendix"
  tmods <- amods
  ttests <- mtests
  out <- "sim-app.pdf"
  th <- 4
  tw <- 7
  ### path
  fobs <- readRDS(file.path(bpath, "fobs.rds"))
  res <- read_sim_results(file.path(bpath, "results"))
  res$oicp[res$oicp == "empty"] <- "Empty"
} else if (setting == "hidden") {
  bpath <- "inst/results/results_hidden"
  tmods <- c(mods, amods)
  ttests <- mtests
  out <- "sim-hidden.pdf"
  th <- 4
  tw <- 10
  ### path
  fobs <- readRDS(file.path(bpath, "fobs.rds"))
  res <- read_sim_results(file.path(bpath, "results"))
  res$oicp[res$oicp == "empty"] <- "Empty"
  res$oicp <- "Empty"
} else if (setting == "link") {
  bpath <- "inst/results/results_link"
  tmods <- c(mods, amods)
  ttests <- mtests
  out <- "sim-link.pdf"
  th <- 4
  tw <- 10
  ### path
  fobs <- readRDS(file.path(bpath, "fobs.rds"))
  res <- read_sim_results(file.path(bpath, "results"))
  res$oicp[res$oicp == "empty"] <- "Empty"
} else if (setting == "wald-extended") {
  ### paths
  bpath <- "inst/results/results_main"
  tmods <- c(mods[-1], amods)
  ttests <- tests
  out <- "sim-wald-extended.pdf"
  th <- 4
  tw <- 10

  fobs <- readRDS(file.path(bpath, "fobs.rds"))
  res <- read_sim_results(file.path(bpath, "results"))
  res$oicp[res$oicp == "empty"] <- "Empty"
  res <- dplyr::filter(res, test == "wald")

  bpath <- "inst/results/results_wald-extended"
  bfobs <- readRDS(file.path(bpath, "fobs.rds"))
  bres <- read_sim_results(file.path(bpath, "results"))
  bres$oicp[bres$oicp == "empty"] <- "Empty"
  bres <- dplyr::filter(bres, test == "wald")
  bres$test <- "wald.test"

  res <- full_join(res, bres)

  bbpath <- "inst/results/results_appendix"
  bbfobs <- readRDS(file.path(bbpath, "fobs.rds"))
  bbres <- read_sim_results(file.path(bbpath, "results"))
  bbres$oicp[bbres$oicp == "empty"] <- "Empty"
  bbres <- dplyr::filter(bbres, test == "wald")
  bbres$test <- "wald"

  res <- full_join(res, bbres)

}

### read
prdat <- res %>%
  mutate(inv = case_when(inv == "1" ~ "Empty", TRUE ~ inv)) %>%
  gather("method", "set", inv, oicp) %>%
  mutate(
    splpaY = str_split(paY, "\\+"),
    splset = str_split(set, "\\+"),
    numpa = unlist(map2(splpaY, splset, ~ as.numeric(length(.intersect(.x, .y)) > 0))),
    nonpa = unlist(map2(splset, splpaY, ~ as.numeric(length(.setdiff(.x[.x != "Empty"], .y)) > 0))),
    empty = as.numeric(set == "Empty"),
    jaccard = unlist(map2(splpaY, splset, ~ length(.intersect(.x, .y)) / length(union(.x, .y)))),
    fwer = unlist(map2(splpaY, splset, ~ as.numeric(length(setdiff(.y[.y != "Empty"], .x)) > 0))),
  )

prdat %>% group_by(n, mod, type, test) %>%
  summarize(jaccard = mean(jaccard), fwer = mean(fwer > 0))

if (setting == "hidden") {
  anY <- paste0("X1+X2+X3", unlist(lapply(fobs$dags, \(x) {
    ret <- names(which(x$dag["Y", c("X4", "X5")] == 0))
    if (identical(ret, character(0)))
      ret <- ""
    else
      ret <- paste0("+", paste0(ret, collapse = "+"))
    ret
  })))
  prdat <- prdat %>% dplyr::select(-paY) %>% group_by(dag) %>% nest() %>% arrange(dag) %>%
    ungroup() %>% mutate(paY = anY) %>% unnest(data)
}

# Vis ---------------------------------------------------------------------

p1 <- vis(tmods, ttests)
ggsave(file.path("inst", "figures", out), p1, height = th, width = tw, scale = 1)
