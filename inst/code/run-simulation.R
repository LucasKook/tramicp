#!/usr/bin/env Rscript
# Simulation setup for reproducing the results in Section 4 and Appendix B2
# Lucas Kook, 2023

set.seed(27) # Seed
parall <- FALSE # Parallel computing?
.libPaths(c("~/tutu/lib", .libPaths()))

### Command line args and defaults
args <- as.numeric(commandArgs(TRUE))
specs <- c("correct", "link", "hidden", "roc", "larger")
spec <- if (noargs <- identical(numeric(0), args)) specs[1] else specs[args[1]]
nsim <- if (noargs) 20 else args[2] # number of repetitions per DAG
ncores <- if (noargs) 20 else min(args[3], nsim) # number of cores for parallel
ndags <- if (noargs) 100 else args[4] # number of (random) DAGs
TEST <- if (noargs) 1 else args[5] # for testing this script
ROW <- if (noargs) NULL else args[6]

# Dependencies ------------------------------------------------------------

library("SimDesign")
library("tidyverse")
devtools::load_all()
devtools::load_all("inst/helpers")
library("colorspace")
theme_set(theme_bw() + theme(legend.position = "top"))

save <- TRUE # Save simulation results in current R session
store <- FALSE # Store simulation results in folder
fixed <- TRUE # Same DAGs for all models

# Params ------------------------------------------------------------------

# Models
mods <- c("binary", "cotram", "weibull", "coxph", "polr", "lm", "colr", "boxcox")
lmods <- c("binary", "cotram", "Survreg", "Coxph", "Polr", "Lm", "Colr", "BoxCox")
names(lmods) <- mods

# Types
types <- c("wald", "residual", "kci", "residual")
ltypes <- c("Wald test", "cor.test", "KCI", "gcm.test")
names(ltypes) <- types

# Test
tests <- c("wald", "cor.test", "KCI", "gcm.test")
names(types) <- tests
ltest <- c("Wald", "cor.test", "KCI", "gcm.test")
names(ltest) <- tests

# ROC-specific args
pkgs <- "magrittr"
tANA <- ANA
if (spec == "roc") {
  mods <- mods[1]
  lmods <- lmods[1]
  tANA <- AUCANA
  pkgs <- c(pkgs, "pROC")
}

# Params
ns <- c(1e2, 3e2, 1e3, 3e3, 1e4) # Sample sizes
if (spec %in% c("hidden", "larger")) ns <- c(ns, 3e4, 1e5)
blfix <- TRUE # fixed baseline transformation
nanc <- 3 # ancestors of Y
panc <- 0.8 # edge probability
if (spec %in% c("hidden", "larger")) panc <- 1
ndec <- 2 # descendants of Y
pdec <- 0.8 # edge probability
nenv <- 2 # number of environments
penv <- 0.8 # edge probability (0 for Y, only applies to an(Y) and de(Y))
tK <- 6 # Dimension of Bernstein
polrK <- 6 # Dimension of ordinal outcome
resp <- "Y" # Name of response var
# Build formula
preds <- paste0("X", 1:(nanc + ndec))
if (spec %in% c("hidden", "larger")) preds <- preds[-1]
fml <- paste0(resp, "~", paste0(preds, collapse = "+"))
env <- "E" # Name of environment var
cfb <- c(-22, 8) # baseline cf
rmc <- FALSE # Whether to remove censoring
stdz <- TRUE # Whether to standardize
sde <- sqrt(10)
errDistAnY <- "normal"
errDistDeY <- "normal"
mixAnY <- 0.1
mixDeY <- 0.1
affect_variance <- FALSE

dags <- if (fixed) {
  lapply(1:ndags, \(iter) {
    tcfx <- .rcfx(nanc, panc, FALSE, sd = sqrt(0.9))
    tcfc <- .rcfx(ndec, pdec, FALSE, sd = sqrt(0.3))
    tcfe <- .rcfx(nanc + ndec, penv, FALSE, sd = sde)
    if (spec == "larger") {
      tcfx[1] <- tcfx[1] * 2.5
      tcfx[-1] <- tcfx[-1] / 5
      tcfc <- tcfc
    }
    random_dag(nenv = nenv, nanc = nanc, ndec = ndec, penv = penv,
               panc = panc, pdec = pdec, cfb = cfb, cfx = tcfx,
               cfe = tcfe, cfc = tcfc, plot_graph = FALSE)
  })
} else NULL

fobs <- list(types = types, fml = fml, resp = resp, env = env, preds = preds,
             blfix = blfix, cfb = cfb, K = tK, polrK = polrK, rmc = rmc,
             stdz = stdz, nsim = nsim, sde = sde, errDistAnY = errDistAnY,
             errDistDeY = errDistDeY, mixAnY = mixAnY, mixDeY = mixDeY,
             dags = dags, mods = mods, tests = tests, spec = spec,
             affect_variance = affect_variance)

if (save) {
  pvec <- c("nanc", nanc, "ndec", ndec, "panc", panc, "pdec", pdec, "penv", penv,
            "stdz", stdz, "sde", round(sde, 2), "errDist", errDistAnY, errDistDeY,
            "spec", spec)
  outdir <- file.path(
    "inst", "results", Sys.Date(), paste0(paste0(
      names(pvec), pvec, collapse = ""), collapse = "_"))
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  write_rds(fobs, file.path(outdir, "fobs.rds"))
}

# Sim ---------------------------------------------------------------------

Design <- tibble(expand_grid(n = ns, mod = mods, dag = 1:ndags, nanc = nanc,
                             panc = panc, ndec = ndec, pdec = pdec, nenv = nenv,
                             penv = penv))

if (!is.null(ROW)) {
  stopifnot(ROW <= nrow(Design))
  Design <- Design[ROW, ]
}

if (TEST) {
  nsim <- 2
  Design <- Design[ROW <- 1, ]
}


# Run ---------------------------------------------------------------------

suppressWarnings(file.remove(list.files(pattern = "SIMDESIGN-TEMPFILE")))
res <- runSimulation(
  design = Design,
  replications = nsim,
  generate = if (fixed) GENFIX else GEN,
  analyse = tANA,
  summarise = NA,
  save = FALSE,
  save_results = save,
  store_results = store,
  debug = "none",
  parallel = parall,
  ncores = ncores,
  fixed_objects = fobs,
  packages = pkgs,
  filename = file.path(outdir, paste0("sim-results", ROW, ".rds")),
  save_details = list(
    safe = TRUE, save_results_dirname = file.path(outdir, paste0("results-row-", ROW)),
    save_seeds_dirname = file.path(outdir, "seeds")
  )
)
