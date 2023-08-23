## ROC-based invariance test
## LK 2023

set.seed(27)
parall <- TRUE # Parallel computing?

# Dependencies ------------------------------------------------------------

library("SimDesign")
library("tidyverse")
library("future.apply")
devtools::load_all()
devtools::load_all("inst/helpers")
library("colorspace")
library("RCIT")
library("pROC")
theme_set(theme_bw() + theme(legend.position = "top"))

save <- TRUE # Save simulation results in current R session
store <- FALSE # Store simulation results in folder
fixed <- TRUE # Same DAGs for all models

spec <- "correct"

# Params ------------------------------------------------------------------

# Models
mods <- "binary"
lmods <- "binary"
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

# Params
ndags <- ifelse(spec == "correct", 100, 50) # Number of DAGs
nsim <- 20 # Number of repetitions
ncores <- min(nsim, 20) # Number of cores
ns <- c(1e2, 3e2, 1e3, 3e3) # , 1e4, 3e4) # Sample sizes
blfix <- TRUE # fixed baseline transformation
nanc <- 3 # ancestors of Y
panc <- 0.8 # edge probability
if (spec == "hidden") panc <- 1
ndec <- 2 # descendants of Y
pdec <- 0.8 # edge probability
nenv <- 2 # number of environments
penv <- 0.8 # edge probability (0 for Y, only applies to an(Y) and de(Y))
tK <- 6 # Dimension of Bernstein
polrK <- 6 # Dimension of ordinal outcome
resp <- "Y" # Name of response var
# Build formula
preds <- paste0("X", 1:(nanc + ndec))
if (spec == "hidden") preds <- preds[-1]
fml <- paste0(resp, "~", paste0(preds, collapse = "+"))
env <- "E" # Name of environment var
cfb <- c(-22, 8) # baseline cf
rmc <- FALSE # Whether to remove censoring
stdz <- TRUE # Whether to standardize
sde <- sqrt(30)
errDistAnY <- "normal"
errDistDeY <- "normal"
mixAnY <- 0.1
mixDeY <- 0.1

dags <- if (fixed) {
  lapply(1:ndags, \(iter) {
    random_dag(nenv = nenv, nanc = nanc, ndec = ndec, penv = penv,
               panc = panc, pdec = pdec, cfb = cfb,
               cfe = .rcfx(nanc + ndec, penv, FALSE, sde),
               plot_graph = FALSE)
  })
} else NULL

fobs <- list(types = types, fml = fml, resp = resp, env = env, preds = preds,
             blfix = blfix, cfb = cfb, K = tK, polrK = polrK, rmc = rmc,
             stdz = stdz, nsim = nsim, sde = sde, errDistAnY = errDistAnY,
             errDistDeY = errDistDeY, mixAnY = mixAnY, mixDeY = mixDeY,
             dags = dags, mods = mods, tests = tests, spec = spec)

if (save) {
  pvec <- c("nanc", nanc, "ndec", ndec, "panc", panc, "pdec", pdec, "penv", penv,
            "stdz", stdz, "sde", round(sde, 2), "errDist", errDistAnY, errDistDeY,
            "spec", spec)
  outdir <- file.path(
    "inst", "results", Sys.Date(), "rdag-all", paste0(
      paste0(names(pvec), pvec, collapse = ""), collapse = "_"))
  if (!dir.exists(outdir))
    dir.create(outdir, recursive = TRUE)

  write_rds(fobs, file.path(outdir, "fobs.rds"))
}

# Sim ---------------------------------------------------------------------

Design <- tibble(expand_grid(n = ns, mod = mods, dag = 1:ndags, nanc = nanc,
                             panc = panc, ndec = ndec, pdec = pdec, nenv = nenv,
                             penv = penv))

# Run ---------------------------------------------------------------------

res <- runSimulation(
  design = Design,
  replications = nsim,
  generate = if (fixed) GENFIX else GEN,
  analyse = AUCANA,
  summarise = NA,
  save_results = save,
  store_results = store,
  debug = "none",
  parallel = parall,
  ncores = ncores,
  fixed_objects = fobs,
  packages = c("magrittr", "tramicp", "pROC"),
  filename = "sim-results.rds",
  save_details = list(
    safe = TRUE, save_results_dirname = file.path(outdir, "results"),
    save_seeds_dirname = file.path(outdir, "seeds")
  )
)
