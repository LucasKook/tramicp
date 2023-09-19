# Experiments with random DAGs

#' Generating data from a random DAG
#'
#' @param n Integer; sample size.
#' @param dag DAG as returned from \code{random_dag()}
#' @param ia Interacting
#' @param standardize Standardize ancestral and descendental data
#' @param plot_model Toggle; plot model
#' @param errDistAnY Error distribution for ancestral DAG
#' @param mixAnY Mixture proportion for ancestral DAG
#' @param errDistDeY Error distribution for descendental DAG
#' @param mixDeY Mixture proportion for deccendental DAG
#'
#' @importFrom pcalg randomDAG rmvDAG possAn
#' @importFrom Rgraphviz plot
#'
#' @return Simulated data from random DAG.
#'
#' @export
#'
#' @examples
#' set.seed(21)
#' nanc <- ndec <- 2
#' dag <- random_dag(nanc = nanc, ndec = ndec, penv = 1)
#' d <- dgp_random_dag(n = 1e3, mod = "colr", dag = dag)
#' attr(d, "paY")
#' attr(d, "oracle_icp")
#' fml <- Y ~ X1 + X2 + X3 + X4
#' ColrICP(fml, data = d, env = ~ E, type = "wald")
#' ColrICP(fml, data = d, env = ~ E, type = "residual")
#'
#' mods <- c("polr", "boxcox", "cotram", "colr", "weibull", "lm")
#' sapply(mods, \(m) dgp_random_dag(mod = m))
#'
#' set.seed(1)
#' d <- dgp_random_dag(mod = "weibull", plot = TRUE, standardize = TRUE)
#' attr(d, "cfb")
#'
dgp_random_dag <- function(
  n = 100, mod = "polr", K = 6, ia = FALSE, rm_censoring = TRUE,
  standardize = FALSE, plot_model = FALSE, errDistAnY = "normal",
  mixAnY = 0.1, errDistDeY = "normal", mixDeY = 0.1, dag = random_dag()
) {

  # Unpack
  nanc <- dag$nanc
  ndec <- dag$ndec
  nenv <- dag$nenv
  pa_dag <- dag$pa_dag
  ch_dag <- dag$ch_dag
  cfr <- dag$cfr
  cfe <- dag$cfe
  cfx <- dag$cfx
  cfb <- dag$cfb
  cfc <- dag$cfc

  ## Environment
  E <- matrix(sample(0:1, n, TRUE), nrow = n, ncol = nanc + ndec)

  ## Ancestor DAG
  ddat <- rmvDAG(n, pa_dag, errDist = errDistAnY, mix = mixAnY)
  ddat <- ddat + E[, 1:nanc] * matrix(dag$cfe[1:nanc], nrow = n,
                                      ncol = nanc, byrow = TRUE)
  if (standardize) ddat <- scale(ddat)
  dat <- data.frame(ddat)

  ## Response
  if (mod == "binary") {
    ncY <- as.numeric(ddat %*% cfx >= rlogis(n))
    dat$Y <- factor(ncY)
  } else if (mod == "slm") {
    ncY <- ddat %*% cfx + rnorm(n)
    dat$Y <- ncY
  } else {
    tm <- tramicp:::.tram_from_name(
      model = mod, cfx = cfx, nvars = colnames(dat), K = K, ia = ia, cfb = cfb
    )
    if (plot_model) plot(tm, newdata = dat, K = 3e2)
    Y <- simulate(tm, newdata = dat)
    ncY <- .R2vec(Y)
    if (mod == "cotram") Y <- ncY <- as.integer(ncY)
    dat$Y <- if (rm_censoring) ncY else Y
    if (mod == "polr") dat$Y <- ordered(ncY, levels = 1:K)
  }

  ## Descendant DAG
  ch_dat <- rmvDAG(n, ch_dag, errDist = errDistDeY, mix = mixDeY)
  ch_dat <- ch_dat + matrix(as.numeric(ncY), nrow = n, ncol = ndec) *
    matrix(cfc, nrow = n, ncol = ndec, byrow = TRUE) + E[, (nanc + 1):(nanc + ndec)] *
    matrix(cfe[-1:-nanc], nrow = n, ncol = ndec) + ddat %*% cfr
  if (standardize) ch_dat <- scale(ch_dat)
  ch_dat <- data.frame(ch_dat)
  colnames(ch_dat) <- paste0("X", (nanc + 1):(nanc + ndec))

  ## Return
  structure(data.frame(dat, ch_dat, E = factor(E[, 1])), mod = mod, pa_dag = pa_dag,
            ch_dag = ch_dag, cfx = cfx, cfc = cfc, cfb = cfb, cfe = cfe,
            paY = dag$paY, oracle_icp = dag$oracle_icp, chE = dag$chE, dag = dag)
}

#' Title
#'
#' @param nenv Number of environments
#' @param nanc Number of potential ancestors
#' @param ndec Number of potential descendants
#' @param penv Edge probability for environment
#' @param panc Edge probability for ancestral graph
#' @param pdec Edge probability for descendental graph
#' @param cfx Coefficients from X to Y
#' @param cfc Coefficients from Y to children
#' @param cfb Coefficients of baseline trafo in case of (log-) linear basis
#' @param cfe Coefficients from E to X
#' @param plot_graph Toggle, plot DAG
#' @param BanYdeY Coefficient matrix for ancestral to descendental variables
#' @param AanYdeY Masking matrix for 0 coefs in \code{BanYdeY}
#'
#' @return List containing random DAG and infos for simulation experiments
#' @export
#'
#' @examples
#' random_dag()
#'
random_dag <- function(
    nenv = 2, nanc = 3, ndec = 3, penv = 0.3, panc = 0.3, pdec = 0.3,
    cfx = .rcfx(nanc, panc, ensure_one = TRUE), cfc = .rcfx(ndec, pdec),
    cfb = c(-6, 4), cfe = .rcfx(nanc + ndec, penv), plot_graph = FALSE,
    BanYdeY = matrix(rnorm(nanc * ndec, sd = sqrt(0.5)), nrow = nanc, ncol = ndec),
    AanYdeY = matrix(sample(0:1, nanc * ndec, replace = TRUE), nrow = nanc, ncol = ndec)
) {

  ## Ancestor DAG
  pa_dag <- randomDAG(nanc, panc)

  ## Descendant DAG
  ch_dag <- randomDAG(ndec, pdec)

  ## Oracle ICP
  paY <- .chkpa(.addx(which(cfx != 0)))
  chE <- .addx(which(cfe != 0))
  invisible(capture.output(
    pAn <- possAn(as(pa_dag, "matrix"), which(cfx != 0), type = "dag")
  ))
  anY <- .addx(pAn)
  anY_int_chE <- intersect(chE, anY)
  if (!identical(anY_int_chE, character(0)) && anY_int_chE[1] != "X") {
    adj <- .dag2adj(pa_dag, nms = .addx(1:nanc))
    tmp <- .ch_from_node(adj, anY_int_chE)
  } else tmp <- "X"
  pa_anY_intersect_chE <- .chkpa(tmp)
  oicp <- .chkpa(intersect(paY, union(chE, pa_anY_intersect_chE)))

  ## Reconstruct full DAG adj
  pdm <- as(pa_dag, "matrix")
  pdc <- as(ch_dag, "matrix")
  colnames(pdm) <- rownames(pdm) <- .addx(1:nanc)
  colnames(pdc) <- rownames(pdc) <- .addx((nanc + 1):(nanc + ndec))
  p2c <- BanYdeY * AanYdeY
  rownames(p2c) <- .addx(1:nanc)
  colnames(p2c) <- .addx((nanc + 1):(nanc + ndec))
  c2p <- t(p2c)
  c2p[] <- 0

  fdag <- rbind(cbind(pdm, p2c), cbind(c2p, pdc))
  fdag <- rbind(cbind(fdag, Y = c(cfx, rep(0, ndec))), Y = c(rep(0, nanc), cfc, 0))
  fdag <- rbind(cbind(fdag, E = 0), E = c(cfe, 0, 0))

  # Plot full DAG
  if (plot_graph) Rgraphviz::plot(as(abs(fdag), "graphNEL"))

  ## Return
  list(pa_dag = pa_dag, ch_dag = ch_dag, cfx = cfx, cfc = cfc, cfb = cfb,
       cfe = cfe, paY = paY, oracle_icp = oicp, chE = chE, dag = fdag,
       cfr = BanYdeY * AanYdeY, nanc = nanc, panc = panc, ndec = ndec,
       pdec = pdec, nenv = nenv, penv = penv)

}

# Helpers -----------------------------------------------------------------

.addx <- function(x) {
  paste0("X", x)
}

.chkpa <- function(x) {
  if (identical(x, "X")|identical(x, character(0)))
    "empty"
  else
    x
}

.ch_from_node <- function(adj, nodes) {
  rownames(adj)[which(adj[nodes, ] == 1 & t(adj[, nodes] == 0))]
}

.dag2adj <- function(dag, nms = paste0("X", dag@nodes)) {
    adj <- t(as(dag, "matrix")) != 0
    rownames(adj) <- colnames(adj) <- nms
    adj
}

.rcfx <- function(n, p, ensure_one = FALSE, sd = sqrt(0.9)) {
  cfx <- rnorm(n, sd = sd)
  idx <- rbinom(n, 1, p)
  if (ensure_one & all(idx == 0))
    idx[sample.int(n, 1)] <- 1
  cfx * idx
}
