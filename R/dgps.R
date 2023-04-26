# Data generating processes

#' Simple data-generating process for illustrating tramicp
#'
#' @param n Sample size
#' @param K Number of outcome classes or order of Bernstein polynomial
#' @param nenv Number of environments
#' @param ge Environment specific effect
#' @param ae Environment specific effect
#' @param bx3 Effect of Y on X3
#' @param mod Type of model
#' @param interacting Toggle baseline interaction with env
#' @param rm_censoring Remove censoring from simulated responses
#' @param cfb Baseline coefs
#' @param cfx Shift coefs
#' @param bx2x1 coef from x2 to x1
#'
#' @return \code{data.frame} with simulated data
#'
#' @details Simulates from X2 -> X1 -> Y -> X3, with E affecting X1, X2, X3, but
#'     not Y.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#'
#' # Weibull
#' d <- dgp_dicp(mod = "weibull", nenv = 3)
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Survreg, w = rep(1, nrow(d)))
#'
#' # Colr
#' d <- dgp_dicp(mod = "colr")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Colr, w = rep(1, nrow(d)))
#'
#' # Binary
#' d <- dgp_dicp(mod = "binary", nenv = 7)
#' glmICP(Y ~ X1 + X2 + X3, data = d, env = ~ E, family = "binomial")
#'
#' # Cotram
#' d <- dgp_dicp(mod = "cotram")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = cotram::cotram, w = rep(1, nrow(d)))
#'
#' # Polr
#' d <- dgp_dicp(mod = "polr")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Polr, w = rep(1, nrow(d)))
#'
#' # Coxph
#' d <- dgp_dicp(mod = "coxph")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Coxph, w = rep(1, nrow(d)))
#'
#' # Normal linear model, either lm or Lm
#' d <- dgp_dicp(mod = "lm")
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Lm, w = rep(1, nrow(d)))
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = lm, w = rep(1, nrow(d)))
#' dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E,
#'      modFUN = Lm, w = rep(1, nrow(d)), type = "wald", test = "wald")
#'
dgp_dicp <- function(
  n = 1e3, K = 6, nenv = 2, bx3 = rnorm(1), ge = rnorm(nenv), ae = rnorm(nenv),
  mod = "polr", interacting = FALSE, rm_censoring = TRUE, cfb = c(-3, 1.35),
  cfx = rnorm(2), bx2x1 = rnorm(1)
) {

  # Environments
  E <- t(sapply(rep(1, n), rmultinom, size = 1, prob = rep(1, nenv)))
  tge <- E %*% ge
  tae <- E %*% ae

  # Exogenous
  ex1 <- rnorm(n, sd = 0.5^2)
  ex2 <- rnorm(n, sd = 0.5^2)
  ex3 <- rnorm(n, sd = 0.5^2)

  # Endogenous
  x2 <- tae + ex2
  x1 <- 0.5 * tge + bx2x1 * x2 + ex1

  # Data frame
  df <- data.frame(
    X1 = x1, X2 = x2, E = factor(apply(E, 1, which.max))
  )

  # Response
  if (mod == "binary") {
    ncY <- as.numeric(cbind(x1, x2) %*% cfx >= rlogis(n))
    Y <- factor(ncY)
    x3 <- tge + c(-0.5, 0.5)[1 + ncY] + ex3
  } else if (mod == "slm") {
    Y <- ncY <- - cfb[1] / cfb[2] + cbind(x1, x2) %*% cfx / cfb[2] + rnorm(n)
    x3 <- tge + Y / sd(Y) + ex3
  } else {
    m <- .tram_from_name(mod, ia = interacting, cfb = cfb, cfx = cfx)
    Y <- simulate(m, newdata = df)
    if (!is.ordered(Y)) {
      if (rm_censoring)
        Y <- .R2vec(Y)
      if (mod == "cotram")
        Y <- as.integer(Y)
      x3 <- .R2vec(Y) * bx3 / (2 * sd(Y)) + tge + ex3
    } else {
      x3 <- as.numeric(Y) / (2 * sd(as.numeric(Y))) * bx3 + tge + ex3
    }
  }

  # Return
  ret <- data.frame(Y = Y, X1 = x1, X2 = x2, X3 = x3,
                    E = factor(apply(E, 1, which.max)))

  structure(ret, "tae" = tae, "tge"= tge, "bx3" = bx3)

}

.R2vec <- as.double
