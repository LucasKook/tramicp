# Faithfulness in additive linear TRAMs
# LK 2023

# DEPs --------------------------------------------------------------------

devtools::load_all()

# DGP ---------------------------------------------------------------------

### DGP with faithfulness violations for cf = 0.5
dgp <- function(n, cf) {
  ### Environment
  E <- sample(c(-1, 1), n, TRUE)
  ### True invariant sets: X1 and (X1, X2)
  X1 <- E + rnorm(n)
  X2 <- E + rnorm(n)
  ### For cf = 0.5, we have cancellation and the empty set becomes invariant
  Y <- qchisq(pnorm(rnorm(n) - 0.5 * X1 + cf * X2), df = 3)
  ### Return
  data.frame(Y = Y, X1 = X1, X2 = X2, E = E)
}

# RUN ---------------------------------------------------------------------

cfs <- 0.5 + -3:3/10
ret <- structure(lapply(cfs, \(bx1) {
  ### Generate data with fixed random seed
  set.seed(1)
  ### and vary coef around 0.5
  d <- dgp(1e4, bx1)
  ### Run tramicp with BoxCox
  BoxCoxICP(Y ~ X1 + X2, data = d, env = ~ E, prob = c(0.001, 0.999),
            verbose = TRUE, test = "gcm.test")
}), names = cfs)

# OUT ---------------------------------------------------------------------

(o1 <- dplyr::bind_rows(lapply(ret, pvalues, which = "set"), .id = "bx1"))
o1$output <- unlist(lapply(ret, \(x) paste0(x$cand, collapse = "+")))

### Generate latex table
knitr::kable(o1, booktabs = TRUE, format = "latex", digits = 3)
