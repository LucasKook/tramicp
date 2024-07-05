### Environment is a parent and faithfulness is violated
### in BoxCox model can lead to non-empty output but requires
### highly tuned coefficients
# LK 2024

set.seed(1123)
library("tram")

g <- \(z) qchisq(pnorm(z), df = 5)
gi <- \(y) qnorm(pchisq(y, df = 5))

# DGP ---------------------------------------------------------------------

dgp <- function(bY, n = 3e3) {
  E <- rnorm(n)
  X <- E + rnorm(n)
  Y <- g(X + rnorm(n))
  SR <- 1/sqrt(2) * (gi(Y) - E)
  nY <- g(sqrt(2) * (bY * E + SR))
  LR <- X - (-0.721 + 0.115 * Y + 0.624 * E)
  nX <- - 0.721 + 0.115 * Y + 0.624 * E + LR
  data.frame(Y = nY, X = nX, E = E)
}

# Run ---------------------------------------------------------------------

cfs <- data.frame(bY = seq(-0.2, 0.2, length.out = 5))
out <- apply(cfs, 1, \(cfx) {
  res <- lapply(1:20, \(iter) {
    d <- dgp(1/sqrt(2) + cfx)
    ret <- try(
      tramicp::BoxCoxICP(Y ~ X, data = d, env = ~ E)$candidate_causal_predictors
    )
    if (inherits(ret, "try-error"))
      return(NA)
    ret
  })
  sum(res == "X", na.rm = TRUE)
})

# Output ------------------------------------------------------------------

knitr::kable(
  as.data.frame(cbind(cfs, output = out)), format = "latex", booktabs = TRUE
)
