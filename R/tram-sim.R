# Generic transformation models for simulation

#' Setting up a general transformation model
#'
#' @param type Response type
#' @param distr Inverse link
#' @param basis Basis type
#' @param K Number of outcome classes or basis functions
#' @param supp Support of response
#' @param bounds Bounds of response
#' @param extrapolate See \code{\link[basefun]{Bernstein_basis}}
#' @param log_first See \code{\link[basefun]{Bernstein_basis}}
#' @param numeric_vars Names of numeric variables
#' @param cfx Shift coefs
#' @param factor_vars Names of factor variables
#' @param factor_levels Levels of factors
#' @param interacting Toggle varying baseline
#' @param ia_fac Scale coef
#' @param cfb Baseline coef
#' @param add See \code{\link[tram]{tram}}
#'
#' @return Object of class \code{"mlt"}
#' @export
#'
#' @examples
#' # Polr model
#' m <- .tram_setup()
#' plot(m, newdata = data.frame(X1 = 0, X2 = 0), type = "trafo", col = 1)
#'
#' # Weibull model
#' m <- .tram_setup(type = "survival", cfb = c(-1,1), basis = "log")
#' plot(m, newdata = data.frame(X1 = 0, X2 = 0), type = "trafo", col = 1)
#'
#' # Cox model
#' m <- .tram_setup(type = "survival", basis = "bernstein")
#' plot(m, newdata = data.frame(X1 = -3, X2 = -3), type = "distribution", col = 1)
#'
#' # Cotram
#' m <- .tram_setup(type = "count", basis = "bernstein")
#' plot(m, newdata = data.frame(X1 = -1, X2 = -1), type = "distribution", col = 1)
#'
#' # BoxCox
#' m <- .tram_setup(type = "continuous", basis = "bernstein", distr = "Normal")
#' plot(m, newdata = data.frame(X1 = -1, X2 = -1), type = "distribution", col = 1)
#'
#' # Small check
#' m <- .tram_setup(cfx = c(0, 0.5), type = "count", basis = "bernstein")
#' nd <- data.frame(X1 = runif(tn <- 1e4), X2 = runif(tn),
#'                  E = factor(1 + rbinom(tn, 1, 0.5)))
#' nd$Y <- as.integer(tramicp:::.R2vec(simulate(m, newdata = nd)))
#' confint(cotram::cotram(Y ~ X1 + X2, data = nd, prob = 0.999))
#'
#' @import variables
#' @import mlt
#'
.tram_setup <- function(
  type = c("ordered", "continuous", "survival", "count"),
  distr = c("Logistic", "MinExtrVal", "MaxExtrVal", "Normal"),
  basis = c("bernstein", "linear", "log"),
  K = 6,
  supp = .supp(type),
  bounds = .bounds(type),
  extrapolate = FALSE,
  log_first = FALSE,
  numeric_vars = c("X1", "X2"),
  cfb = c(-2, -1),
  cfx = c(0.6, -0.3),
  factor_vars = "E",
  factor_levels = c(1, 2),
  interacting = FALSE,
  ia_fac = 0.7,
  add = c(0, 0)
) {

  # Match args
  type <- match.arg(type)
  distr <- match.arg(distr)
  basis <- match.arg(basis)

  # Response basis
  if (type == "ordered")
    y <- ordered_var(name = "Y", levels = 1:K)
  else
    y <- numeric_var(name = "Y", bounds = bounds, support = supp, add = add)
  ybasis <- .rbasis(y = y, type = type, basis = basis, distr = distr, K = K,
                    interacting = interacting, extrapolate = extrapolate,
                    log_first = log_first, cfb = cfb)
  coefs <- ybasis$coefs

  # Covariate basis
  Xs <- do.call("c", lapply(numeric_vars, numeric_var))

  # non-prop
  if (interacting) {
    ia <- factor_var(factor_vars, levels = factor_levels)
    iab <- as.basis(~ 0 + E, data = ia, remove_intercept = TRUE)
    coefs <- c(coefs, ia_fac * coefs)
  } else {
    iab <- NULL
  }

  m <- ctm(
    response = ybasis[["basis"]],
    interacting = iab,
    shifting = as.basis(
      as.formula(paste0("~", paste0(numeric_vars, collapse = "+"))), data = Xs,
      remove_intercept = TRUE
    ),
    todistr = distr, data = df
  )

  ncfs <- names(coef(m))
  cfa <- c(coefs, cfx)
  names(cfa) <- ncfs
  coef(m) <- cfa

  m

}

# Helpers

.supp <- function(type) {
  switch(
    type, "ordered" = NA, "continuous" = c(-10, 10),
    "survival" = c(1.1 * .Machine$double.eps, 10),
    "count" = c(1.1 * .Machine$double.eps, 100)
  )
}

.bounds <- function(type) {
  switch(
    type, "ordered" = NA, "continuous" = c(-Inf, Inf),
    "survival" = c(1.1 * .Machine$double.eps, Inf),
    "count" = c(1.1 * .Machine$double.eps, Inf),
  )
}

#' @import basefun
.rbasis <- function(y, type, basis, distr, K, interacting, extrapolate,
                    log_first, cfb,
                    prmin = 0.0002, # ifelse(type == "ordered", 0.05, 0.0002),
                    prmax = 1 - prmin) {
  if (type == "ordered")
    ybasis <- as.basis(~ 0 + Y, data = y)
  else
    ybasis <- switch(
      basis,
      "linear" = as.basis(~ Y, data = y),
      "log" = log_basis(y),
      "bernstein" = Bernstein_basis(y, order = K - 1, ui = "increasing",
                                    extrapolate = extrapolate,
                                    log_first = log_first)
    )

  fam <- .distr(distr)

  coefs <- fam$q(seq(prmin, prmax, length.out = K))
  if (type == "ordered")
    coefs <- fam$q(seq(prmin, prmax, length.out = K + 1))[-1]
  if (basis %in% c("linear", "log"))
    coefs <- cfb

  list(basis = ybasis, coefs = coefs)

}

.tram_from_name <- function(model = c("polr", "colr", "cotram", "coxph",
                                      "weibull", "lm", "boxcox"),
                            nvars = c("X1", "X2"), cfb = c(-3, 1.35),
                            cfx = c(0.6, -0.3), ia = FALSE, K = 6, tadd = 50) {
  switch(
    match.arg(model),
    "polr" = .tram_setup(type = "ordered", distr = "Logistic", K = K,
                         numeric_vars = nvars, cfb = cfb, cfx = cfx,
                         interacting = ia),
    "colr" = .tram_setup(type = "continuous", distr = "Logistic", K = K,
                         basis = "bernstein", numeric_vars = nvars, cfb = cfb,
                         cfx = cfx, interacting = ia, add = c(-tadd, tadd)),
    "cotram" = .tram_setup(type = "count", distr = "Logistic", K = K,
                           basis = "bernstein", numeric_vars = nvars, cfb = cfb,
                           cfx = cfx, interacting = ia, add = c(0, tadd)),
    "coxph" = .tram_setup(type = "survival", distr = "MinExtrVal", K = K,
                          basis = "bernstein", numeric_vars = nvars, cfb = cfb,
                          interacting = ia, extrapolate = FALSE, cfx = cfx,
                          log_first = FALSE, add = c(0, tadd)),
    "weibull" = .tram_setup(type = "survival", distr = "MinExtrVal", K = K,
                            basis = "log", numeric_vars = nvars, cfb = cfb,
                            cfx = cfx, interacting = ia, extrapolate = FALSE,
                            log_first = FALSE, add = c(0, tadd),
                            bounds = c(-Inf, Inf)),
    "lm" = .tram_setup(type = "continuous", distr = "Normal", K = K,
                       basis = "linear", numeric_vars = nvars, cfb = cfb,
                       cfx = cfx, interacting = ia, add = c(-2, 2)),
    "boxcox" = .tram_setup(type = "continuous", distr = "Normal", K = K,
                           basis = "bernstein", numeric_vars = nvars, cfb = cfb,
                           cfx = cfx, interacting = ia, add = c(-tadd, tadd)),
  )
}


# from {mlt} {
.Normal <- function()
  list(parm = function(x) NULL,
       p = pnorm, d = dnorm, q = qnorm,
       ### see also MiscTools::ddnorm
       dd = function(x) -dnorm(x = x) * x,
       ddd = function(x) dnorm(x = x) * (x^2 - 1),
       dd2d = function(x) -x,
       call = ".Normal",
       name = "normal")

.Exponential <- function()
  list(parm = function(x) NULL,
       p = pexp, d = dexp, q = qexp,
       dd = function(x) -dexp(x = x),
       ddd = function(x) dexp(x = x),
       dd2d = function(x) -1,
       call = ".Exponential",
       name = "exponential")

.Logistic <- function()
  list(parm = function(x) NULL,
       p = plogis, d = dlogis, q = qlogis,
       dd = function(x) {
         ex <- exp(x)
         (ex - ex^2) / (1 + ex)^3
       },
       ddd = function(x) {
         ex <- exp(x)
         (ex - 4 * ex^2 + ex^3) / (1 + ex)^4
       },
       dd2d = function(x) {
         ex <- exp(x)
         (1 - ex) / (1 + ex)
       },
       call = ".Logistic",
       name = "logistic")

### Gompertz distribution
.MinExtrVal <- function()
  list(parm = function(x) NULL,
       p = function(x, lower.tail = TRUE, log.p = FALSE) {
         ### p = 1 - exp(-exp(x))
         ret <- exp(-exp(x))
         if (log.p) {
           if (lower.tail)
             return(log1p(-ret))
           return(-exp(x))
         }
         if (lower.tail)
           return(1 - exp(-exp(x)))
         return(ret)
       },
       q = function(p) log(-log1p(- p)),
       d = function(x, log = FALSE) {
         ret <- x - exp(x)
         if (!log) return(exp(ret))
         ret
       },
       dd = function(x) {
         ex <- exp(x)
         (ex - ex^2) / exp(ex)
       },
       ddd = function(x) {
         ex <- exp(x)
         (ex - 3*ex^2 + ex^3) / exp(ex)
       },
       dd2d = function(x)
         1 - exp(x),
       call = ".MinExtrVal",
       name = "minimum extreme value")

### Gumbel distribution
.MaxExtrVal <- function()
  list(parm = function(x) NULL,
       p = function(x, lower.tail = TRUE, log.p = FALSE) {
         ### p = exp(-exp(-x))
         if (log.p) {
           if (lower.tail)
             return(-exp(-x))
           return(log1p(-exp(-exp(-x))))
         }
         if (lower.tail)
           return(exp(-exp(-x)))
         1 - exp(-exp(-x))
       },
       q = function(p) -log(-log(p)),
       d = function(x, log = FALSE) {
         ret <- - x - exp(-x)
         if (!log) return(exp(ret))
         ret
       },
       dd = function(x) {
         ex <- exp(-x)
         exp(-ex - x) * (ex - 1)
       },
       ddd = function(x) {
         ex <- exp(-x)
         exp(-x - ex) * (ex - 1)^2 - exp(-ex - 2 * x)
       },
       dd2d = function(x)
         exp(-x) - 1,
       call = ".MaxExtrVal",
       name = "maximum extreme value")

.distr <- function(which = c("Normal", "Logistic", "MinExtrVal", "MaxExtrVal",
                             "Exponential")) {
  which <- match.arg(which)
  do.call(paste(".", which, sep = ""), list())
}

# }
