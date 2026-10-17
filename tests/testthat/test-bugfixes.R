test_that("predictor p-values match terms exactly", {
  pv <- c("Empty" = 0.01, "X1" = 0.01, "X10" = 0.9, "X1+X10" = 0.01)
  ipv <- tramicp:::.indiv_pvals(c("X1", "X10"), pv, alpha = 0.05)
  expect_equal(unname(ipv), c(0.9, 0.01))
})

test_that("combined test p-value is at most one", {
  set.seed(123)
  r <- rnorm(100)
  e <- factor(rep(1:2, each = 50))
  tst <- tramicp:::.combined_test(r, e, NULL)
  expect_lte(tst$p.value, 1)
})

test_that("custom test works with multivariable environments", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "boxcox")
  ctest <- \(r, e, controls) list(p.value = 0.5)
  expect_no_error(dicp(Y ~ X1, d, ~ E + X2, modFUN = "BoxCox", test = ctest,
                       verbose = FALSE))
})

test_that("Wald test for empty set uses all environments", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "boxcox")
  d$E2 <- factor(sample(1:3, nrow(d), TRUE))
  res <- dicp(Y ~ X1, d, ~ E + E2, modFUN = "BoxCox", type = "wald",
              verbose = FALSE)
  expect_equal(nrow(res$tests[[1]]$linfct), 3)
})

test_that("greedy ICP keeps last test and handles failed tests", {
  set.seed(2)
  d <- dgp_dicp(n = 500, mod = "boxcox")
  res <- dicp(Y ~ X1 + X2, d, ~ E, modFUN = "BoxCox", type = "wald",
              verbose = FALSE, greedy = TRUE)
  expect_true(all(c("X1", "X2") %in% names(pvalues(res, "set"))))

  ### Wald test fails for larger models, resulting in NA p-values
  ctrl <- dicp_controls("wald", "wald", vcov = function(object) {
    if (length(coef(object)) > 3) stop("vcov failed")
    vcov(object)
  })
  expect_no_error(suppressMessages(capture.output(
    res <- dicp(Y ~ X1 + X2, d, ~ E, modFUN = "BoxCox", type = "wald",
                controls = ctrl, verbose = FALSE, greedy = TRUE),
    type = "message")))
  expect_type(res$candidate_causal_predictors, "character")
})

test_that("invariant_sets ignores NA p-values", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "boxcox")
  res <- dicp(Y ~ X1 + X2, d, ~ E, modFUN = "BoxCox", type = "wald",
              verbose = FALSE)
  res$set_pvals[2] <- NA
  expect_false(anyNA(invariant_sets(res)))
})

test_that("Wald test only uses environment coefficients", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "boxcox")
  d$AGE <- rnorm(nrow(d))
  res <- dicp(Y ~ AGE, d, ~ E, modFUN = "BoxCox", type = "wald",
              verbose = FALSE)
  expect_setequal(rownames(res$tests[[2]]$linfct), c("E2", "AGE:E2"))
})

test_that("failed Wald test for empty set is labelled", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "boxcox")
  ctrl <- dicp_controls("wald", "wald",
                        vcov = function(object) stop("vcov failed"))
  capture.output(res <- dicp(Y ~ X1, d, ~ E, modFUN = "BoxCox", type = "wald",
                             controls = ctrl, verbose = FALSE),
                 type = "message")
  expect_equal(names(pvalues(res, "set")), c("Empty", "X1"))
})

test_that("GCM regresses environment on all transformed terms", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "lm")
  set.seed(1)
  r1 <- lmICP(Y ~ I(X1) + I(X2), d, ~ E, verbose = FALSE)
  set.seed(1)
  r2 <- lmICP(Y ~ X1 + X2, d, ~ E, verbose = FALSE)
  expect_equal(r1$tests[[4]]$rEX, r2$tests[[4]]$rEX)
})

test_that("ranger residuals on new data use new response", {
  set.seed(123)
  d <- dgp_dicp(n = 100, mod = "lm")
  train <- d[1:50, ]
  test <- d[51:100, ]
  rf <- tramicp:::RANGER(Y ~ X1, data = train)
  expect_equal(residuals(rf, newdata = test),
               unname(test$Y - predict(rf, data = test)$predictions))
  rf0 <- tramicp:::RANGER(Y ~ 1, data = train)
  expect_equal(residuals(rf0, newdata = test), test$Y - mean(train$Y))
})

test_that("dgp_dicp respects K", {
  set.seed(123)
  d <- dgp_dicp(n = 300, mod = "polr", K = 3)
  expect_equal(nlevels(d$Y), 3)
})

test_that("predictor p-values are one if all sets are rejected", {
  set.seed(3)
  d <- dgp_dicp(n = 1e3, mod = "boxcox")
  res <- dicp(Y ~ X3, d, ~ E, modFUN = "BoxCox", type = "wald",
              verbose = FALSE)
  expect_true(all(pvalues(res, "set") < 0.05))
  expect_equal(unname(pvalues(res, "predictor")), 1)
  expect_equal(res$candidate_causal_predictors, "Empty")
})
