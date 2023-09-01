
test_that("RANGER works as inteded", {
  ### Continuous Y
  # RANGER with empty set
  set.seed(1)
  tn <- 1e2
  d <- dgp_dicp(n = tn, mod = "lm")
  rf <- RANGER(Y ~ 1, data = d)
  expect_length(residuals(rf), tn)
  expect_length(residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5]), 5)
  # RANGER with non-empty set
  rf <- RANGER(Y ~ X1, data = d)
  res <- residuals(rf)
  reso <- d$Y - predict(rf, data = d)$predictions
  expect_equal(res, reso)
  resn <- residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5])
  resno <- d$Y[1:5] - predict(rf, data = d[1:5, ])$predictions
  expect_equal(resn, resno)
  # RANGER GCM
  expect_no_error(dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = "RANGER"))

  ### Binary Y
  # RANGER with empty set
  d <- dgp_dicp(n = tn, mod = "binary")
  rf <- RANGER(Y ~ 1, data = d)
  expect_length(residuals(rf), tn)
  expect_length(residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5]), 5)
  # RANGER with non-empty set
  rf <- RANGER(Y ~ X1, data = d)
  res <- residuals(rf)
  reso <- as.numeric(d$Y) - 1 - predict(rf, data = d)$predictions[, 2]
  expect_equal(res, reso)
  resn <- residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5])
  resno <- as.numeric(d$Y[1:5]) - 1 - predict(rf, data = d[1:5, ])$predictions[, 2]
  expect_equal(resn, resno)
  # RANGER GCM
  expect_no_error(dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = "RANGER"))

  ### Ordinal Y
  # RANGER with empty set
  d <- dgp_dicp(n = tn, mod = "polr")
  rf <- RANGER(Y ~ 1, data = d)
  expect_length(residuals(rf), tn)
  expect_length(residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5]), 5)
  # RANGER with non-empty set
  rf <- RANGER(Y ~ X1, data = d)
  res <- residuals(rf)
  reso <- as.numeric(d$Y) - predict(rf, data = d)$predictions %*% seq_len(max(d$Y))
  expect_equal(res, reso)
  resn <- residuals(rf, newdata = d[1:5, ], newy = d$Y[1:5])
  resno <- as.numeric(d$Y[1:5]) - predict(rf, data = d[1:5, ])$predictions %*% seq_len(max(d$Y))
  expect_equal(resn, resno)
  # RANGER GCM
  expect_no_error(dicp(Y ~ X1 + X2 + X3, data = d, env = ~ E, modFUN = "RANGER"))
})
