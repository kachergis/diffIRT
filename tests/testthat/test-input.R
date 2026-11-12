test_that("invalid data give clear errors", {
  s <- sim_data(50, 5, "D", seed = 41)
  x2 <- s$x; x2[1, 1] <- 2
  expect_error(quiet_fit(s$rt, x2), "coded 0")
  rt2 <- s$rt; rt2[1, 1] <- -.5
  expect_error(quiet_fit(rt2, s$x), "positive")
  rt3 <- s$rt; rt3[2, ] <- NA
  expect_error(quiet_fit(rt3, s$x), "no observed responses")
  expect_error(quiet_fit(s$rt[, -1], s$x), "same size")
})

test_that("anova warns when the two fits use different data", {
  skip_on_cran()
  a <- sim_data(150, 5, "D", seed = 42); b <- sim_data(150, 5, "D", seed = 43)
  f0 <- quiet_fit(a$rt, a$x, constrain = "ai.equal"); f1 <- quiet_fit(a$rt, a$x); g1 <- quiet_fit(b$rt, b$x)
  expect_silent(utils::capture.output(anova(f0, f1)))
  expect_warning(utils::capture.output(tryCatch(anova(f0, g1), error = function(e) NULL)), "different")
})
