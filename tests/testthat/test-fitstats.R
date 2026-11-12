test_that("RTFit expected cell probabilities agree with a large simulation from the model", {
  s <- sim_data(300, 5, "D", seed = 71)
  f <- quiet_fit(s$rt, s$x)
  r <- RTFit(f, bins = 5)
  expect_s3_class(r, "RTFit")
  expect_equal(dim(r$items), c(5, 4))
  e <- r$cells[[2]]$exp
  expect_equal(sum(e), 1, tolerance = 1e-8)
  expect_equal(e[1:5] + e[6:10], rep(.2, 5), tolerance = 2e-3)   # bins are predicted quintiles
  # simulate 100000 persons from the fitted parameters and classify item 2
  set.seed(1)
  sim <- diffIRT:::simPersons(f$par.log, 5, 1, 1e5)
  b <- findInterval(sim$rt[, 2], r$edges[[2]]) + 1
  mc <- tabulate(sim$x[, 2] * 5 + b, 10) / 1e5
  expect_lt(max(abs(mc - e)), 4e-3)
})

test_that("RTFit computes the covariance of the estimates when the fit has no standard errors", {
  s <- sim_data(200, 4, "D", seed = 72)
  f1 <- quiet_fit(s$rt, s$x, se = TRUE)
  f0 <- f1; f0$hessian <- NULL
  expect_equal(diffIRT:::fitVcov(f0), vcov(f1), tolerance = 1e-6)
  expect_equal(RTFit(f0)$overall, RTFit(f1)$overall, tolerance = 1e-5)
})

test_that("RTFit detects an item with contaminated response times, and not the others", {
  skip_on_cran()
  s <- sim_data(500, 8, "Q", seed = 73, miss = .1)
  i <- sample(500, 150); s$rt[i, 3] <- s$rt[i, 3] + .3
  r <- RTFit(quiet_fit(s$rt, s$x, model = "Q"))
  expect_lt(r$items$p[3], .001)
  expect_lt(r$overall[["p"]], .05)
  expect_true(all(r$items$df <= 19 & r$items$df >= 17))   # near-singular directions are dropped
  expect_output(print(r), "Overall")
})

test_that("RTFit p-values are roughly uniform under the model", {
  skip_on_cran()
  p <- unlist(lapply(1:6, function(k) {
    s <- sim_data(400, 6, "D", seed = 80 + k, miss = .1)
    RTFit(quiet_fit(s$rt, s$x))$items$p
  }))
  expect_gt(mean(p), .35); expect_lt(mean(p), .65)
  expect_lt(mean(p < .05), .15)
})

test_that("PersonFit is reproducible, handles missing patterns, and flags inattentive persons", {
  s <- sim_data(300, 10, "D", seed = 74, miss = .1)
  g <- 1:15
  m <- is.na(s$x)
  s$rt[g, 6:10] <- 2 * s$rt[g, 6:10]; s$x[g, 6:10] <- rbinom(75, 1, .5); s$x[m] <- NA
  f <- quiet_fit(s$rt, s$x)
  a <- PersonFit(f, nsim = 300, seed = 1); b <- PersonFit(f, nsim = 300, seed = 1)
  expect_identical(a, b)
  pf <- a$persons
  expect_equal(nrow(pf), 300)
  expect_true(all(pf$p > 0 & pf$p <= 1))
  expect_equal(pf$items, rowSums(!is.na(f$score)))
  expect_gt(mean(pf$p[g] < .05), .5)
  expect_lt(mean(pf$p[-g] < .05), .12)
  expect_output(print(a), "Lowest p-values")
})
