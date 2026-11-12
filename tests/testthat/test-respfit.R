test_that("RespFit returns a finite Mr statistic", {
  s <- sim_data(200, 6, "D", seed = 51)
  f <- quiet_fit(s$rt, s$x)
  r <- RespFit(f, order = 2)
  expect_s3_class(r, "RespFit")
  expect_true(is.finite(r$Mr)); expect_gt(r$df, 0)
})

test_that("RespFit warns about missing data", {
  s <- sim_data(200, 6, "D", seed = 52, miss = .05)
  f <- quiet_fit(s$rt, s$x)
  expect_warning(RespFit(f, order = 2), "missing")
})

test_that("RespFit moments and covariance agree with full pattern enumeration", {
  # brute force for 5 items: probabilities of all 2^5 response patterns by quadrature
  s <- sim_data(200, 5, "Q", seed = 53)
  f <- quiet_fit(s$rt, s$x, model = "Q")
  r <- RespFit(f, order = 2)
  J <- 5; pl <- f$par.log; A <- f$A; W <- f$W; nq <- length(A)
  pats <- as.matrix(expand.grid(rep(list(0:1), J)))
  sA <- sqrt(exp(pl[3 * J + 1])); sV <- sqrt(exp(pl[3 * J + 2]))
  pp <- numeric(nrow(pats))
  for (a in 1:nq) for (b in 1:nq) {
    p <- 1 / (1 + exp(-(exp(A[a] * sA) / exp(pl[1:J])) * (exp(A[b] * sV) / exp(pl[J + 1:J]))))
    pp <- pp + W[a] * W[b] * apply(pats, 1, function(x) prod(p^x * (1 - p)^(1 - x)))
  }
  X <- diffIRT:::fitItems(J, 2)[-1, ]
  pi_bf <- apply(X, 1, function(S) sum(pp[apply(pats[, S == 1, drop = FALSE], 1, all)]))
  expect_equal(r$Z$pred, round(pi_bf, 3))
})

test_that("RespFit is well calibrated under the true model and the df count identified parameters", {
  skip_on_cran()
  # df: moments minus identified combinations (D: 2J+1, Q: J+1 without constraints)
  for (m in c("D", "Q")) {
    s <- sim_data(300, 6, m, seed = 54)
    r <- RespFit(quiet_fit(s$rt, s$x, model = m), 2)
    expect_equal(r$df, 6 + 15 - if (m == "D") 13 else 7, info = m)
  }
  # mean of Mr over replications close to df (chi-square)
  Mr <- sapply(1:20, function(k) { s <- sim_data(400, 6, "D", seed = 100 + k)
    RespFit(quiet_fit(s$rt, s$x), 2)$Mr })
  expect_lt(abs(mean(Mr) - 8), 3 * sqrt(2 * 8 / 20) + 1)
})

test_that("RespFit handles 20+ items quickly", {
  skip_on_cran()
  s <- sim_data(300, 22, "D", seed = 55)
  f <- quiet_fit(s$rt, s$x)
  t <- system.time(r <- RespFit(f, 2))[["elapsed"]]
  expect_lt(t, 30)
  expect_equal(nrow(r$Z), 22 + choose(22, 2))
  expect_true(is.finite(r$Mr))
})

test_that("RespFit works with constraints", {
  s <- sim_data(300, 6, "D", seed = 56)
  r <- RespFit(quiet_fit(s$rt, s$x, constrain = "ai.equal"), 2)
  expect_true(is.finite(r$Mr))
  expect_equal(r$df, 21 - 8)                     # 1 common a + 6 v + 2 omegas, minus the free drift scale
})
