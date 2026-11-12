test_that("diffIRT recovers known parameters (D and Q)", {
  skip_on_cran()
  for (m in c("D", "Q")) {
    s <- sim_data(400, 10, m, seed = 21)
    f <- quiet_fit(s$rt, s$x, model = m)
    J <- 10
    expect_equal(f$conv, 0)
    expect_lt(max(abs(f$par[1:J] / s$ai - 1)), .15)                 # a
    expect_gt(cor(f$par[J + 1:J], s$vi), .95)                         # v
    expect_lt(max(abs(f$par[2 * J + 1:J] - s$ter)), .05)              # Ter
    expect_lt(abs(f$par[3 * J + 1] - .3), .06)                        # omega_gamma (unbiased with adaptive quadrature)
  }
})

test_that("standard errors: vcov is solve(H/2) and SEs are its delta-method transform", {
  data(brightness, package = "diffIRT")
  x <- brightness[, 1:12]; rt <- brightness[, 13:24]
  cb <- c(rep(1, 6), rep(2, 6), 3:8, 3:8, rep(9, 12), 0, 10)
  f <- quiet_fit(rt, x, constrain = cb, start = c(rep(NA, 36), 0, NA), se = TRUE)
  V <- vcov(f)
  expect_equal(V, solve(.5 * f$hessian))
  k <- match(1:10, cb); nit <- 12
  jac <- rep(1, 3 * nit + 2)                                          # d reported / d internal
  jac[c(1:nit, (2 * nit + 1):(3 * nit))] <- f$par[c(1:nit, (2 * nit + 1):(3 * nit))]
  jac[3 * nit + 1:2] <- f$par[3 * nit + 1:2] / 2
  expect_equal(f$std.err[k], abs(jac[k]) * sqrt(diag(V)), tolerance = 1e-8)
  expect_equal(f$std.err[37], -999)                                   # fixed parameter
  expect_true(all(is.finite(f$std.err)))
})

test_that("constraint shortcuts and the old (non-adaptive, numeric-gradient) settings run", {
  skip_on_cran()
  s <- sim_data(200, 6, "D", seed = 22)
  for (cn in c("ai.equal", "vi.equal", "ter.equal", "scale.equal")) {
    f <- quiet_fit(s$rt, s$x, constrain = cn)
    expect_equal(f$conv, 0, info = cn)
  }
  f1 <- quiet_fit(s$rt, s$x, control = list(adaptive = FALSE, grad = "numeric", eps = .01))
  f2 <- quiet_fit(s$rt, s$x, control = list(adaptive = FALSE))
  expect_lt(max(abs(f1$par - f2$par)), .02)
})

test_that("Ter starting values are below every item's fastest response", {
  for (seed in 1:6) {
    s <- sim_data(300, 10, "D", seed = seed); J <- 10
    st <- diffIRT:::get.start(s$rt, s$x, 1:(3 * J + 2), rep(NA, 3 * J + 2), 1)
    expect_true(all(exp(st[2 * J + 1:J]) < apply(s$rt, 2, min)), info = seed)
  }
})

test_that("a fit that would stop with Ter collapsed to 0 is rescued", {
  skip_on_cran()
  # with the original starting values, Ter[3] here started above the fastest RT and collapsed
  set.seed(3); s <- simdiff(500, 20, model = "D")
  f <- quiet_fit(s$rt, s$x, se = TRUE)
  expect_gt(min(f$par[41:60] / s$ter), .9)
  expect_false(any(is.na(f$std.err)))
})
