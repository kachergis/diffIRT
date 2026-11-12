test_that("non-adaptive factest scores are stationary points of LLfacA / LLfacV", {
  for (m in c("D", "Q")) {
    s <- sim_data(60, 8, m, seed = 31, miss = .1)
    f <- quiet_fit(s$rt, s$x, model = m)
    fs <- factest(f, control = list(adaptive = FALSE))   # LLfacA/V use the non-adaptive grid
    rt <- f$rt; x <- f$score; rt[is.na(rt) | is.na(x)] <- -999; x[is.na(x)] <- -999
    pl <- f$par.log; J <- 8; N <- 60; mod <- f$model
    run <- function(fun, P) .C(fun, as.double(rt), as.double(x), as.integer(N), as.integer(J), pl[1:J], pl[J + 1:J],
      pl[3 * J + 1], pl[3 * J + 2], pl[2 * J + 1:J], as.double(P), as.integer(mod), f$A, f$W, as.integer(length(f$A)),
      f$control$eps, -1, double(N), PACKAGE = "diffIRT")[[17]]
    h <- 1e-3
    la <- log(fs[, 1])
    gA <- (log(run("LLfacA", la + h)) - log(run("LLfacA", la - h))) / (2 * h)
    lv <- if (m == "D") fs[, 2] else log(fs[, 2]); tv <- function(z) if (m == "D") z else exp(z)
    gV <- (log(run("LLfacV", tv(lv + h))) - log(run("LLfacV", tv(lv - h)))) / (2 * h)
    expect_lt(max(abs(gA)), .01, label = paste(m, "gamma"))
    expect_lt(max(abs(gV)), .01, label = paste(m, "theta"))
  }
})

test_that("factest standard errors are finite and positive", {
  s <- sim_data(80, 8, "D", seed = 32, miss = .1)
  f <- quiet_fit(s$rt, s$x)
  fs <- factest(f, se = TRUE)
  expect_equal(dim(fs), c(80, 4))
  expect_true(all(is.finite(fs)))
  expect_true(all(fs[, c(2, 4)] > 0))
})

test_that("factest with adaptive grids matches a 61-point reference grid", {
  for (m in c("D", "Q")) {
    s <- sim_data(200, 15, m, seed = 33)
    f <- quiet_fit(s$rt, s$x, model = m)
    ref <- f; q <- quad(61); ref$A <- q$A; ref$W <- q$W; ref$aq <- NULL
    fr <- factest(ref, se = TRUE)
    fa <- factest(f, se = TRUE)
    fn <- factest(f, se = TRUE, control = list(adaptive = FALSE))
    expect_lt(max(abs(fa[, c(1, 3)] - fr[, c(1, 3)])), 2e-3, label = m)
    expect_lt(max(abs(fa[, c(2, 4)] / fr[, c(2, 4)] - 1)), .02, label = m)
    expect_gt(max(abs(fn[, 3] - fr[, 3])), 5 * max(abs(fa[, 3] - fr[, 3])))  # non-adaptive is worse
  }
})

test_that("factest falls back to the non-adaptive grid for non-adaptive fits", {
  s <- sim_data(80, 8, "D", seed = 34)
  f <- quiet_fit(s$rt, s$x, control = list(adaptive = FALSE))
  expect_null(f$aq)
  fs <- factest(f, se = TRUE)
  expect_true(all(is.finite(fs)))
})
