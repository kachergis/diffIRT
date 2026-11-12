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

test_that("EAP scores and posterior SDs match a brute-force posterior on a fine grid", {
  for (m in c("D", "Q")) {
    s <- sim_data(150, 8, m, seed = 35, miss = .1)
    f <- quiet_fit(s$rt, s$x, model = m)
    e <- factest(f, se = TRUE, type = "EAP")
    expect_equal(colnames(e), c("gamma[p]", "gamma[p].se", "theta[p]", "theta[p].se"))
    expect_equal(factest(f, type = "EAP"), e[, c(1, 3)])
    rt <- f$rt; x <- f$score; rt[is.na(rt) | is.na(x)] <- -999; x[is.na(x)] <- -999
    pl <- f$par.log; J <- 8; G <- 161; u <- seq(-6, 6, length = G)
    for (i in 1:10) {
      zA <- f$aq$muA[i] + f$aq$sA[i] * rep(u, G); zV <- f$aq$muV[i] + f$aq$sV[i] * rep(u, each = G)
      ll <- diffIRT:::llPoints(pl, rt[i, , drop = FALSE], x[i, , drop = FALSE], 1, J, f$model, f$control$eps,
                               matrix(zA, 1), matrix(zV, 1))
      lw <- ll - zA^2 / (2 * exp(pl[3 * J + 1])) - zV^2 / (2 * exp(pl[3 * J + 2])); w <- exp(lw - max(lw)); w <- w / sum(w)
      g <- exp(zA); th <- if (m == "Q") exp(zV) else zV
      ref <- c(sum(w * g), sqrt(sum(w * g^2) - sum(w * g)^2), sum(w * th), sqrt(sum(w * th^2) - sum(w * th)^2))
      expect_equal(unname(e[i, ]), ref, tolerance = 1e-5, label = paste(m, i))
    }
  }
})

test_that("EAP works for non-adaptive fits and when a variance is fixed at 0", {
  s <- sim_data(100, 6, "D", seed = 36)
  f <- quiet_fit(s$rt, s$x, control = list(adaptive = FALSE))
  expect_true(all(is.finite(factest(f, se = TRUE, type = "EAP"))))
  J <- 6; con <- c(1:(3 * J + 1), 0); st <- c(rep(NA, 3 * J + 1), 0)
  f0 <- quiet_fit(s$rt, s$x, constrain = con, start = st)
  e <- factest(f0, type = "EAP")
  expect_true(all(is.finite(e)))
  expect_equal(unname(e[, 2]), rep(0, 100))                # theta_p fixed at 0
})

test_that("plausible values have the right shape, are reproducible, and vary around the EAP", {
  s <- sim_data(200, 8, "Q", seed = 37)
  f <- quiet_fit(s$rt, s$x, model = "Q")                   # no se: item parameters held fixed
  set.seed(1); pv <- factest(f, type = "PV", npv = 20)
  set.seed(1); pv2 <- factest(f, type = "PV", npv = 20)
  expect_identical(pv, pv2)
  expect_length(pv, 20)
  expect_false(attr(pv, "items"))
  expect_equal(dim(pv[[1]]), c(200, 2)); expect_equal(colnames(pv[[1]]), c("gamma[p]", "theta[p]"))
  expect_true(all(sapply(pv, function(v) all(is.finite(v) & v > 0))))
  e <- factest(f, type = "EAP", se = TRUE)
  m <- Reduce(`+`, pv) / 20
  z <- (m[, 2] - e[, 3]) / (e[, 4] / sqrt(20))            # PV means vs posterior means
  expect_lt(abs(mean(z)), 4 / sqrt(200)); expect_lt(abs(sd(z) - 1), .2)
  expect_error(factest(f, type = "PV", control = list(pv.items = TRUE)), "se = TRUE")
})

test_that("plausible values reproduce the estimated population variances", {
  skip_on_cran()
  s <- sim_data(400, 10, "D", seed = 38, miss = .1)
  f <- quiet_fit(s$rt, s$x, se = TRUE)
  set.seed(2); pv <- factest(f, type = "PV", npv = 20, control = list(pv.items = FALSE))
  v <- rowMeans(sapply(pv, function(p) c(var(log(p[, 1])), var(p[, 2]))))
  expect_equal(v, exp(f$par.log[31:32]), tolerance = .05)
  map <- factest(f)
  expect_lt(var(log(map[, 1])), .9 * exp(f$par.log[31]))  # point estimates are shrunk
  set.seed(3); pvi <- factest(f, type = "PV", npv = 2)
  expect_true(attr(pvi, "items"))                         # default with se = TRUE
})
