test_that("predicted RT density is finite, non-negative and integrates to 1, even on a fine grid", {
  for (m in c("D", "Q")) {
    s <- sim_data(200, 5, m, seed = 61)
    f <- quiet_fit(s$rt, s$x, model = m)
    pl <- f$par.log; J <- 5; q <- quad(61)
    for (j in 1:2) {
      dens <- function(y, A = f$A, W = f$W) diffIRT:::mdiffIRT(y, pl[j], pl[J + j], pl[2 * J + j], pl[3 * J + 1], pl[3 * J + 2],
                                                              A, W, f$model, f$control$eps)
      y <- seq(exp(pl[2 * J + j]) + 1e-6, 3 * max(s$rt[, j]), length = 300)
      d61 <- dens(y, q$A, q$W)
      expect_true(all(is.finite(d61)) && all(d61 >= 0), label = paste(m, j))       # was Inf / -999
      expect_equal(integrate(dens, exp(pl[2 * J + j]), Inf, subdivisions = 500)$value, 1, tolerance = 1e-3)
      expect_lt(max(abs(dens(y) - d61)) / max(d61), 2e-3)                         # nq = 7 is accurate here
    }
  }
})

test_that("QQdiff runs and returns quantiles", {
  s <- sim_data(150, 4, "D", seed = 62)
  f <- quiet_fit(s$rt, s$x)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  r <- QQdiff(f, items = 1:2, quant = 20)
  expect_equal(dim(r$qexp), c(20, 2))
  expect_true(all(is.finite(r$qexp)))
  expect_true(all(diff(r$qexp[, 1]) > 0))
})

test_that("predicted quantiles and CDF agree with a high-precision reference", {
  for (m in c("D", "Q")) {
    s <- sim_data(300, 4, m, seed = 63)
    f <- quiet_fit(s$rt, s$x, model = m)
    pl <- f$par.log; J <- 4; j <- 2
    args <- list(pl[j], pl[J + j], pl[2 * J + j], pl[3 * J + 1], pl[3 * J + 2])
    dens <- function(y) do.call(diffIRT:::mdiffIRT, c(list(y), args, list(f$A, f$W, f$model, f$control$eps)))
    lo <- exp(pl[2 * J + j]); mx <- 2 * max(s$rt[, j])
    P <- function(y) integrate(dens, lo, y, rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000)$value
    pr <- c(.05, .2, .5, .8, .95)
    ref <- sapply(pr, function(p) uniroot(function(y) P(y) - p, c(lo + 1e-9, mx), tol = 1e-12)$root)
    q <- do.call(diffIRT:::qdiffIRT, c(list(pr), args, list(mx, f$A, f$W, f$model, f$control$eps)))
    expect_lt(max(abs(q - ref)), 1e-6, label = m)
    y <- lo + c(.05, .3, 1, 2)
    cdf <- do.call(diffIRT:::pdiffIRT, c(list(y), args, list(f$A, f$W, f$model, f$control$eps)))
    expect_lt(max(abs(cdf - sapply(y, P))), 1e-7, label = m)
  }
})

test_that("qdiffIRT extends maxRT when the CDF has not reached the requested probability", {
  s <- sim_data(150, 4, "D", seed = 64)
  f <- quiet_fit(s$rt, s$x)
  pl <- f$par.log; J <- 4
  q <- function(mx) diffIRT:::qdiffIRT(c(.5, .95), pl[1], pl[J + 1], pl[2 * J + 1], pl[3 * J + 1], pl[3 * J + 2],
                                       mx, f$A, f$W, f$model, f$control$eps)
  expect_equal(q(exp(pl[2 * J + 1]) + .05), q(2 * max(s$rt[, 1])), tolerance = 1e-7)
})

test_that("QQdiff defaults to at most 100 quantiles and is fast for large N", {
  skip_on_cran()
  s <- sim_data(1000, 3, "D", seed = 65)
  f <- quiet_fit(s$rt, s$x)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_equal(nrow(QQdiff(f, items = 1, plot = 1)$qexp), 100)
  t <- system.time(r <- QQdiff(f, items = 1:3, quant = 1000))[["elapsed"]]
  expect_lt(t, 5)
  expect_true(all(diff(r$qexp[, 1]) > 0))
})
