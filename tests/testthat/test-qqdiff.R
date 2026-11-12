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
