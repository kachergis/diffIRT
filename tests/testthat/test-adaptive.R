test_that("adaptive rule with the grid at mu = 0, s = sigma equals the non-adaptive rule", {
  for (m in c("D", "Q")) {
    s <- sim_data(100, 6, m, seed = 6, miss = .1); J <- 6; p <- true_par(s, m)
    o <- objective(s, m)
    aq0 <- list(muA = rep(0, 100), sA = rep(sqrt(exp(p[3 * J + 1])), 100),
                muV = rep(0, 100), sV = rep(sqrt(exp(p[3 * J + 2])), 100))
    expect_equal(o$L(p, aq0), o$L(p, NULL), tolerance = 1e-10)
  }
})

test_that("adaptive nq = 7 is far more accurate than non-adaptive nq = 7", {
  s <- sim_data(200, 10, "D", seed = 7); J <- 10; p <- true_par(s, "D")
  ref <- objective(s, "D", nq = 61)$L(p)
  o <- objective(s, "D", nq = 7)
  aq <- NULL; for (k in 1:3) aq <- o$moments(p, aq)
  expect_lt(abs(o$L(p, aq) - ref), .05)
  expect_gt(abs(o$L(p) - ref), 10 * abs(o$L(p, aq) - ref))
})

test_that("a variance fixed at zero gives finite results with adaptive quadrature", {
  s <- sim_data(100, 6, "D", seed = 8); J <- 6; p <- true_par(s, "D"); p[3 * J + 1] <- -746
  o <- objective(s, "D")
  aq <- o$moments(p)
  expect_true(is.finite(o$L(p, aq)))
  expect_true(all(is.finite(o$G(p, aq))))
})
