test_that("analytic gradient matches central differences (D and Q, missing data, constraints)", {
  for (m in c("D", "Q")) for (miss in c(0, .1)) {
    s <- sim_data(150, 6, m, seed = 3, miss = miss); J <- 6
    p_full <- true_par(s, m) + rnorm(3 * J + 2, 0, .02)
    # free parameters
    o <- objective(s, m)
    expect_lt(max(grad_rel_error(o$L, o$G(p_full), p_full)), 1e-4)
    # equal a's, Ter[1] fixed, common variance
    cons <- c(rep(1, J), 1 + 1:J, 0, J + 1 + 1:(J - 1), 2 * J + 1, 2 * J + 1)
    start <- rep(NA, 3 * J + 2); start[2 * J + 1] <- p_full[2 * J + 1]
    o2 <- objective(s, m, constrain = cons, start = start)
    p <- sapply(1:max(cons), function(k) p_full[which(cons == k)[1]])
    expect_lt(max(grad_rel_error(o2$L, o2$G(p), p)), 1e-4)
  }
})

test_that("adaptive-quadrature gradient matches central differences, including capped grids", {
  for (m in c("D", "Q")) {
    s <- sim_data(150, 6, m, seed = 4, miss = .1); J <- 6
    p <- true_par(s, m) + rnorm(3 * J + 2, 0, .02)
    o <- objective(s, m)
    aq <- o$moments(p)
    expect_lt(max(grad_rel_error(function(q) o$L(q, aq), o$G(p, aq), p)), 1e-4)
    aqc <- aq; h <- 1:75                                      # grids wider than the prior: capped
    aqc$sA[h] <- 3 * sqrt(exp(p[3 * J + 1])); aqc$sV[h] <- 3 * sqrt(exp(p[3 * J + 2]))
    expect_lt(max(grad_rel_error(function(q) o$L(q, aqc), o$G(p, aqc), p)), 1e-4)
  }
})

test_that("numeric (forward-difference) gradient is close to the analytic one", {
  s <- sim_data(150, 6, "D", seed = 5); p <- true_par(s, "D")
  o <- objective(s, "D")
  ga <- o$G(p); gn <- o$G(p, grad = "numeric")
  expect_lt(max(abs(ga - gn) / pmax(1, abs(ga))), 1e-2)
})
