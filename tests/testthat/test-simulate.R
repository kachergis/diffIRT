test_that("simdiff with default arguments reproduces earlier versions exactly", {
  fx <- readRDS(test_path("fixtures", "simdiff_r.rds"))     # made with the diffIRT 1.5 simdiff
  for (m in c("D", "Q")) {
    set.seed(2024); expect_identical(simdiff(30, 4, model = m), fx[[paste0(m, "_default")]], label = m)
    g <- fx[[paste0(m, "_given")]]
    set.seed(8); expect_identical(do.call(simdiff, g$args), g$out, label = m)
  }
})

test_that("simdiff omega, rho and miss behave as specified", {
  set.seed(1)
  # (very wide person distributions can make the rejection sampler fail for extreme person-item pairs)
  s <- simdiff(3000, 2, model = "D", omega = c(.35, .6), rho = .5)
  expect_equal(sd(log(s$gamma)), .35, tolerance = .05)
  expect_equal(sd(s$theta), .6, tolerance = .05)
  expect_equal(cor(log(s$gamma), s$theta), .5, tolerance = .1)
  q <- simdiff(3000, 2, model = "Q", rho = -.4)
  expect_equal(cor(log(q$gamma), log(q$theta)), -.4, tolerance = .1)
  m <- simdiff(1000, 5, miss = .2)
  expect_equal(mean(is.na(m$rt)), .2, tolerance = .1)
  expect_identical(is.na(m$rt), is.na(m$x))
  expect_warning(simdiff(5, 2, ai = c(1, 1), vi = c(0, 0), gamma = rep(1, 5), theta = rep(0, 5), ter = c(.5, .5), rho = .3), "ignored")
  expect_error(simdiff(5, 2, rho = 2), "rho")
})

# a minimal stand-in for a fitted object with known parameters
fake_fit <- function(ai, vi, ter, oA, oV, model) {
  q <- quad(21)
  list(nit = length(ai), par = c(ai, vi, ter, oA, oV), model = if (model == "D") 1 else 2, A = q$A, W = q$W)
}

test_that("predict() matches Monte Carlo expectations (person-level and population)", {
  for (m in c("D", "Q")) {
    ai <- c(.6, .9); vi <- if (m == "D") c(-.6, .3) else c(1.2, 2); ter <- c(.5, .8)
    obj <- fake_fit(ai, vi, ter, .3, if (m == "D") .5 else .3, m)
    # person level: everyone has gamma = 1.2, theta = .4 (D) / 1.3 (Q)
    th <- if (m == "D") .4 else 1.3
    set.seed(11); s <- simdiff(6000, 2, ai = ai, vi = vi, gamma = rep(1.2, 6000), theta = rep(th, 6000), ter = ter, model = m)
    p <- diffIRT:::predict.diffIRT(obj, gamma = 1.2, theta = th)
    expect_equal(as.vector(p$prob), colMeans(s$x), tolerance = .03, label = paste(m, "prob"))
    expect_equal(as.vector(p$rt), colMeans(s$rt), tolerance = .03, label = paste(m, "rt"))
    # population level
    set.seed(12); s <- simdiff(8000, 2, ai = ai, vi = vi, ter = ter, model = m,
      gamma = exp(.3 * rnorm(8000)), theta = if (m == "D") .5 * rnorm(8000) else exp(.3 * rnorm(8000)))
    pp <- diffIRT:::predict.diffIRT(obj)
    expect_equal(pp$prob, unname(colMeans(s$x)), tolerance = .03, label = paste(m, "population prob"))
    expect_equal(pp$rt, unname(colMeans(s$rt)), tolerance = .03, label = paste(m, "population rt"))
  }
})

test_that("simulate() is reproducible, keeps the missing pattern, and fitted() has the right shape", {
  s <- sim_data(150, 5, "D", seed = 71, miss = .1)
  f <- quiet_fit(s$rt, s$x)
  a <- simulate(f, nsim = 2, seed = 5); b <- simulate(f, nsim = 2, seed = 5)
  expect_identical(a, b)
  expect_length(a, 2)
  expect_identical(is.na(a[[1]]$rt), is.na(f$rt) | is.na(f$score))
  expect_false(anyNA(simulate(f, nsim = 1, seed = 5, keep.missing = FALSE)[[1]]$rt))
  ft <- fitted(f)
  expect_equal(dim(ft$prob), c(150, 5)); expect_equal(dim(ft$rt), c(150, 5))
  expect_true(all(ft$prob > 0 & ft$prob < 1))
  expect_true(all(sweep(ft$rt, 2, f$par[11:15]) > 0))      # expected RT exceeds Ter
  expect_equal(nrow(predict(f)), 5)
})

test_that("data from simulate() recover the fitted parameters", {
  skip_on_cran()
  s <- sim_data(400, 8, "D", seed = 72)
  f <- quiet_fit(s$rt, s$x)
  d <- simulate(f, nsim = 1, seed = 9)[[1]]
  g <- quiet_fit(d$rt, d$x)
  expect_lt(max(abs(g$par[1:8] / f$par[1:8] - 1)), .15)
  expect_gt(cor(g$par[9:16], f$par[9:16]), .95)
})
