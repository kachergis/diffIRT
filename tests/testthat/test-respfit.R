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
