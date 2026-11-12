threads_used <- function() .C("diffIRT_threads_used", integer(1), PACKAGE = "diffIRT")[[1]]
max_threads <- function() .C("diffIRT_threads", 0L, integer(1), PACKAGE = "diffIRT")[[2]]

test_that("setThreads sets and returns the OpenMP thread count", {
  old <- diffIRT:::setThreads(NULL)          # NULL: query only
  expect_equal(max_threads(), old)
  if (old < 2) skip("compiled without OpenMP, or only one core")
  prev <- diffIRT:::setThreads(2)
  expect_equal(prev, old)
  expect_equal(threads_used(), 2)
  diffIRT:::setThreads(old)
  expect_equal(threads_used(), old)
  expect_error(diffIRT:::setThreads(1.5), "threads")
  expect_error(diffIRT:::setThreads(-1), "threads")
  expect_error(diffIRT:::setThreads(c(1, 2)), "threads")
})

test_that("control$threads is used inside diffIRT() and factest() and restored on exit", {
  skip_on_cran()
  old <- max_threads()
  if (old < 2) skip("compiled without OpenMP, or only one core")
  s <- sim_data(80, 6, "D", seed = 61)
  rec <- new.env(); rec$seen <- integer(0)
  trace("likelihood", bquote(assign("seen", c(.(rec)$seen, .C("diffIRT_threads_used", integer(1), PACKAGE = "diffIRT")[[1]]), envir = .(rec))),
        where = asNamespace("diffIRT"), print = FALSE)
  on.exit(untrace("likelihood", where = asNamespace("diffIRT")), add = TRUE)
  f1 <- quiet_fit(s$rt, s$x, se = TRUE, control = list(threads = 1))
  expect_true(length(rec$seen) > 0)
  expect_true(all(rec$seen == 1))
  expect_equal(max_threads(), old)            # restored

  rec$seen <- integer(0)
  f2 <- quiet_fit(s$rt, s$x, se = TRUE, control = list(threads = 2))
  expect_true(length(rec$seen) > 0)
  expect_true(all(rec$seen == 2))
  expect_equal(max_threads(), old)

  # results do not depend on the thread count
  expect_identical(f1$par.log, f2$par.log)
  expect_identical(f1$totLL, f2$totLL)
  expect_identical(f1$hessian, f2$hessian)

  # factest inherits the fit's setting, can override it, and restores the default
  expect_equal(f1$control$threads, 1)
  fs1 <- factest(f1, se = TRUE)
  fs2 <- factest(f1, se = TRUE, control = list(threads = 3))
  expect_identical(fs1, fs2)
  expect_equal(max_threads(), old)
})

test_that("an invalid control$threads gives an error before fitting", {
  s <- sim_data(30, 4, "D", seed = 62)
  expect_error(suppressWarnings(capture.output(diffIRT(s$rt, s$x, control = list(threads = "a")))), "threads")
})
