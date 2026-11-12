test_that("LLdiff reproduces the original diffIRT C code (per subject)", {
  fx <- readRDS(test_path("fixtures", "orig_ll.rds"))
  for (nm in names(fx)) {
    f <- fx[[nm]]; N <- nrow(f$rt); J <- ncol(f$rt)
    r <- .C("LLdiff", as.double(f$rt), as.double(f$x), as.integer(N), as.integer(J), as.double(f$ai), as.double(f$vi),
            as.double(f$sd2A), as.double(f$sd2V), as.double(f$ter), as.integer(f$model), as.double(f$A), as.double(f$W),
            as.integer(f$nq), as.double(f$eps), as.double(-1), double(N), PACKAGE = "diffIRT")
    expect_equal(r[[16]], f$LL, tolerance = 1e-10, info = nm)
    expect_equal(r[[15]], f$ou, tolerance = 1e-10, info = nm)
  }
})

test_that("density is 0, quickly, when (y - Ter)/a^2 is 0 or NaN (former ARM64 hang)", {
  f <- readRDS(test_path("fixtures", "orig_ll.rds"))$D7; N <- nrow(f$rt); J <- ncol(f$rt)
  run <- function(ai, ter) .C("LLdiff", as.double(f$rt), as.double(f$x), as.integer(N), as.integer(J), as.double(ai),
     as.double(f$vi), as.double(f$sd2A), as.double(f$sd2V), as.double(ter), as.integer(1), as.double(f$A), as.double(f$W),
     as.integer(f$nq), as.double(f$eps), as.double(-1), double(N), PACKAGE = "diffIRT")
  t <- system.time({
    r1 <- run(replace(f$ai, 1, -800), f$ter)               # a overflows -> yy = 0
    r2 <- run(f$ai, replace(f$ter, 1, log(100)))           # Ter beyond every RT -> yy < 0
  })[["elapsed"]]
  expect_lt(t, 5)
  expect_true(is.finite(r1[[15]])); expect_true(is.finite(r2[[15]]))
})

test_that("results do not depend on the number of OpenMP threads", {
  skip_on_cran()
  rscript <- file.path(R.home("bin"), "Rscript")
  script <- tempfile(fileext = ".R"); out <- c(tempfile(), tempfile())
  writeLines(c(
    "suppressMessages(library(diffIRT))",
    sprintf("f <- readRDS('%s')$Q7", normalizePath(test_path("fixtures", "orig_ll.rds"))),
    "N <- nrow(f$rt); J <- ncol(f$rt); P <- 3*J+2",
    "ll <- .C('LLdiff', f$rt, f$x, N, J, f$ai, f$vi, f$sd2A, f$sd2V, f$ter, f$model, f$A, f$W, as.integer(f$nq), 1e-6, -1, double(N))[[16]]",
    "g <- .C('DERdiffAn', f$rt, f$x, N, J, 1:P, as.integer(P), f$ai, f$vi, f$sd2A, f$sd2V, f$ter, f$model, f$A, f$W, as.integer(f$nq), 1e-6, 1e-7, -1, double(N), double(P))[[20]]",
    "saveRDS(list(ll, g), commandArgs(TRUE)[1])"), script)
  env <- paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
  system2(rscript, c(script, out[1]), env = c(env, "OMP_NUM_THREADS=1"))
  system2(rscript, c(script, out[2]), env = c(env, "OMP_NUM_THREADS=4"))
  expect_identical(readRDS(out[1]), readRDS(out[2]))
})
