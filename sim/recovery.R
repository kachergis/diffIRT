#!/usr/bin/env Rscript
# Parameter recovery for diffIRT (installed version).
# Usage: [NQ=7] [ADAPT=1] Rscript sim/recovery.R <reps> <out.rds> [N values...]
# Item parameters are drawn once per model and held fixed; person parameters
# are redrawn every replication from the simdiff() defaults:
#   log gamma_p ~ N(0, .3^2); D: theta_p ~ N(0, .5^2); Q: log theta_p ~ N(0, .3^2).
Sys.setenv(OMP_NUM_THREADS = "1")          # parallelise over replications instead
suppressMessages({ library(diffIRT); library(parallel) })
args <- commandArgs(TRUE)
reps <- as.integer(args[1]); out <- args[2]
Ns   <- if (length(args) > 2) as.integer(args[-(1:2)]) else c(250, 500, 1000)
nit  <- 20
nq   <- as.integer(Sys.getenv("NQ", "7"))
adapt <- Sys.getenv("ADAPT", "0") == "1"

set.seed(2026)
items <- list(
  D = list(ai = exp(runif(nit, -1, 0)), vi = runif(nit, -1.5, .5), ter = runif(nit, .5, 1.5), sdA = .3, sdV = .5),
  Q = list(ai = exp(runif(nit, -1, 0)), vi = exp(runif(nit, 0, 1)), ter = runif(nit, .5, 1.5), sdA = .3, sdV = .3))

one <- function(model, N, rep) {
  set.seed(1e5 * (model == "Q") + 1e3 * N + rep)
  it <- items[[model]]
  gam <- rlnorm(N, 0, it$sdA)
  the <- if (model == "D") rnorm(N, 0, it$sdV) else rlnorm(N, 0, it$sdV)
  s <- simdiff(N, nit, ai = it$ai, vi = it$vi, gamma = gam, theta = the, ter = it$ter, model = model)
  truth <- c(it$ai, it$vi, it$ter, it$sdA, it$sdV)
  msg <- character(0)
  t <- system.time(f <- tryCatch(withCallingHandlers(
         { o <- capture.output(r <- diffIRT(s$rt, s$x, model = model, se = TRUE, control = list(nq = nq, adaptive = adapt))); msg <- o; r },
         warning = function(w) invokeRestart("muffleWarning")),
       error = function(e) conditionMessage(e)))[["elapsed"]]
  if (is.character(f)) return(list(model = model, N = N, rep = rep, error = f, time = t))
  fs <- tryCatch(suppressWarnings(factest(f, se = TRUE)), error = function(e) NULL)
  list(model = model, N = N, rep = rep, error = NA, time = t,
       truth = truth, est = f$par, se = f$std.err, LL = f$totLL,
       restarted = any(grepl("Restarted", msg)), nr_fail = f$nr_fail,
       fs = if (is.null(fs)) NULL else list(est = fs[, c(1, 3)], se = fs[, c(2, 4)], true = cbind(gam, the)))
}

grid <- expand.grid(rep = seq_len(reps), N = Ns, model = c("D", "Q"), stringsAsFactors = FALSE)
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(k) one(grid$model[k], grid$N[k], grid$rep[k]),
                mc.cores = max(1, detectCores() - 2), mc.preschedule = FALSE)
saveRDS(list(items = items, res = res, grid = grid, nq = nq, adaptive = adapt), out)
cat("done:", nrow(grid), "fits in", format(round(Sys.time() - t0, 1)), "\n")
