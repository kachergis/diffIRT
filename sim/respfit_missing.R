#!/usr/bin/env Rscript
# Calibration of RespFit Mr under MCAR missingness (issue #11): 100 replications,
# N = 500, 8 items, 10% and 20% missing, new pairwise moments vs the old
# "missing = incorrect" moments. Usage: Rscript sim/respfit_missing.R
suppressMessages({library(diffIRT); library(parallel)})
oldMr <- function(f, order = 2) {   # previous behaviour: missing counted as incorrect
  x <- f$score; x[is.na(x)] <- 0; g <- f; g$score <- x; suppressWarnings(RespFit(g, order))$Mr }
for (m in c("D", "Q")) for (miss in c(.1, .2)) {
  res <- mclapply(1:100, function(k) {
    set.seed(1000 + k); s <- simdiff(500, 8, model = m)
    i <- sample(length(s$x), round(miss * length(s$x))); s$x[i] <- NA; s$rt[i] <- NA
    capture.output(f <- diffIRT(s$rt, s$x, model = m, control = list(threads = 1)))
    if (f$conv != 0) return(c(NA, NA, NA))
    r <- RespFit(f, 2); c(r$Mr, r$df, oldMr(f)) }, mc.cores = 12)
  z <- do.call(rbind, res); z <- z[complete.cases(z), , drop = FALSE]; df <- z[1, 2]
  cat(sprintf("%s miss=%.1f  n=%d df=%d | new: mean %.1f var %.1f rej %.2f KS p %.2f | old: mean %.1f rej %.2f\n",
    m, miss, nrow(z), df, mean(z[,1]), var(z[,1]), mean(z[,1] > qchisq(.95, df)), ks.test(z[,1], "pchisq", df)$p.value,
    mean(z[,3]), mean(z[,3] > qchisq(.95, df))))
}
