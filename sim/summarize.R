#!/usr/bin/env Rscript
# Summarize sim/recovery.R output. Usage: Rscript sim/summarize.R sim/recovery_results.rds
z <- readRDS(commandArgs(TRUE)[1]); res <- z$res; nit <- 20
ok <- Filter(function(r) is.na(r$error), res)
cat("Fits:", length(res), " errors:", length(res) - length(ok), "\n")
for (r in Filter(function(r) !is.na(r$error), res)) cat("  error:", r$model, r$N, r$rep, r$error, "\n")
grp <- list(a = 1:nit, v = nit + 1:nit, ter = 2*nit + 1:nit, omega_gamma = 3*nit + 1, omega_theta = 3*nit + 2)
rows <- list()
for (m in c("D", "Q")) for (N in sort(unique(sapply(ok, `[[`, "N")))) {
  rr <- Filter(function(r) r$model == m && r$N == N, ok); if (!length(rr)) next
  E <- sapply(rr, `[[`, "est"); S <- sapply(rr, `[[`, "se"); T <- rr[[1]]$truth
  S[S == -999] <- NA
  for (g in names(grp)) { k <- grp[[g]]
    e <- E[k, , drop = FALSE]; s <- S[k, , drop = FALSE]; tr <- T[k]
    bias <- rowMeans(e - tr); rmse <- sqrt(rowMeans((e - tr)^2)); sdest <- apply(e, 1, sd)
    cover <- rowMeans(abs(e - tr) <= 1.96 * s, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(model = m, N = N, param = g,
      mean_bias = mean(bias), mean_abs_bias = mean(abs(bias)), mean_rel_bias = mean(abs(bias / tr)), mean_rmse = mean(rmse),
      se_over_sd = mean(rowMeans(s, na.rm = TRUE) / sdest), coverage95 = mean(cover),
      min_coverage = min(cover), nan_se = sum(is.na(s))) }
}
tab <- do.call(rbind, rows); options(width = 200); num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], signif, 3)
cat("\nStructural parameters (averaged over the items of each type)\n"); print(tab, row.names = FALSE)

cat("\nFit diagnostics\n")
for (m in c("D", "Q")) for (N in sort(unique(sapply(ok, `[[`, "N")))) {
  rr <- Filter(function(r) r$model == m && r$N == N, ok); if (!length(rr)) next
  cat(sprintf("  %s N=%4d: %3d fits, restarted %2d, any NaN SE %2d, subjects with intractable LL %d, median time %.1f s (1 thread)\n",
      m, N, length(rr), sum(sapply(rr, `[[`, "restarted")), sum(sapply(rr, function(r) any(is.na(r$se)))),
      sum(sapply(rr, `[[`, "nr_fail")), median(sapply(rr, `[[`, "time"))))
}

cat("\nFactor scores (factest, EAP-style MAP with SEs)\n")
for (m in c("D", "Q")) for (N in sort(unique(sapply(ok, `[[`, "N")))) {
  rr <- Filter(function(r) r$model == m && r$N == N && !is.null(r$fs), ok); if (!length(rr)) next
  f <- function(j) { c(cor = mean(sapply(rr, function(r) cor(r$fs$est[, j], r$fs$true[, j]))),
    rmse = mean(sapply(rr, function(r) sqrt(mean((r$fs$est[, j] - r$fs$true[, j])^2)))),
    cover = mean(sapply(rr, function(r) mean(abs(r$fs$est[, j] - r$fs$true[, j]) <= 1.96 * r$fs$se[, j], na.rm = TRUE)))) }
  g <- f(1); h <- f(2)
  cat(sprintf("  %s N=%4d: gamma cor %.3f rmse %.3f cover %.3f | theta cor %.3f rmse %.3f cover %.3f\n", m, N, g[1], g[2], g[3], h[1], h[2], h[3]))
}
