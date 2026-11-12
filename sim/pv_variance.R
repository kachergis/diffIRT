#!/usr/bin/env Rscript
# Do plausible values reproduce the population variances? (issue #8)
# 100 replications per model, N = 300, 10 items, 10% MCAR missing; true omega from the
# simdiff() defaults (log gamma: .3; D theta: .5; Q log theta: .3). Compares the variance
# across persons of MAP, EAP and plausible-value scores (on the latent scale) with omega^2.
# Usage: Rscript sim/pv_variance.R
suppressMessages({ library(diffIRT); library(parallel) })
one <- function(model, k) {
  set.seed(2000 + k)
  s <- simdiff(300, 10, model = model, miss = .1)
  capture.output(f <- diffIRT(s$rt, s$x, model = model, se = TRUE, control = list(threads = 1)))
  if (f$conv != 0) return(NULL)
  lt <- function(v) cbind(log(v[, 1]), if (model == "Q") log(v[, 2]) else v[, 2])
  map <- lt(factest(f)); eap <- lt(factest(f, type = "EAP"))
  pv  <- suppressWarnings(factest(f, type = "PV", npv = 5))
  pvv <- rowMeans(sapply(pv, function(v) apply(lt(v), 2, var)))
  c(map = apply(map, 2, var), eap = apply(eap, 2, var), pv = pvv, est = exp(f$par.log[31:32]))
}
for (m in c("D", "Q")) {
  res <- mclapply(1:100, function(k) tryCatch(one(m, k), error = function(e) conditionMessage(e)), mc.cores = 12)
  bad <- !sapply(res, is.numeric)
  if (any(bad)) { cat("errors:\n"); print(table(unlist(res[bad]))) }
  r <- do.call(rbind, res[!bad])
  if (anyNA(r)) cat("replications with NA:", which(apply(is.na(r), 1, any)), "\n")
  true <- c(.3^2, if (m == "D") .5^2 else .3^2)
  cat(sprintf("%s (%d reps)  true omega^2: %.4f %.4f\n", m, nrow(r), true[1], true[2]))
  for (w in c("map", "eap", "pv", "est")) {
    cols <- grep(paste0("^", w), colnames(r))
    cat(sprintf("  %-4s mean var: gamma %.4f (se %.4f)  theta %.4f (se %.4f)\n", w,
                mean(r[, cols[1]], na.rm = TRUE), sd(r[, cols[1]], na.rm = TRUE) / sqrt(nrow(r)), mean(r[, cols[2]], na.rm = TRUE), sd(r[, cols[2]], na.rm = TRUE) / sqrt(nrow(r))))
  }
}
