#!/usr/bin/env Rscript
# Calibration and power of RTFit() and PersonFit() (issue #10).
# 100 replications per condition, N = 500, 10 items, 10% MCAR missing, D and Q models.
# Null: data from the model. Item misfit: item 1's RT delayed by 0.3 s for a random 30% of
# persons. Person misfit (5% of persons, on items 6-10): inattentive (RT doubled, random
# responses) or rapid guessing (decision time U(0, .05), random responses).
# Usage: Rscript sim/fit_statistics.R
suppressMessages({ library(diffIRT); library(parallel) })
N <- 500; J <- 10; g <- 1:25; half <- 6:10
gen <- function(model, k, cond) {
  set.seed(3000 + k)
  s <- simdiff(N, J, model = model, miss = .1)
  if (cond == "item") { i <- sample(N, .3 * N); s$rt[i, 1] <- s$rt[i, 1] + .3 }
  if (cond == "inattentive") { s$rt[g, half] <- 2 * s$rt[g, half]; s$x[g, half] <- rbinom(25 * 5, 1, .5) }
  if (cond == "guessing") { s$rt[g, half] <- matrix(runif(25 * 5, 0, .05), 25) + rep(s$ter[half], each = 25)
                            s$x[g, half] <- rbinom(25 * 5, 1, .5) }
  s$rt[is.na(s$x)] <- NA
  s
}
one <- function(model, k, cond) {
  s <- gen(model, k, cond)
  capture.output(f <- diffIRT(s$rt, s$x, model = model, se = TRUE, control = list(threads = 1)))
  if (f$conv != 0) return(NULL)
  r <- RTFit(f)
  pf <- PersonFit(f, nsim = 500, seed = k)$persons$p
  list(item.p = r$items$p, overall.p = r$overall[["p"]], person.p = pf)
}
for (model in c("D", "Q")) for (cond in c("null", "item", "inattentive", "guessing")) {
  res <- mclapply(1:100, function(k) tryCatch(one(model, k, cond), error = function(e) conditionMessage(e)), mc.cores = 12)
  bad <- !sapply(res, is.list)
  if (any(bad)) { cat("errors:\n"); print(table(unlist(res[bad]))) }
  res <- res[!bad]
  ip <- do.call(rbind, lapply(res, `[[`, "item.p")); op <- sapply(res, `[[`, "overall.p")
  pp <- do.call(rbind, lapply(res, `[[`, "person.p"))
  cat(sprintf("%s %-11s (%d reps)  RTFit rejected at .05: item 1 %.2f, items 2-%d %.3f, overall %.2f | ",
              model, cond, length(res), mean(ip[, 1] < .05), J, mean(ip[, -1] < .05), mean(op < .05)))
  cat(sprintf("PersonFit flagged: aberrant %.2f, others %.3f | KS p (items 2-%d, null p uniform) %.2f\n",
              mean(pp[, g] < .05), mean(pp[, -g] < .05), J, suppressWarnings(ks.test(as.vector(ip[, -1]), "punif")$p.value)))
}
