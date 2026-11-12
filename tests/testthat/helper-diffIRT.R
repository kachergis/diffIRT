# Shared helpers for the diffIRT tests.

quad <- function(nq) {
  Q <- statmod::gauss.quad(nq, kind = "hermite")
  list(A = Q$nodes * sqrt(2), W = Q$weights / sqrt(pi), nq = nq)
}

# Simulate data; optionally make a fraction of responses missing.
sim_data <- function(N, J, model = "D", seed = 1, miss = 0) {
  set.seed(seed)
  s <- simdiff(N, J, model = model)
  if (miss > 0) {
    k <- sample(N * J, round(miss * N * J))
    s$rt[k] <- NA; s$x[k] <- NA
  }
  s
}

# Internal-scale parameter vector (log a, v or log v, log Ter, log omega^2) at the true values.
true_par <- function(s, model = "D", sdA = .3, sdV = if (model == "D") .5 else .3) {
  c(log(s$ai), if (model == "D") s$vi else log(s$vi), log(s$ter), log(sdA^2), log(sdV^2))
}

# Objective (-2LL) and gradient as functions of the free parameters, via the package's
# own wrappers (so constraints are handled exactly as in diffIRT()).
objective <- function(s, model = "D", nq = 7, eps = 1e-6, constrain = NULL, start = NULL, aq = NULL) {
  J <- ncol(s$rt); N <- nrow(s$rt)
  rt <- s$rt; x <- s$x
  rt[is.na(rt) | is.na(x)] <- -999; x[is.na(x)] <- -999
  if (is.null(constrain)) constrain <- 1:(3 * J + 2)
  if (is.null(start)) start <- rep(NA, 3 * J + 2)
  q <- quad(nq)
  args <- list(rt = rt, x = x, N = N, nit = J, model = if (model == "D") 1 else 2, A = q$A, W = q$W,
               nq = nq, eps = eps, delta = 1e-7, start = start, constrain = constrain)
  list(L = function(p, aq_ = aq) do.call(diffIRT:::likelihood, c(list(p), args, list(aq = aq_))),
       G = function(p, aq_ = aq, grad = "analytic") do.call(diffIRT:::deriv, c(list(p), args, list(aq = aq_, grad = grad))),
       moments = function(p, aq_ = NULL) do.call(diffIRT:::aqMoments,
                   c(list(p), args[c("rt", "x", "N", "nit", "model", "A", "W", "nq", "eps", "start", "constrain")], list(aq = aq_))),
       args = args)
}

# Central differences; returns, per parameter, the smallest relative error over several step sizes
# (the truncated series makes -2LL jump very slightly, which can spoil any single step size).
grad_rel_error <- function(L, g, p, hs = c(1e-3, 1e-4, 1e-5)) {
  errs <- sapply(hs, function(h) {
    cd <- sapply(seq_along(p), function(j) { a <- p; a[j] <- a[j] + h; b <- p; b[j] <- b[j] - h; (L(a) - L(b)) / (2 * h) })
    abs(g - cd) / pmax(1, abs(cd))
  })
  apply(matrix(errs, length(p)), 1, min)
}

quiet_fit <- function(...) { out <- utils::capture.output(f <- diffIRT(...)); f }
